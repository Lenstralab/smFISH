#!/usr/local/bin/ipython3 -i

#########################################################################################################
#########################################################################################################
import os
import sys
import re
import numpy as np
import scipy
import pandas
from scipy import optimize
from scipy import ndimage
import copy as copy2
import yaml
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
from lmfit.models import GaussianModel
import scipy.stats
from scipy.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.morphology import binary_erosion
from datetime import datetime
import math
import rtree
from parfor import parfor, parpool
from tiffwrite import tiffwrite, IJTiffFile
from tllab_common.wimread import imread as imr
from tllab_common.transforms import Transform
from tllab_common.findcells import findcells
from tllab_common.misc import save_roi, ipy_debug
import psutil
from time import time
from tqdm.auto import tqdm
from PIL import Image, ImageDraw, ImageFont
from numbers import Number
from collections.abc import Iterable

from warnings import filterwarnings, warn
filterwarnings('ignore', 'loadtxt: Empty input file:')

if __package__ is None or __package__ == '':
    import spotDetection_functions as sdf
    import misc
    import _version
else:
    from . import spotDetection_functions as sdf
    from . import misc
    from . import _version


#########################################################################################################
#########################################################################################################
if not '__file__' in locals(): #when executed using execfile
    import inspect
    __file__ = inspect.getframeinfo(inspect.currentframe()).filename


#########################################################################################################
### Defining functions to make nuclear, cellular and TS histograms
#########################################################################################################
def histNuc(distr, bins, freqAxis, xCenter, exclBinsCount, color, outName):
    hist_count = np.histogram(distr, bins=bins - .5, density=True)[0]
    hist_count_norm = hist_count[exclBinsCount:] / sum(hist_count[exclBinsCount:])
    bins = bins[exclBinsCount:-1]
    fig = plt.figure()
    plt.bar(bins, hist_count_norm)
    if freqAxis != "None": plt.ylim([0, freqAxis])
    plt.title("Number of RNA in the nucleus, " + color)
    plt.xlabel('Nr RNA/nucleus')
    plt.ylabel('Frequency')
    plt.text(xCenter, max(hist_count_norm),
             "Mean nuclear count = " + str('{0:.2f}'.format((np.mean(distr)))) + " +/- " + str(
                 '{0:.2f}'.format((np.std(distr) / np.sqrt(len(distr))))) + "\nNumber of cells analyzed = " + str(len(distr)),
             fontsize=12, horizontalalignment='center', verticalalignment='top')
    plt.savefig(outName, format='pdf')
    plt.close(fig)

def histCell(distr, bins, freqAxis, xCenter, exclBinsCount, color, outName):
    hist_count = np.histogram(distr, bins=bins - .5, density=True)[0]
    hist_count_norm = hist_count[exclBinsCount:] / sum(hist_count[exclBinsCount:])
    bins = bins[exclBinsCount:-1]
    fig = plt.figure()
    plt.bar(bins, hist_count_norm)
    if freqAxis != "None": plt.ylim([0, freqAxis])
    plt.title("Number of RNA in the cell, " + color)
    plt.xlabel('Nr RNA/cell')
    plt.ylabel('Frequency')
    ##
    plt.text(xCenter, max(hist_count_norm),
             "Mean cellular count = " + str('{0:.2f}'.format((np.mean(distr)))) + " +/- " + str(
                 '{0:.2f}'.format((np.std(distr) / np.sqrt(len(distr))))) + "\nNumber of cells analyzed = " + str(len(distr)),
             fontsize=12, horizontalalignment='center', verticalalignment='top')
    plt.savefig(outName, format='pdf')
    plt.close(fig)

def histTSvals(distr, bins, exclBins):
    h = np.histogram(distr, bins=bins - .5, density=True)[0]
    normalize = lambda x: x / sum(x[exclBins:])
    hNorm = h / sum(h[exclBins:])
    return hNorm


def FitBurstModel(histTSvals, bins, exclBins):
    P = lambda x: x[0] ** bins / (x[0] + 1) ** (bins + x[1]) * scipy.special.binom(x[1] + bins - 1, x[1] - 1)
    Pnorm = lambda x: (lambda p: p / sum(p[exclBins:]))(P(x))
    [n, aT] = optimize.fmin(lambda x: sum((Pnorm(x) - histTSvals)[exclBins:] ** 2), [4, 160. / 70])
    Pres = P([n, aT])
    Pres *= sum(histTSvals[exclBins:]) / sum(Pres[exclBins:])
    return (Pres, n, aT)

def FitPoissonModel(histTSvals, bins, exclBins):
    PPoiss = lambda x: x ** bins / scipy.special.gamma(1 + bins)
    PPoissnorm = lambda x: (lambda p: p / sum(p[exclBins:]))(PPoiss(x))
    cT = optimize.fmin(lambda x: sum((PPoissnorm(x) - histTSvals)[exclBins:] ** 2), [.1])
    PPoissres = PPoiss(cT[0])
    PPoissres *= sum(histTSvals[exclBins:]) / sum(PPoissres[exclBins:])
    return (PPoissres, cT)

def histTS(distr, hist, bins, Pres, n, aT, PPoissres, cT, freqAxis, xCenter, color, params, outName):
    fig = plt.figure()
    plt.bar(bins, hist, label="FISH")
    if params['FitBurst'] == 1:
        plt.plot(Pres, color='red',
                 label="Bursting model n = " + str(round(n, 3)) + " RNAs/burst\naT = " + str(round(aT, 3)))
    if params['FitPoisson'] == 1:
        plt.plot(PPoissres, color='green', label="Poisson model Mean = " + str(round(cT[0], 3)))
    plt.legend(loc='center right')
    ymax = 1.1 * max(hist)
    if freqAxis != misc.none():
        ymax = freqAxis
    plt.ylim([0, ymax])
    plt.title("Number of nascent RNA at the TS, " + color)
    plt.xlabel('Nr nascent RNA at TS')
    plt.ylabel('Frequency')
    plt.text(xCenter, 0.98 * ymax, "Mean nr nascent RNA = " + str('{0:.2f}'.format((np.mean(distr)))) + " +/- " + str(
        '{0:.2f}'.format((np.std(distr) / np.sqrt(len(distr))))) + "\nNote: Does not include the " + str(
        range(0, params['exclBins'])) + " bins", fontsize=12, horizontalalignment='center', verticalalignment='top')
    plt.savefig(outName, format='pdf')
    plt.close(fig)

def writestatsfile(outFilePart, color, distrCell, distrNuc, distrTS, params):
    statsfile = open(outFilePart + "_statistics_" + color + ".txt", "w+")
    statsfile.write("Number of cells analyzed: " + str(len(distrCell)))
    statsfile.write("\nCellular count " + color + ": " + str('{0:.2f}'.format(np.mean(distrCell))) + " +/- " + str(
        '{0:.2f}'.format(np.std(distrCell) / np.sqrt(len(distrCell)))))
    statsfile.write("\nNuclear count " + color + ": " + str('{0:.2f}'.format((np.mean(distrNuc)))) + " +/- " + str(
        '{0:.2f}'.format(np.std(distrNuc) / np.sqrt(len(distrNuc)))))
    statsfile.write("\nNote:These counts do not include the " + str(range(0, params['exclBinsCount'])) + " bins")
    statsfile.write("\nTxpnsite average " + color + ": " + str('{0:.2f}'.format(np.mean(distrTS))) + " +/- " + str(
        '{0:.2f}'.format(np.std(distrTS) / np.sqrt(len(distrTS)))))
    statsfile.write("\nNote:The TS average does not include the " + str(range(0, params['exclBins'])) + " bins")
    if params['includeMultTS'] == 1:
        statsfile.write("\nThe multiple txpnsite data was calculated with TS threshold " + str(params['thresholdTS']))
    elif params['includeMultTS'] == 2:
        statsfile.write("\nThe multiple txpnsite data was calculated with " + str(params['nrTS2Include']) + " TSs per cell")
    statsfile.write(
        "\nNote: errors reported above are standard errors assuming a normal distribution, but the distribution is not normal")
    statsfile.close()


#########################################################################################################
#########################################################################################################

def calculate_general_parameters(params):
    params['fileListIn'] = [os.path.join(params['folderIn'], file) for file in params['fileListIn']]

    if params['useControlForEmptyCells'] == 1 and params['channelsToAnalyze'] == [0, 1]:
        if params['controlColor'] == 1:
            params['channelsToAnalyze'] = [1, 0]

    os.makedirs(params['analysisfolder'], exist_ok=True)
    dim = params['localize']['dimension']
    if isinstance(dim, str):
        params['localize']['dimension'] = int(dim[0])
    return params


def get_exp_info(pathIn, params):
    paramsexp = {}

    #### make file list and output directory
    paramsexp['date'] = re.findall('(?:\d[_-]?){7}\d', pathIn)[0]
    paramsexp['expName'] = os.path.basename(pathIn.rstrip(os.path.sep))
    paramsexp['pathOut'] = os.path.join(params['outputfolder'], f"{paramsexp['date']}_{paramsexp['expName']}")
    paramsexp['outfilePart1'] = os.path.join(paramsexp['pathOut'],  f"{paramsexp['date']}_{paramsexp['expName']}")
    os.makedirs(paramsexp['pathOut'], exist_ok=True)
    paramsexp['zSlices'] = None if params['zSlices'] is None else \
        range(params['zSlices']) if isinstance(params['zSlices'], Number) else range(*params['zSlices'])

    return paramsexp


def get_filelist(pathIn, params):
    lfn = os.listdir(pathIn)
    if "display_and_comments.txt" in lfn: lfn.remove("display_and_comments.txt")
    if "._display_and_comments.txt" in lfn: lfn.remove("._display_and_comments.txt")
    if ".DS_Store" in lfn: lfn.remove(".DS_Store")
    if "._.DS_Store" in lfn: lfn.remove("._.DS_Store")
    lfn.sort()
    if "excluded" in lfn: lfn.remove("excluded")
    if params['processAll'] == 0:
        lfn = lfn[0:params['nrFilesToAnalyze']]
    return lfn

def get_turret_filters(pathIn, lfn, paramsexp):
    allFiles = os.listdir(os.path.join(pathIn, lfn[0]))
    if "metadata.txt" in allFiles: allFiles.remove("metadata.txt")
    for element in allFiles:
        if element[0] == '.': allFiles.remove(element)

    splitFiles = [x.split("_") for x in allFiles]
    Turretfilter = []
    for i in range(len(splitFiles)):
        if len(splitFiles[i]) > 4:
            Turretfilter.append("_".join(splitFiles[i][2:-1]))
        else:
            Turretfilter.extend(splitFiles[i][2:-1])
    paramsexp['ChNames'] = sorted(set(Turretfilter))
    return paramsexp


###########################################################################################################
##### Does max projection of images
###########################################################################################################

def max_projection(pathIn, lfn, paramsexp):
    # channels are saved in alphabetical order, which probably is different from the order in the original data
    print("Making maximum intensity projections")
    path = os.path.join(paramsexp['pathOut'], paramsexp['date'])
    for fn in lfn:
        outfile = path + '_{}_max.tif'.format(os.path.splitext(fn)[0])
        if not os.path.exists(outfile):
            with imr(os.path.join(pathIn, fn)) as im:
                with IJTiffFile(outfile, (im.shape[2], 1, im.shape[4]), dtype='uint16', pxsize=im.pxsize,
                                deltaz=im.deltaz) as tif:
                    for c, C in enumerate(sorted(range(im.shape[2]), key=lambda i: im.cnamelist[i])):
                        for t in range(im.shape[4]):
                            tif.save(im.max(C, paramsexp['zSlices'], t), c, 0, t)


############################################################################################################
##### Runs cell profiler to get cell and nucleus masks
############################################################################################################

def run_cellprofiler(lfn, params, paramsexp):
    print("Making nucleus and cell masks")
    shutil.copyfile(params['PipelineCellProfiler'], paramsexp['pathOut'] + params['PipelineCellProfiler'])
    for filenr in range(1, (len(lfn) + 1)):
        if not os.path.exists(get_tif(paramsexp['pathOut'] + paramsexp['date'] + "_" + lfn[(filenr - 1)] + "_max_cell_mask.tif")):
            os.system("cellprofiler -p " + re.escape(params['PipelineCellProfiler']) + " -c -i " + re.escape(
                paramsexp['pathOut']) + " -o " + re.escape(paramsexp['pathOut']) + " -f " + str(filenr) + " -l " + str(filenr))


############################################################################################################
#### Runs findCells ipython script to get cell and nucleus masks
############################################################################################################

def find_cells(lfn, params, paramsexp):
    path = os.path.join(paramsexp['pathOut'], paramsexp['date'])
    if params['channelCellMask'] == "Cy3":
        channelCell = 0
    elif params['channelCellMask'] == "Cy5":
        channelCell = 1
    else:
        channelCell = 0
        print("params['channelCellMask'] is not defined correctly, taking channel 0")
    if params['channelsToAnalyze'] == [0]:
        channelNuc = 1
    elif params['channelsToAnalyze'] == [1]:
        channelNuc = 2
    elif params['channelsToAnalyze'] == [0, 1] or params['channelsToAnalyze'] == [1, 0]:
        channelNuc = 2
    else:
        raise Exception("no valid channels defined")
    if params['threshold'] == 'None': threshold = None
    else: threshold = params['threshold']
    if params['thresholdnuc'] =="None": thresholdnuc = None
    else: thresholdnuc = params['thresholdnuc']
    for fn in tqdm(lfn, desc="Making nucleus and cell masks"):
        if not get_tif(path + "_{}_max_cell_mask.tif".format(fn)) or not get_tif(path + "_{}_max_nucleus_mask.tif"):
            with imr(path + '_{}_max.tif'.format(fn)) as im:
                cells, nuclei = findcells(im(channelCell), im(channelNuc),
                    cellcolormask=None, ccdist=params['ccdist'], threshold=threshold, thresholdnuc=thresholdnuc,
                    removeborders=params['removeborders'], thres=params['thres'], thresnuc=params['thresnuc'],
                    smooth=params['smooth'], smoothnuc=params['smoothnuc'], minfeatsize=params['minfeatsize'],
                    minfeatsizenuc=params['minfeatsizenuc'], dilate=params['dilate'], dilatenuc=params['dilatenuc'])
            tiffwrite(path + '_{}_max_nucleus_mask.tif'.format(fn), nuclei.astype('uint16'), colormap='glasbey',
                      pxsize=im.pxsize, deltaz=im.deltaz)
            tiffwrite(path + '_{}_max_cell_mask.tif'.format(fn), cells.astype('uint16'), colormap='glasbey',
                      pxsize=im.pxsize, deltaz=im.deltaz)

############################################################################################################
#### runs threshold optimization
############################################################################################################

def run_optimize_threshold(pathIn, lfn, params, paramsexp, colors=None):
    if colors is None:
        colors = 'cy3', 'cy5'
    for channel in params['channelsToAnalyze']:
        color = colors[channel]
        path = os.path.join(paramsexp['pathOut'], paramsexp['date'])

        spotslist = spotslistlocalmax = np.zeros(50)
        if params.get('focianalysis'):  # focianalysis
            if params['localize']['threshold']['method'].upper() == 'SD':
                tbins = 1
                threshvals = np.arange(1., 51., tbins)  # multiplied by SD later
            else:
                tbins = 25
                threshvals = range(100, 1350, tbins)
        else:
            tbins = 2
            threshvals = list(range(0, 40, tbins)) + list(range(40, 120, tbins * 2)) + list(range(120, 200, tbins * 4))
        locMaxValAll = []

        for fn in tqdm(lfn, desc="Optimizing Threshold " + color, disable=len(lfn)<5):
            if params['localize']['dimension'] == 2:
                with imr(path + "_{}_max.tif".format(fn), dtype=float) as im:
                    imBpass = sdf.bpass(im(channel), *(np.array(params['localize']['bpassLimits']) *
                                                       params['localize']['psfPx']))
            elif params['localize']['dimension'] == 3:
                with imr(os.path.join(pathIn, fn), dtype=float) as im:
                    imchannel = sorted(range(im.shape[2]), key=lambda i: im.cnamelist[i])[channel]
                    imBpass = sdf.bpass3D(im[imchannel, paramsexp['zSlices']].transpose((3, 0, 1, 2, 4)).squeeze(),
                                          *(np.array(params['localize']['bpassLimits']) * params['localize']['psfPx']),
                                          *(np.array(params['localize']['bpassLimits']) * params['localize']['psfPxZ']))
            else:
                raise ValueError("parameter 'localize.dimension' not defined correctly")

            locMax = ndimage.maximum_filter(imBpass, 3)
            locMax = (locMax == imBpass) * 1  # finds were in bandpass image the maximum value was, rest will be 0
            locMax[0, ...], locMax[-1, ...], locMax[:, 0, ...], locMax[:, -1, ...] = 0, 0, 0, 0
            if locMax.ndim == 3:
                locMax[:, :, 0], locMax[:, :, -1] = 0, 0

            locMaxCoo = tuple(np.array(np.where(locMax)))  # finds coordinates
            locMaxVal = imBpass[locMaxCoo]
            locMaxValAll.extend(locMaxVal)
            if params.get('focianalysis') and params['localize']['threshold']['method'].upper() == 'SD':
                f = np.std(imBpass)
            else:
                f = 1
            # this is the most time intensive step
            spotslist += [ndimage.label(imBpass > thresh * f)[1] for thresh in threshvals]
            spotslistlocalmax += [np.sum(locMaxVal > thresh * f) for thresh in threshvals]

        locMaxValAll = np.asarray(locMaxValAll)
        thresh = get_threshold(params, color, imBpass)

        if params.get('focianalysis'):  # focianalysis
            file = os.path.splitext(lfn[0])[0]
            loc_max_pdf = f"{path}_{file}_pixel_values_local_maxima_{color}_" \
                          f"{params['localize']['threshold']['method']}.pdf"
            thresh_opt_pdf = f"{path}_{file}_threshold_optimization_{color}_" \
                             f"{params['localize']['threshold']['method']}.pdf"
            thresh_txt = f"{path}_{file}_threshold_optimization_{color}_" \
                         f"{params['localize']['threshold']['method']}.txt"
        else:
            loc_max_pdf = f"{path}_{paramsexp['expName']}_threshold_optimization_{color}_" \
                          f"{params['localize']['dimension']}D.pdf"
            thresh_opt_pdf = f"{path}_{paramsexp['expName']}_pixel_values_local_maxima_{color}_" \
                             f"{params['localize']['dimension']}D.pdf"
            thresh_txt = f"{path}_{paramsexp['expName']}_thresholds_{color}.txt"

        fig = plt.figure()
        plt.plot(threshvals, spotslist)
        plt.plot(threshvals, spotslistlocalmax)
        if params.get('focianalysis') and params['localize']['threshold']['method'].upper() == 'SD':
            plt.xlabel(f'Threshold SD {color}')
        else:
            plt.xlabel(f'Threshold {color}')
        plt.ylabel('Number of spots detected')
        plt.yscale('log')
        if thresh is not None:
            plt.axvline(thresh, c='r', ls='-', label='Your threshold')
        plt.legend(loc='upper center')
        plt.savefig(thresh_opt_pdf, format='pdf')
        plt.close(fig)

        np.savetxt(thresh_txt, spotslist, delimiter="\t")

        fig = plt.figure(figsize=(8, 3))
        h, x = np.histogram(np.log10(np.maximum(1e-12, locMaxValAll)), bins=np.r_[-1:np.log10(locMaxValAll.max()):.01])
        plt.xscale('log')
        plt.xlabel('Pixel intensity')
        plt.yscale('log')
        plt.ylabel('Count')
        plt.axvline(np.median(locMaxValAll), ls='--', label='median')
        plt.axvline(np.mean(locMaxValAll), ls='-', label='mean')
        plt.axvline(np.mean(locMaxValAll) + np.var(locMaxValAll) ** .5, ls=':', label='mean + sd')
        if thresh is not None:
            plt.axvline(thresh, c='r', ls='-', label='Your threshold')
        plt.text(1, max(h), "Mean + SD = " + str(
            '{0:.2f}'.format(np.mean(locMaxValAll) + np.var(locMaxValAll) ** .5)) + "\nMean = " + str(
            '{0:.2f}'.format(np.mean(locMaxValAll))) + "\nSD = " + str(
            '{0:.2f}'.format(np.var(locMaxValAll) ** .5)) + "\nYour threshold = " + str(thresh), fontsize=12,
                 horizontalalignment='center', verticalalignment='top')
        plt.plot(10 ** x[:-1], h + 1)
        plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.))
        plt.tight_layout()
        plt.savefig(loc_max_pdf, format='pdf')
        plt.close(fig)

############################################################################################################
#### runs localize
############################################################################################################


def rm_duplicates_idxs(points, r):
    """ Remove duplicates """
    if len(points) == 0:
        return np.array([])
    elif len(points) == 1:
        return np.array((0,))

    def dist(a, b):
        return math.hypot(*[n - m for n, m in zip(a, b)])
    keep = []
    prop = rtree.index.Property()
    prop.dimension = points.shape[1]
    index = rtree.index.Index(properties=prop)
    for i, p in enumerate(points):
        nearby = index.intersection([q - r for q in p] + [q + r for q in p])
        if all(dist(p, points[j]) >= r for j in nearby):
            keep.append(i)
            index.insert(i, (*p, *p))
    return np.array(keep)


def get_threshold(params, color, imBpass):
    if params['localize']['threshold']['method'].upper() == 'SD':
        return params['localize']['threshold'][color] * np.std(imBpass)
    else:
        return params['localize']['threshold'][color]


def filter_localizations(res, filter, pxsize, psfPx, shape, deltaz=1, psfPxZ=1, zSlices=None, verbose=False):
    if not 'z' in res:
        res['z'], res['dz'], res['Oz'], res['dOz'], res['z_ini'] = 0, 0, 0, 0, 0
    res['dx_um'] = res['dx'] * pxsize
    res['dy_um'] = res['dy'] * pxsize
    res['dz_um'] = res['dz'] * deltaz
    res['dr_um'] = np.sqrt(res['dx_um'] ** 2 + res['dy_um'] ** 2 + res['dz_um'] ** 2)

    if l := res.shape[0]:
        q = []
        if len(shape) == 2:
            q.append(f'-0.5<=y<={shape[0] - 0.5} & -0.5<=x<={shape[1] - 0.5}')
        elif len(shape) == 3:
            q.append(f'-0.5<=y<={shape[1] - 0.5} & -0.5<=x<={shape[2] - 0.5}')
            if zSlices is None:
                q.append(f'-0.5<=z<={shape[0] - 0.5}')
            elif isinstance(zSlices, Iterable):
                q.append(f'{max(0, zSlices[0]) - 0.5}<=z<={min(zSlices[1] + 1, shape[0]) - 0.5}')
            else:
                q.append(f'-0.5<=z<={min(zSlices + 1, shape[0]) - 0.5}')
        if (i := filter.get('distance')) != misc.none():
            q.append(f'((x-x_ini)/{psfPx}) ** 2 + ((y-y_ini)/{psfPx}) ** 2 + ((z-z_ini)/{psfPxZ}) ** 2 < {i ** 2}')
        if (i := filter.get('intensity')) != misc.none():
            q.append(f'Ii>{i}')
        if (i := filter.get('dintensity')) != misc.none():
            q.append(f'dIi/abs(Ii)<{i}')
        if (i := filter.get('peak_intensity')) != misc.none():
            q.append(f'Ip>{i}')
        if (i := filter.get('dpeak_intensity')) != misc.none():
            q.append(f'dIp/abs(Ip)<{i}')
        if (i := filter.get('precision')) != misc.none():
            q.append(f'dr_um<{i}')
        if (i := filter.get('offset')) != misc.none():
            q.append(f'O>{i}')
        if (i := filter.get('doffset')) != misc.none():
            q.append(f'dO/abs(O)<{i}')
        if (i := filter.get('tilt')) != misc.none():
            q.append(f'Ox<{i}&Oy<{i}&Oz<{i}')
        if (i := filter.get('dtilt')) != misc.none():
            q.append(f'dOx/abs(Ox)<{i}&dOy/abs(Oy)<{i}&dOz/abs(Oz)<{i}')
        if (i := filter.get('R2')) != misc.none():
            q.append(f'R2>{i}')
        if (i := filter.get('R2_peak')) != misc.none():
            q.append(f'R2_peak>{i}')
        if (i := filter.get('R2_background')) != misc.none():
            q.append(f'R2_background>{i}')
        if (i := filter.get('s')) != misc.none():
            q.append(f'{i[0]}<s<{i[1]}')
        if (i := filter.get('ds')) != misc.none():
            q.append(f'ds<{i}')
        if (i := filter.get('sz')) != misc.none():
            q.append(f'{i[0]}<sz<{i[1]}')
        if (i := filter.get('dsz')) != misc.none():
            q.append(f'dsz<{i}')
        if q:
            res = res.query('&'.join(q)).copy()
        if (s := filter.get('separation')) != misc.none():
            idxs = rm_duplicates_idxs(res[['x', 'y', 'z']].to_numpy(float) / (psfPx, psfPx, psfPxZ), s)
            res = res.iloc[idxs].copy()
        if verbose:
            print(f'Filter is removing {l - res.shape[0]} of {l} spots.')
    return res


def add_cell_as_column(df, path, fn):
    file = os.path.splitext(fn)[0]
    nuc_file = get_tif(f"{path}_{file}_max_nucleus_mask.tif")
    cell_file = get_tif(f"{path}_{file}_max_cell_mask.tif")
    if nuc_file:
        with imr(nuc_file, dtype=int) as im:
            nuc_mask = im(0)
        try:
            df['nucleus'] = [nuc_mask[int(np.round(row['y'])), int(np.round(row['x']))] for _, row in df.iterrows()]
        except IndexError:
            df['nucleus'] = 0
    else:
        warn(f'Nucleus mask for {fn} not found, adding nucleus number to results skipped.')
    if cell_file:
        with imr(cell_file, dtype=int) as im:
            cell_mask = im(0)
        try:
            df['cell'] = [cell_mask[int(np.round(row['y'])), int(np.round(row['x']))] for _, row in df.iterrows()]
        except IndexError:
            df['cell'] = 0
    else:
        warn(f'Cell mask for {fn} not found, adding cell number to results skipped')

    if cell_file and nuc_file:
        # Warn for inconsistencies between cell and nucleus masks
        c, n = set(cell_mask.flatten()), set(nuc_mask.flatten())
        if not c == n:
            if c - n:
                warn(f'Cell(s) {c - n} in {path}_{file} do(es) not have a nucleus.')
            if n - c:
                warn(f'Nucle(us/i) {n - c} in {path}_{file} do(es) not have a cell.')
    return df


def run_localize(pathIn, lfn, params, paramsexp, channel_names=None):
    if channel_names is None:
        channel_names = ('cy3', 'cy5')

    def output_files(path, fn, color):
        file = os.path.splitext(fn)[0]
        if params.get('focianalysis'):  # focianalysis
            fnTxt = f"{path}_{file}_loc_results_{color}_{params['localize']['threshold']['method']}.txt"
            fnTif = f"{path}_{file}_loc_results_{color}_{params['localize']['threshold']['method']}.tif"
            fnRoi = f"{path}_{file}_loc_results_{color}_{params['localize']['threshold']['method']}.roi"
        else:
            fnTxt = f'{path}_{file}_loc_results_{color}.txt'
            fnTif = f'{path}_{file}_loc_results_{color}.tif'
            fnRoi = f'{path}_{file}_loc_results_{color}.roi'
        return fnTxt, fnTif, fnRoi

    print(f"Running Localize using {params['localize']['dimension']}D")

    par = yaml.dump(params['localize'])
    for channel in params['channelsToAnalyze']:
        color = channel_names[channel]
        path = os.path.join(paramsexp['pathOut'], paramsexp['date'])

        if params['localize']['dimension'] == 2:  ### localize in 2D
            for fn in lfn:
                fnTxt, fnTif, fnRoi = output_files(path, fn, color)
                if not os.path.exists(fnTxt):
                    # open image
                    with imr(os.path.join(pathIn, fn), dtype=float) as im:
                        pxsize = im.pxsize
                    # save 6 frame image: original image, imBpass, imBinary (areas above threshold), local maxima,
                    # maxima above threshold and fitted spots. These are useful for optimizing the threshold.
                    with IJTiffFile(fnTif, (6, 1, 1), dtype='uint16', pxsize=im.pxsize) as tif:
                        # TODO: check transforms, also fit first, then transform
                        with imr(f'{path}_{fn}_max.tif',
                                 transform=os.path.exists(os.path.join(pathIn, "beads.czi")), dtype=float) as im:
                            shape = im.shape[:2]
                            # bandpass image to find spots
                            imBpass = sdf.bpass(im(channel), *(np.array(params['localize']['bpassLimits']) *
                                                               params['localize']['psfPx']))
                            tif.save(im(channel).astype('uint16'), 0, 0, 0)
                            tif.save(imBpass.astype('uint16'), 1, 0, 0)
                            thresh = get_threshold(params, color, imBpass)

                            # find local maxima
                            locMax = ndimage.maximum_filter(imBpass, 3)
                            # finds were in bandpass image the maximum value was, rest will be 0
                            locMax = (locMax == imBpass) * 1
                            locMax[0, :], locMax[-1, :], locMax[:, 0], locMax[:, -1] = 0, 0, 0, 0
                            tif.save(locMax, 3, 0, 0)
                            locMaxCoo = tuple(np.array(np.where(locMax)))  # finds coordinates
                            locMaxVal = imBpass[locMaxCoo]
                            yG, xG, valG = tuple(np.array(locMaxCoo + (locMaxVal,))[:, np.where(locMaxVal > thresh)])
                            cooGuess = np.r_[xG, yG].T

                            # for optimizing threshold
                            imGuess = np.zeros(shape, dtype='uint16')
                            for a in cooGuess:
                                imGuess[int(a[1]), int(a[0])] = 1  #
                            tif.save(imGuess, 4, 0, 0)
                            tif.save((imBpass > thresh), 2, 0, 0)  # Binary image after fixed threshold

                            # generate an ROI and its coordinates on demand
                            def ROIgen(im, cooGuess, winSize):
                                for cg in cooGuess:
                                    window, roi = sdf.crop[im, cg[1] - winSize:cg[1] + winSize,
                                                           cg[0] - winSize:cg[0] + winSize]
                                    yield roi.astype(float), window[::-1, 0], cg - window[::-1, 0]

                            @parfor(ROIgen(im(channel), cooGuess, params['localize']['winSize']), (params,),
                                    length=len(cooGuess), nP=8, bar=False)
                            def fun(ROIgen, params):
                                psfPx = params['localize']['psfPx']
                                distance = params['localize']['filter'].get('distance')
                                max_distance = None if distance is None else distance / psfPx
                                fit_params = params['localize']['fit']
                                fit0 = sdf.GaussianMaskFit(*ROIgen, psfPx, params['localize']['winSize'],
                                                           max_distance=max_distance, optimize_width=False,
                                                           correct_fit_window=params['localize']['correct_fit_window'],
                                                           refine_fit=params['localize']['refine_fit'],
                                                           max_sigma=fit_params.get('max_sigma'),
                                                           fit_weight_sigma=fit_params.get('fit_weight_sigma'))
                                if params['localize']['optimize_width']:
                                    fit1 = sdf.GaussianMaskFit(*ROIgen, psfPx, params['localize']['winSize'],
                                                               max_distance=max_distance, optimize_width=True,
                                                               reset_if_not_fit=False)
                                    return fit1.data if fit1.fit else fit0.data
                                else:
                                    return fit0.data

                        # concatenate results. If no spots are detected, use empty dataframe.
                        fitResults = pandas.concat(fun, ignore_index=True) if fun else sdf.GaussianMaskFit.empty_data(2)
                        fitResults = add_cell_as_column(fitResults, path, fn)
                        fitResults.to_csv(fnTxt.replace('.txt', '_unfiltered.txt'), sep='\t', index=False)
                        fitResults = filter_localizations(fitResults, params['localize']['filter'], pxsize,
                                                          params['localize']['psfPx'], shape, verbose=True)

                        # Save the results in a text file
                        fitResults.to_csv(fnTxt, sep='\t', index=False)

                        # Save the results in a roi file
                        roi_coo_guess = pandas.DataFrame(cooGuess, columns=('x', 'y'))
                        roi_coo_guess['C'] = 4
                        roi_fit_results = fitResults[['x', 'y']].copy()
                        roi_fit_results['C'] = 5
                        save_roi(fnRoi, pandas.concat((roi_fit_results, roi_coo_guess)),
                                 (*im.shape[1:], 3, im.shape[0], 1), 'xyC__')

                        # Compute and save image results as a tif file
                        imSpots = np.zeros(shape, dtype='uint16')  # empty image
                        imSpots[round(fitResults['y']).astype(int), round(fitResults['x']).astype(int)] = 1
                        tif.save(imSpots, 5, 0, 0)

                    print(f'*** Found {len(fitResults)} spots in channel {color}. ***')
                    print(f'Results saved in {fnTxt} and {fnTif}.')
                    print('\n\n')

                    par += f'\nImage {fn} in channel {color} has threshold: {get_threshold(params, color, imBpass)}, ' \
                           f'found {len(fitResults)} spots'
                    if len(fitResults):
                        par += f", of which {sum(fitResults['fit']):.0f} " \
                               f"({100 * sum(fitResults['fit']) / len(fitResults):.2f}%) from fit."
                    else:
                        par += '.'

        else:  ### localize in 3D
            for fn in lfn:
                fnTxt, fnTif, fnRoi = output_files(path, fn, color)
                if not os.path.exists(fnTxt):
                    with imr(os.path.join(pathIn, fn), dtype=float) as jm:
                        im = jm[channel, :, 0].squeeze((2, 4)).transpose(2, 0, 1)
                        z_slices = jm.get_czt(0, paramsexp['zSlices'], 0)[1]
                        pxsize, deltaz = jm.pxsize, jm.deltaz
                    with IJTiffFile(fnTif, (3, len(z_slices), 1), dtype='uint16', pxsize=pxsize, deltaz=deltaz) as tif:
                        for z, zim in enumerate(im[z_slices]):
                            tif.save(zim, 0, z, 0)
                        imBpass = sdf.bpass3D(im,
                                        *(np.array(params['localize']['bpassLimits']) * params['localize']['psfPx']),
                                        *(np.array(params['localize']['bpassLimits']) * params['localize']['psfPxZ']))

                        # find local maxima
                        locMax = ndimage.maximum_filter(imBpass, 3)
                        # finds were in bandpass image the maximum value was, rest will be 0
                        locMax = (locMax == imBpass) * 1
                        locMax[0, :, :], locMax[-1, :, :] = 0, 0
                        locMax[:, 0, :], locMax[:, -1, :] = 0, 0
                        locMax[:, :, 0], locMax[:, :, -1] = 0, 0
                        locMaxCoo = tuple(np.array(np.where(locMax)))  # finds coordinates
                        locMaxVal = imBpass[locMaxCoo]
                        thresh = get_threshold(params, color, imBpass)

                        zG, yG, xG, valG = tuple(np.array(locMaxCoo + (locMaxVal,))[:, np.where(locMaxVal > thresh)])
                        cooGuess = np.r_[zG, yG, xG].T  # for 2D x and y should be swapped!

                        zDepth = (params['localize']['winSizeZ'] - 1) // 2

                        intersect = [value for value in zG[0] if value in set(z_slices[:zDepth] + z_slices[-zDepth:])]

                        # generate an ROI and its coordinates on demand
                        def ROIgen(im, cooGuess, winSize, winSizeZ):
                            for cg in cooGuess:
                                window, roi = sdf.crop[im, cg[0] - winSizeZ:cg[0] + winSizeZ,
                                                       cg[1] - winSize:cg[1] + winSize, cg[2] - winSize:cg[2] + winSize]
                                yield roi.astype(float), window[:, 0], cg - window[:, 0]

                        # fit each spot with 3D gaussian with tilted plane.
                        # The GaussianMaskFit3D function is described in spotDetection_Functions.py.
                        # Write spots to parameter fitResults.
                        @parfor(ROIgen(im, cooGuess, params['localize']['winSize'], params['localize']['winSizeZ']),
                                (params,), length=len(cooGuess), nP=8, bar=False)
                        def fun(ROIgen, params):
                            psfPx, psfPxZ = params['localize']['psfPx'], params['localize']['psfPxZ']
                            distance = params['localize']['filter'].get('distance')
                            if distance is not None:
                                max_distance, max_distance_z = distance / psfPx, distance / psfPxZ
                            else:
                                max_distance, max_distance_z = None, None
                            fit_params = params['localize']['fit']
                            fit0 = sdf.GaussianMaskFit(*ROIgen, psfPx, params['localize']['winSize'],
                                                       psfPxZ, params['localize']['winSizeZ'],
                                                       max_distance=max_distance, max_distance_z=max_distance_z,
                                                       optimize_width=False,
                                                       correct_fit_window=params['localize']['correct_fit_window'],
                                                       refine_fit=params['localize']['refine_fit'],
                                                       max_sigma=fit_params.get('max_sigma'),
                                                       max_sigma_z=fit_params.get('max_sigma_z'),
                                                       fit_weight_sigma=fit_params.get('fit_weight_sigma'))
                            if params['localize']['optimize_width']:
                                fit1 = sdf.GaussianMaskFit(*ROIgen, psfPx, params['localize']['winSize'],
                                                           psfPxZ, params['localize']['winSizeZ'],
                                                           max_distance=max_distance, max_distance_z=max_distance_z,
                                                           optimize_width=True, reset_if_not_fit=False)
                                return fit1.data if fit1.fit else fit0.data
                            else:
                                return fit0.data

                        # concatenate results. If no spots are detected, use empty dataframe.
                        fitResults = pandas.concat(fun, ignore_index=True) if fun else sdf.GaussianMaskFit.empty_data(3)
                        fitResults = add_cell_as_column(fitResults, path, fn)
                        fitResults.to_csv(fnTxt.replace('.txt', '_unfiltered.txt'), sep='\t', index=False)
                        fitResults = filter_localizations(fitResults, params['localize']['filter'], pxsize,
                                                          params['localize']['psfPx'], im.shape, deltaz,
                                                          params['localize']['psfPxZ'], (z_slices[0], z_slices[-1]),
                                                          verbose=True)

                        # Save the results in a text file
                        fitResults.to_csv(fnTxt, sep='\t', index=False)

                        # Save the results in a roi file
                        roi_coo_guess = pandas.DataFrame(cooGuess, columns=('z', 'y', 'x'))
                        roi_coo_guess['C'] = 1
                        roi_fit_results = fitResults[['x', 'y', 'z']].copy()
                        roi_fit_results['C'] = 2
                        save_roi(fnRoi, pandas.concat((roi_fit_results, roi_coo_guess)),
                                 (*im.shape[1:], 3, im.shape[0], 1), 'xyCz_')

                        # for optimizing threshold
                        imGuess = np.zeros(im.shape, dtype='uint16')
                        for a in cooGuess:
                            imGuess[int(a[0]), int(a[1]), int(a[2])] = 1  #
                        for z, zim in enumerate(imGuess[z_slices]):
                            tif.save(zim, 1, z, 0)

                        # Compute and save image results as a tif file
                        imSpots = np.zeros(im.shape, dtype='uint16')  # empty image
                        imSpots[round(fitResults['z']).astype(int), round(fitResults['y']).astype(int),
                                round(fitResults['x']).astype(int)] = 1
                        for z, zim in enumerate(imSpots[z_slices]):
                            tif.save(zim, 2, z, 0)

                    print(f'*** Found {len(fitResults)} spots in channel {color}. ***')
                    print(f'Results saved in {fnTxt} and {fnTif}.')
                    if intersect:
                        print(f'Warning: {len(intersect)} out of {len(zG[0])} spots were at outer frames of z focus '
                              'and may not have been fit!')
                    print('\n\n')

                    par += f'\nImage {fn} in channel {color} has threshold: {get_threshold(params, color, imBpass)}, ' \
                           f'found {len(fitResults)} spots'
                    if len(fitResults):
                        par += f", of which {sum(fitResults['fit']):.0f} "\
                               f"({100 * sum(fitResults['fit']) / len(fitResults):.2f}%) from fit. "\
                               f"{len(intersect)} spots where at outer frames of z focus."
                    else:
                        par += '.'

        # save parameters localize
        parameterfile = f"{os.path.join(paramsexp['pathOut'], paramsexp['date'])}_localize_parameters.txt"
        with open(parameterfile, 'w') as file:
            file.write(par)


def get_tif(file):  # tif or tiff
    name, ext = os.path.splitext(file)
    if not ext in ('.tif', '.tiff'):
        name += ext
    for ext in ('.tif', '.tiff'):
        if os.path.exists(name + ext):
            return name + ext
    else:
        return None


############################################################################################################
#### makes outlines on cell and nucleus and boxes around localized spots
############################################################################################################

def make_mask_image(lfn, params, paramsexp, colors=None):
    print("Combining mask with image")
    if colors is None:
        colors = 'cy3', 'cy5'
    for fn in lfn:
        file = os.path.splitext(fn)[0]
        path = os.path.join(paramsexp['pathOut'], paramsexp['date'])
        if not os.path.exists(f"{path}_{file}_mask_cell+nuc+spots.tif"):
            maximagefile = f"{path}_{file}_max.tif"
            if params.get('focianalysis'):  # focianalysis
                locfile = [f"{path}_{file}_loc_results_{colors[c]}_{params['localize']['threshold']['method']}.txt"
                           for c in params['channelsToAnalyze']]
            else:
                locfile = [f"{path}_{file}_loc_results_{colors[c]}.txt" for c in params['channelsToAnalyze']]
            nucleusfile = get_tif(f"{path}_{file}_max_nucleus_mask.tif")
            cellfile = get_tif(f"{path}_{file}_max_cell_mask.tif")

            outfileMask = f"{path}_{file}_mask_cell+nuc+spots.tif"
            outfileRoi = f"{path}_{file}_mask_cell+nuc+spots.roi"
            with imr(maximagefile, dtype='uint16') as im, \
                imr(nucleusfile, dtype='uint16') as imn, imr(cellfile, dtype='uint16') as imc, \
                IJTiffFile(outfileMask, (im.shape[2] + len(locfile) + 2, 1, 1), dtype='uint16', pxsize=im.pxsize,
                           deltaz=im.deltaz) as tif:
                for c in range(im.shape[2]):  # copy frames
                    tif.save(im(c), c, 0, 0)

                for i, lbl in enumerate((imc(0), imn(0))):  # draw cell and nucleus outlines
                    mask = np.zeros(im.shape[:2])
                    lbl_n = np.unique(lbl)
                    for cell in lbl_n[lbl_n > 0]:
                        mask += (lbl == cell).astype(int) - binary_erosion(lbl == cell, np.ones((3, 3)))
                    tif.save(mask, im.shape[2] + i, 0, 0)

                datas = []
                for c, l in enumerate(locfile, im.shape[2] + 2):  # draw dots for localisations
                    data = pandas.read_table(l)
                    Pos = np.zeros(im.shape[:2])
                    Pos[round(data['y']).astype(int), round(data['x']).astype(int)] = 1
                    tif.save(Pos, c, 0, 0)
                    data['C'] = c
                    if not data.empty:
                        datas.append(data)
                if datas:
                    save_roi(outfileRoi, pandas.concat(datas), im.shape[:2] + tif.shape, 'xyC__')


############################################################################################################
#### makes montage of cell mask, nuclear mask and max projection
############################################################################################################

def make_montage(path, lfn, outfile, suffix, scale=4, colormap=None):
    def resize(im, scale):
        T = Transform()
        T.origin = [0, 0]
        T.parameters = (scale, 0, 0, scale, 0, 0)
        return T.frame(im)[:im.shape[0]//scale, :im.shape[1]//scale]

    outfile += f'_{suffix}.tif'
    posY = [int(fn.split("Pos")[-1].split("_")[1]) for fn in lfn]
    posX = [int(fn.split("Pos")[-1].split("_")[2]) for fn in lfn]
    xMax = max(posX)
    xMin = min(posX)
    yMax = max(posY)
    yMin = min(posY)

    if not os.path.exists(outfile):
        im = [imr(path + f'_{fn}_{suffix}.tif', dtype='int16') for fn in lfn]
        shapeIm = im[0].shape[:2]
        nchannels = im[0].shape[2]
        colors = ('red', 'lime', 'blue', 'magenta', 'cyan', 'yellow')[:nchannels] + ('white',)
        shape = ((yMax - yMin + 1) * shapeIm[0], (xMax - xMin + 1) * shapeIm[1])
        kwargs = {'colors': colors} if colormap is None else {'colormap': colormap}
        with IJTiffFile(outfile, (nchannels + (colormap is None), 1, 1), 'uint16', pxsize=im[0].pxsize/scale,
                        deltaz=im[0].deltaz, **kwargs) as tif:
            for c in range(nchannels):
                montage = np.zeros(shape, dtype='uint16')
                for i, (fm, x, y) in enumerate(zip(im, posX, posY)):
                    jm = np.flip(np.rot90(fm(c), k=-1, axes=(0, 1)), axis=1)
                    montage[(y - yMin) * shapeIm[0]:(y - yMin + 1) * shapeIm[0],
                            (x - xMin) * shapeIm[1]:(x - xMin + 1) * shapeIm[1]] = jm
                if not scale == 1:
                    montage = resize(montage, scale)
                if colormap is not None:
                    montage[::shapeIm[0] // scale, :] = 10
                    montage[:, ::shapeIm[1] // scale] = 10
                tif.save(montage, c, 0, 0)
            if colormap is None:
                text = Image.new('1', montage.shape[::-1])
                draw = ImageDraw.Draw(text)
                for x, y in zip(posX, posY):
                    draw.text(((x - xMin + 0.1) * shapeIm[1] / scale, (y - yMin + 0.1) * shapeIm[0] / scale),
                       f'{y}_{x}', 1, ImageFont.truetype(os.path.join(os.path.dirname(__file__), 'arial.ttf'), size=20))
                text = np.array(text)
                text[::shapeIm[0] // scale, :] = 1
                text[:, ::shapeIm[1] // scale] = 1
                plt.imshow(text)
                tif.save((65535 * text).astype('uint16'), nchannels, 0, 0)
        for fm in im:
            fm.close()


############################################################################################################
#### calculate transcription sites
############################################################################################################
def CellCycle(lfn, params, paramsexp):
    path = os.path.join(paramsexp['pathOut'], paramsexp['date'])

    @parfor(lfn, (path,), desc='Calculating cell intensities', terminator=imr.kill_vm)
    def fun(fn, path):
        maxFile = path + "_{}_max.tif".format(fn)
        nucleusFile = get_tif(path + "_{}_max_nucleus_mask.tif".format(fn))
        with imr(maxFile) as im:
            img_dapi = im(-1)
        with imr(nucleusFile, dtype='uint16') as im:
            mask_nuc = im(0)
        nuclei = list(np.unique(mask_nuc))
        if np.min(nuclei) == 0: nuclei.remove(0)
        if (len(nuclei) !=0 and np.max(nuclei) == 65535): nuclei.remove(65535)
        return [np.sum((mask_nuc == n) * img_dapi) for n in nuclei]

    integrated_intensity_cell = np.array(sum(fun, []))  # list of integrated nuclear intensities of all cells
    integrated_intensity_cell = integrated_intensity_cell[misc.outliers(integrated_intensity_cell)]

    print("Plotting initial nuclear intensities histogram")
    np.savetxt(paramsexp['outfilePart1'] + "_dapi_intensity_histogram.txt", integrated_intensity_cell)
    plt.figure()
    plt.hist(integrated_intensity_cell, bins=50, density=True, color='dimgray')
    plt.xlabel('Nucleus Intensity')
    plt.ylabel('Frequency')
    plt.title('Histogram of Nuclear DAPI Intensities')
    plt.savefig(paramsexp['outfilePart1'] + "_dapi_intensity_histogram.pdf", format='pdf')

    if params['fitNucHistogram'] == 1:
        print("Fitting nuclear intensities histogram")
        if params['removeOutliers'] != 'None':
            integrated_intensity_cell = integrated_intensity_cell[integrated_intensity_cell < params['removeOutliers']]
        maxDAPI = integrated_intensity_cell.max()
        minDAPI = integrated_intensity_cell.min()
        weights, means, sds = misc.gaussian_unmix(integrated_intensity_cell, 2)
        idx = np.argsort(means)
        weights, means, sds = weights[idx], means[idx], sds[idx]
        meanG1cells, meanG2cells = means
        stdevG1cells, stdevG2cells = sds

        sigmaG1cells = params['sigmaG1cells']
        sigmaG2cells = params['sigmaG2cells']
        G1cellswithbuds = [minDAPI, meanG1cells+(sigmaG1cells[1]*stdevG1cells)]
        G1cells = [meanG1cells-(sigmaG1cells[0]*stdevG1cells), meanG1cells+(sigmaG1cells[1]*stdevG1cells)]
        Scells = [meanG1cells+(sigmaG1cells[1]*stdevG1cells)+0.1, meanG2cells-(sigmaG2cells[0]*stdevG2cells)-0.1]
        G2cells = [meanG2cells-(sigmaG2cells[0]*stdevG2cells), meanG2cells+(sigmaG2cells[1]*stdevG2cells)]
        with open(paramsexp['outfilePart1'] + "_dapi_thresholds_nrcells.txt", "w+") as threshfile:
            threshfile.write("Thresholds: ")
            threshfile.write("\nG1cellswithbuds: " + str(G1cellswithbuds))
            threshfile.write("\nG1cells: " + str(G1cells))
            threshfile.write("\nScells: " + str(Scells))
            threshfile.write("\nG2cells: " + str(G2cells))
            threshfile.write("\nminDAPI: " + str(minDAPI))
            threshfile.write("\nmaxDAPI: " + str(maxDAPI))
            threshfile.write("\n")
            threshfile.write("\nNumber of Cells: ")
            threshfile.write("\nG1cellswithbuds: " + str(len([i for i in integrated_intensity_cell if G1cellswithbuds[0] <= i < G1cellswithbuds[1]])))
            threshfile.write("\nG1cells: " + str(len([i for i in integrated_intensity_cell if G1cells[0] <= i < G1cells[1]])))
            threshfile.write("\nScells: " + str(len([i for i in integrated_intensity_cell if Scells[0] <= i < Scells[1]])))
            threshfile.write("\nG2cells: " + str(len([i for i in integrated_intensity_cell if G2cells[0] <= i < G2cells[1]])))
        np.save(paramsexp['outfilePart1'] + "_dapi_thresholds.npy", [G1cellswithbuds, G1cells, Scells, G2cells])

        print("Plotting final nuclear intensity histogram")
        fig2 = plt.figure()
        ax = fig2.add_subplot(111)

        x = np.linspace(minDAPI, maxDAPI, 1000)
        g = np.sum([misc.gaussian(x, A, mu, s) for A, mu, s in zip(weights, means, sds)], 0)
        plt.hist(integrated_intensity_cell, bins=50, density=True, color='dimgray')
        plt.plot(x, g, 'k')
        plt.axvline(G1cells[0], ls='--', color='b', label='G1cells')
        plt.axvline(G1cells[1], ls='--', color='b')
        plt.axvline(G2cells[0], ls='-', color='darkgreen', label='G2cells')
        plt.axvline(G2cells[1], ls='-', color='darkgreen')
        plt.legend()
        plt.text(1, 0.5, 'z-score: \nG1cells = {}\nG2cells = {}'.format(sigmaG1cells, sigmaG2cells),
                 horizontalalignment='right', verticalalignment='top', transform = ax.transAxes)
        plt.xlabel('Nucleus Intensity')
        plt.ylabel('Frequency')
        plt.title('Histogram of Nuclear Intensities')
        plt.xlim(0, maxDAPI)
        plt.savefig(paramsexp['outfilePart1'] + "_dapi_intensity_histogram_final.pdf", format='pdf')

    if params['makeMasks'] == 1:
        if not os.path.exists(path + "_{}_max_cell_mask_threshold.tif".format(lfn[-1]))\
                or not os.path.exists(path + "_{}_max_nucleus_mask_threshold.tif".format(lfn[-1])):

            CellCycleStages = params['CellCycleStages']
            if params['manualThreshold'] in (None, 'None'):
                for c in CellCycleStages:
                    if not os.path.exists(os.path.join(paramsexp['pathOut'], c)):
                        os.makedirs(os.path.join(paramsexp['pathOut'], c))
                thresholds = np.load(paramsexp['outfilePart1'] + "_dapi_thresholds.npy")
                order = ['G1cellswithbuds', 'G1cells', 'Scells', 'G2cells']
            else:
                for c in CellCycleStages:
                    if not os.path.exists(os.path.join(paramsexp['pathOut'], 'ManualThreshold', c)):
                        os.makedirs(os.path.join(paramsexp['pathOut'], 'ManualThreshold', c))
                thresholds = params['manualThreshold']
                if not isinstance(thresholds[0], list):
                    thresholds = [thresholds]
                order = CellCycleStages

            print("Saving nucleus and cell masks of {} classified cells".format(CellCycleStages))

            def fun(fn, get_tif, path, CellCycleStages, thresholds, order):
                maxFile = path + "_{}_max.tif".format(fn)
                nucleusFile = get_tif(path + "_{}_max_nucleus_mask.tif".format(fn))
                cellFile = get_tif(path + "_{}_max_cell_mask.tif".format(fn))
                with imr(maxFile) as im:
                    img_dapi = im(-1)
                    pxsize, deltaz, timeinterval = im.pxsize, im.deltaz, im.timeinterval
                with imr(nucleusFile, dtype='uint16') as im:
                    mask_nuc = im(0)
                with imr(cellFile, dtype='uint16') as im:
                    mask_cell = im(0)
                nuclei = list(np.unique(mask_nuc))
                if np.min(nuclei) == 0: nuclei.remove(0)
                if (len(nuclei) !=0 and np.max(nuclei) == 65535): nuclei.remove(65535)
                I = [np.sum((mask_nuc == n) * img_dapi) for n in nuclei]
                res = []
                for c in CellCycleStages:
                    threshold = thresholds[order.index(c)]
                    nmask = np.zeros(mask_nuc.shape, dtype='uint16')
                    cmask = np.zeros(mask_cell.shape, dtype='uint16')
                    for n, i in zip(nuclei, I):
                        if threshold[0] <= i < threshold[1]:
                            nmask += n * (mask_nuc == n)
                            cmask += n * (mask_cell == n)
                    res.append((nmask, cmask, pxsize, deltaz, timeinterval))
                return res

            with parpool(fun, (get_tif, path, CellCycleStages, thresholds, order)) as pool:
                for fn in lfn:
                    pool[fn] = fn  # adding work to the parallel pool
                for _ in tqdm(lfn, desc='Writing cell cycle labels'):
                    fn, images = pool.get_newest()  # taking results from pool in fifo order and save as tifs
                    for c, (mask_nuc, mask_cell, pxsize, deltaz, timeinterval) in zip(CellCycleStages, images):
                        if params['manualThreshold'] in (None, 'None'):
                            p = os.path.join(paramsexp['pathOut'], c, paramsexp['date'])
                        else:
                            p = os.path.join(paramsexp['pathOut'], 'ManualThreshold', c, paramsexp['date'])
                        tiffwrite(p + "_{}_max_nucleus_mask_threshold.tif".format(fn),
                                  mask_nuc.astype('uint16'), colormap='glasbey', pxsize=pxsize, deltaz=deltaz)
                        tiffwrite(p + "_{}_max_cell_mask_threshold.tif".format(fn),
                                  mask_cell.astype('uint16'), colormap='glasbey', pxsize=pxsize, deltaz=deltaz)


############################################################################################################
#### calculate transcription sites
############################################################################################################

def get_common_path_out(params, paramsexp, CellCycleStages, c):
    if params['CellCycle'] == 0:
        return os.path.join(paramsexp['pathOut'], f"{paramsexp['date']}_{paramsexp['expName']}")
    elif params['CellCycle'] == 1 and params['manualThreshold'] is None:
        return os.path.join(paramsexp['pathOut'], CellCycleStages[c], f"{paramsexp['date']}_{paramsexp['expName']}")
    elif params['CellCycle'] == 1 and params['manualThreshold'] is not None:
        return os.path.join(paramsexp['pathOut'], 'ManualThreshold', CellCycleStages[c],
                            f"{paramsexp['date']}_{paramsexp['expName']}")


def calculate_TS(lfn, params, paramsexp):
    # TODO: use localize file and pandas to do this much simpler
    print("Calculating transcription sites")

    CellCycleStages = params['CellCycleStages']
    if params['CellCycle'] == 0: Z = 1
    elif params['CellCycle'] == 1 and params['manualThreshold'] is None: Z = len(CellCycleStages)
    elif params['CellCycle'] == 1 and params['manualThreshold'] is not None: Z = 1

    path = os.path.join(paramsexp['pathOut'], paramsexp['date'])

    for c in range(Z):
        imageTSCy3 = []
        imageTSCy5 = []

        ##### create list of dictionaries (frames) with 0 for each nucleus
        ld_zero = [{} for i in range(0, len(lfn))]
        ld_zero5 = [{} for i in range(0, len(lfn))]

        for fi, fn in enumerate(lfn):
            if params['MaskTS'] == "nucleus":
                if params['CellCycle'] == 0:
                    nucleusfile = get_tif(path + "_{}_max_nucleus_mask.tif".format(fn))
                elif params['CellCycle'] == 1 and params['manualThreshold'] is None:
                    nucleusfile = os.path.join(paramsexp['pathOut'], CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_nucleus_mask_threshold.tif".format(fn))
                elif params['CellCycle'] == 1 and params['manualThreshold'] is not None:
                    nucleusfile = os.path.join(paramsexp['pathOut'], 'ManualThreshold', CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_nucleus_mask_threshold.tif".format(fn))
            elif params['MaskTS'] == "cell":
                if params['CellCycle'] == 0:
                    nucleusfile = get_tif(path + "_{}_max_cell_mask.tif".format(fn))
                elif params['CellCycle'] == 1 and params['manualThreshold'] is None:
                    nucleusfile = os.path.join(paramsexp['pathOut'], CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_cell_mask_threshold.tif".format(fn))
                elif params['CellCycle'] == 1 and params['manualThreshold'] is not None:
                    nucleusfile = os.path.join(paramsexp['pathOut'], 'ManualThreshold', CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_nucleus_mask_threshold.tif".format(fn))
            with imr(nucleusfile) as im:
                nuclist = np.unique(im(0))
            # nucleus = pil.open(nucleusfile)
            # nuclist = np.unique(np.array(nucleus.getdata()).reshape(nucleus.size[::-1]))
            for nuc in nuclist[1:]:
                ld_zero[fi][nuc] = 0
                ld_zero5[fi][nuc] = np.array([0, 0, 0, 0, 0])

        ##### loop over frames
        for channel in params['channelsToAnalyze']:
            if channel == 0: color = "cy3"
            if channel == 1: color = "cy5"

            ldCountNuc = copy2.deepcopy(ld_zero)  # parameter file for number of spots in nucleus
            ldCountCell = copy2.deepcopy(ld_zero)  # parameter file for number of spots in cytoplasm (excluding nucleus)
            ldCountBoth = copy2.deepcopy(
                ld_zero)  # parameter file for number of spots in entire cell (both cytoplasm and nucleus)
            ld = copy2.deepcopy(ld_zero)  # intensity for brightest TS
            ldPos = copy2.deepcopy(ld_zero5)  # for brightest TS with position information
            ldNucTS = copy2.deepcopy(ld_zero)  # for all nucl spots (filter later for TS)
            ldNucTSPos = copy2.deepcopy(ld_zero5)  # for all nucl spots (filter later for TS) wiht position information
            ldCytoPos = copy2.deepcopy(
                ld_zero5)  # for all spots in  cytoplams (without nucleus) with position information

            dataNucleus = []  # parameter file for spot intensities in nucleus
            dataCell = []  # parameter file for spot intensities in cytoplasm (excluding nucleus)

            for fi, fn in enumerate(tqdm(lfn, desc='Computing TS')):
                locfile = path + "_{}_loc_results_{}.txt".format(fn, color)
                if params['CellCycle'] == 0:
                    nucleusfile = get_tif(path + "_{}_max_nucleus_mask.tif".format(fn))
                    cellfile = get_tif(path + "_{}_max_cell_mask.tif".format(fn))
                elif params['CellCycle'] == 1 and params['manualThreshold'] is None:
                    nucleusfile = os.path.join(paramsexp['pathOut'], CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_nucleus_mask_threshold.tif".format(fn))
                    cellfile = os.path.join(paramsexp['pathOut'], CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_cell_mask_threshold.tif".format(fn))
                elif params['CellCycle'] == 1 and params['manualThreshold'] is not None:
                    nucleusfile = os.path.join(paramsexp['pathOut'], 'ManualThreshold', CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_nucleus_mask_threshold.tif".format(fn))
                    cellfile = os.path.join(paramsexp['pathOut'], 'ManualThreshold', CellCycleStages[c],
                                                paramsexp['date'] + "_{}_max_cell_mask_threshold.tif".format(fn))

                ##### load data
                data = pandas.read_table(locfile)
                with imr(nucleusfile, dtype=float) as imn, imr(cellfile, dtype=float) as imc:
                    nucleusarray = imn(0)
                    cellarray = imc(0)

                # select for cell and nucleus spots
                for _, row in data.iterrows():
                    x = int(round(row['x']))
                    y = int(round(row['y']))
                    Ixyzfit = row[['I', 'x', 'y', 'z', 'fit']].to_numpy(float)
                    nucnr = nucleusarray[y, x]
                    cellnr = cellarray[y, x]
                    if nucnr > 0 and nucnr != 32768 and nucnr != 65535.0:
                        ldCountNuc[fi][nucnr] += 1  ### add count to ldCountNuc
                        ldCountBoth[fi][nucnr] += 1  ### add count to ldCountBoth
                        if row['fit'] > 0:  ### check if fit was not forced
                            dataNucleus.append(row['I'])  ### add intensity to dataNucleus
                        if params['includeMultTS'] > 0:
                            ldNucTS[fi][nucnr] = np.vstack((ldNucTS[fi][nucnr], row['I']))
                            ldNucTSPos[fi][nucnr] = np.vstack((ldNucTSPos[fi][nucnr], Ixyzfit))
                        if row['I'] > ld[fi][nucnr]:
                            ld[fi][nucnr] = row['I']  ### replace intensity in ld if spot has higher intensity
                            ldPos[fi][nucnr] = Ixyzfit  ### replace intensity and x,y location, fit
                    elif cellnr > 0 and cellnr != 32768 and nucnr != 65535.0:  ### if spot is in cytoplasm but not in nucleus
                        ldCountCell[fi][cellnr] += 1  ### add count to ldCountCell
                        ldCountBoth[fi][cellnr] += 1  ### add count to ldCountBoth
                        ldCytoPos[fi][cellnr] = np.vstack((ldCytoPos[fi][cellnr], Ixyzfit))  ### add cyto spots to ldCytoPos, also if fit is 0
                        if row['fit'] > 0:  ### check if fit is not forced
                            dataCell.append(row['I'])  ### add intensity to dataCell
                        if params['MaskTS'] == "cell":
                            if params['includeMultTS'] > 0:
                                ldNucTS[fi][cellnr] = np.vstack((ldNucTS[fi][cellnr], row[0]))  ### add intensity cell spot
                                ldNucTSPos[fi][cellnr] = np.vstack((ldNucTSPos[fi][cellnr], Ixyzfit))
                            if row['I'] > ld[fi][cellnr]:
                                ld[fi][cellnr] = row['I']  ### replace intensity in ld if spot has higher intensit
                                ldPos[fi][cellnr] = Ixyzfit  ### replace intensity and x,y location
                    else:
                        continue

            ### warning, the unfiltered data includes spots that are not fit, and includes data for cells that have no RNAs. So this unfiltered data is not using the parameters params['ExclEmpty']. This unfiltered data is useful if you want to compare Cy3 and Cy5 data from the same cell, for example during the params['CalcSingleCellCorr']. If you want filtering, this occurs in the steps below, and the results will be saved as "_CountNuc+Cytoplasm_"+color+".txt" and "_intensity_brightest_txpnsite_"+color+".txt".

            common_path = get_common_path_out(params, paramsexp, CellCycleStages, c)

            # save unfiltered data
            countUnfiltered = []
            for d in range(0, int(len(ldCountBoth))): countUnfiltered.extend(ldCountBoth[d].values())

            dataTxpnSiteUnfiltered = []
            for d in range(0, int(len(ld))): dataTxpnSiteUnfiltered.extend(ld[d].values())
            np.savetxt(common_path + "_CountNuc+Cytoplasm_" + color + "_nofilter.txt", countUnfiltered, fmt='%.8e', delimiter='\t')
            np.savetxt(common_path + "_intensity_brightest_txpnsite_" + color + "_nofilter.txt", dataTxpnSiteUnfiltered, fmt='%.8e', delimiter='\t')

            dataTxpnSiteUnfilteredPos = []
            for d in range(0, int(len(ldPos))): dataTxpnSiteUnfilteredPos.extend(ldPos[d].values())
            np.savetxt(common_path + "_intensity_brightest_txpnsite_xyfit_" + color + "_nofilter.txt", dataTxpnSiteUnfilteredPos, fmt='%.8e', delimiter='\t')

            ### filter spots

            ### remove unfitted spots from TS distribution, (only for single TS, for multiple TS see code below)
            for fn in range(0, len(lfn)):
                for nuc in ldPos[fn].keys():
                    if ldPos[fn][nuc][-1] == 0:
                        ld[fn][nuc] = []
                        ldPos[fn][nuc] = []
            #                    if params['includeMultTS'] > 0:
            #                        for TS in range(0,np.r_[ldNucTS[fn][nuc]].shape[0]):
            #                            if ldNucTSPos[fn][nuc][TS][3] == 0:
            #                                ldNucTS[fn][nuc][TS] = NaN
            #                                ldNucTSPos[fn][nuc][TS] = [NaN,NaN,NaN,NaN]

            ### remove empty cells from TS list and count list (use either own channel or other channel)
            if params['ExclEmpty'] == 1:
                for fn in range(0, len(lfn)):
                    for nuc in ld[fn].keys():
                        if ldCountBoth[fn][nuc] == 0:
                            ld[fn][nuc] = []
                            ldPos[fn][nuc] = []
                            if params['includeMultTS'] > 0:
                                ldNucTS[fn][nuc] = []
                                ldNucTSPos[fn][nuc] = []
                            ldCountNuc[fn][nuc] = []
                            ldCountBoth[fn][nuc] = []
                            ldCountCell[fn][nuc] = []
            elif params['useControlForEmptyCells'] == 1:
                if channel == params['controlColor']:
                    CountControl = ldCountBoth
                for fn in range(0, len(lfn)):
                    for nuc in ldNucTS[fn].keys():
                        if CountControl[fn][nuc] == 0:
                            ld[fn][nuc] = []
                            ldPos[fn][nuc] = []
                            if params['includeMultTS'] > 0:
                                ldNucTS[fn][nuc] = []
                                ldNucTSPos[fn][nuc] = []
                            ldCountNuc[fn][nuc] = []
                            ldCountBoth[fn][nuc] = []
                            ldCountCell[fn][nuc] = []

            cytoMean = []
            cytoSD = []
            cytoNr = []
            cytoDist = []
            cytoDistMin = []
            distAll = []
            spotAll = []
            cytoNrAll = []
            for fn in range(0, len(lfn)):
                for cyto in ldCytoPos[fn].keys():
                    if np.sum(ldCountCell[fn][cyto]) > 1:
                        cytoMean.append(sum(ldCytoPos[fn][cyto][1:, 0]) / np.r_[ldCytoPos[fn][cyto]].shape[0] - 1)
                        cytoSD.append(np.std(ldCytoPos[fn][cyto][1:, 0]))
                        cytoNr.append(np.r_[ldCytoPos[fn][cyto]].shape[0] - 1)
                        distCell = []
                        distCellMin = []
                        for i in range(1, np.r_[ldCytoPos[fn][cyto]].shape[0]):
                            distSpot = []
                            for j in range(1, np.r_[ldCytoPos[fn][cyto]].shape[0]):
                                if i != j and ldCytoPos[fn][cyto][i, -1] > 0:
                                    distCell.append(((ldCytoPos[fn][cyto][i, 1] - ldCytoPos[fn][cyto][j, 1]) ** 2 + (
                                                ldCytoPos[fn][cyto][i, 2] - ldCytoPos[fn][cyto][j, 2]) ** 2) ** 0.5)
                                    distSpot.append(((ldCytoPos[fn][cyto][i, 1] - ldCytoPos[fn][cyto][j, 1]) ** 2 + (
                                                ldCytoPos[fn][cyto][i, 2] - ldCytoPos[fn][cyto][j, 2]) ** 2) ** 0.5)
                            if distSpot != []:
                                distAll.append(min(distSpot))
                                distCellMin.append(min(distSpot))
                                spotAll.append(ldCytoPos[fn][cyto][i, 0])
                                cytoNrAll.append(np.r_[ldCytoPos[fn][cyto]].shape[0] - 1)
                        try: cytoDist.append(sum(distCell) / np.r_[distCell].shape[0])  # avg distance between spots per cell
                        except ZeroDivisionError: cytoDist.append(999999.)
                        try: cytoDistMin.append(sum(distCellMin) / np.r_[distCellMin].shape[
                            0])  # average of minimal distance between spots per cell
                        except ZeroDivisionError: cytoDistMin.append(999999.)

            tableCyto = np.c_[cytoMean, cytoSD, cytoNr, cytoDist, cytoDistMin]  # mean, SD, nr spots, avg distance, avg min distance per cell) ### Note, this data is not filtered for fitted spots
            np.savetxt(common_path + "_cyto_table_info_mean+SD+nrspots+avgdist+avgmindist_per_cell_" + color + ".txt", tableCyto, fmt='%.8e', delimiter='\t')
            tabledistAll = np.c_[distAll, spotAll, cytoNrAll]  # min distance per spot, intensity spot
            np.savetxt(common_path + "_cyto_table_info_mindist+intensity+nrspotsincell_per_spot_" + color + ".txt", tabledistAll, fmt='%.8e', delimiter='\t')

            #### define normalizing spots
            normSpots = []
            if params['AreaNorm'] == "cyto":
                normSpots = dataCell
            elif params['AreaNorm'] == "nucleus":
                normSpots = dataNucleus
            elif params['AreaNorm'] == "cell":
                for d in dataCell: normSpots.append(d)
                for d in dataNucleus: normSpots.append(d)
            elif params['AreaNorm'] == "mindist":
                idx = [a for a, b in enumerate(distAll) if b >10]
                for i in idx: normSpots.append(spotAll[i])

            ### save cytoplasmic intensity distribution
            dataCellMedian = np.median(np.array(normSpots))
            maxbin = 8 * dataCellMedian
            nobins = 100
            binvalues = np.zeros(nobins)
            for i in range(0, nobins):
                binvalues[i] = maxbin * i / nobins
            histo = np.histogram(normSpots, binvalues, density=True)
            bincenters = np.zeros(len(binvalues) - 1)
            for i in range(0, nobins - 1):
                bincenters[i] = 0.5 * (binvalues[i] + binvalues[i + 1])

            mod = GaussianModel()
            pars = mod.guess(histo[0], x=bincenters)
            out = mod.fit(histo[0], pars, x=bincenters)
            normvalue = out.params['center'].value

            if params['useFit'] == 0:
                cal = "median"
            elif params['useFit'] == 1:
                cal = "fit"

            fig = plt.figure()
            plt.bar(bincenters, histo[0], width=0.5 * binvalues[1])
            plt.plot(bincenters, out.best_fit)
            plt.axvline(x=normvalue, color='red')
            plt.title("Fit=" + str('{0:.1f}'.format(normvalue)) + ", median=" + str(
                '{0:.1f}'.format(dataCellMedian)) + ", calibration used: " + cal)
            plt.xlabel('Cytoplasmic spot intensity')
            plt.ylabel('Frequency')
            plt.xlim([0, 8 * dataCellMedian])
            plt.savefig(common_path + "_cyto_intensity_distribution_" + color + ".pdf", format='pdf')
            plt.close(fig)

            if params['useFit'] == 0:
                normCyt = dataCellMedian
            else:
                normCyt = normvalue

            np.savetxt(common_path + "_cyto_normalization_value_" + color + ".txt", np.array([normCyt]), fmt='%.8e',
                    delimiter='\t')

            ### normalize Txpn site
            dataNucleus = np.array(dataNucleus)

            dataTxpnSite = []

            for d in range(0, int(len(ld))): dataTxpnSite.extend(ld[d].values())
            dataTxpnSite = [x for x in dataTxpnSite if isinstance(x, (int, float))]

            dataNuclNorm = dataNucleus / normCyt
            dataTxpnSiteNorm = dataTxpnSite / normCyt
            dataTxpnSiteNormUnfiltered = dataTxpnSiteUnfiltered / normCyt

            # To normalize transcription sites with a user-defined value
            if params['useNormValue'] == 1 and params['includeMultTS'] == 0:
                if channel == 0:
                    normCyt = float(params['NormValueCy3'])
                elif channel == 1:
                    normCyt = float(params['NormValueCy5'])
                dataTxpnSiteNormValue = [x / normCyt for x in dataTxpnSite]
                dataTxpnSiteNormValueUnfiltered = [x / normCyt for x in dataTxpnSiteUnfiltered]
                np.savetxt(common_path + "_normalized_intensity_brightest_txpnsite_NormValue" + color + ".txt",
                        dataTxpnSiteNormValue, fmt='%.8e', delimiter='\t')
                np.save(common_path + "_normalized_intensity_brightest_txpnsite_NormValue" + color + ".npy",
                     dataTxpnSiteNormValue)
                np.savetxt(common_path + "_normalized_intensity_brightest_txpnsite_NormValue" + color + "_nofilter.txt",
                        dataTxpnSiteNormValueUnfiltered, fmt='%.8e', delimiter='\t')
                np.save(common_path + "_normalized_intensity_brightest_txpnsite_NormValue" + color + "_nofilter.npy",
                     dataTxpnSiteNormValueUnfiltered)

            ### filter and define TS for multiple TSs
            if params['includeMultTS'] > 0:
                dataTxpnSiteMult = []
                countTS = []
                dataTxpnSiteMultNorm = []
                for fn in range(0, len(ldNucTS)):
                    for nuc in ldNucTS[fn].keys():
                        TSlist = np.array(ldNucTS[fn][nuc]).reshape(-1, 1).ravel().tolist()
                        if len(TSlist) > 1: TSlist.sort()
                        if ldNucTSPos[fn][nuc] != []:
                            ldNucTSPos[fn][nuc].view('i8,i8,i8,i8,i8').sort(order=['f0'], axis=0)
                            if params['includeMultTS'] == 2:
                                if params['nrTS2Include'] > len(TSlist):
                                    for addzero in range(len(TSlist), params['nrTS2Include']): TSlist.append(0.)
                                    for addzero2 in range(ldNucTSPos[fn][nuc].reshape(-1, 5).shape[0], params['nrTS2Include']):
                                        ldNucTSPos[fn][nuc] = np.vstack((np.array([0, 0, 0, 0, 0]), ldNucTSPos[fn][nuc]))
                                fittedSpots = []
                                for fittedSpot in range(len(TSlist) - params['nrTS2Include'], len(TSlist)):
                                    if ldNucTSPos[fn][nuc][fittedSpot][-1] != 0:
                                        dataTxpnSiteMult.append(TSlist[fittedSpot])
                                        fittedSpots.append(fittedSpot)
                                    elif (ldNucTSPos[fn][nuc][fittedSpot][-1] == 0) & (
                                            ldNucTSPos[fn][nuc][fittedSpot,] == np.array([0, 0, 0, 0, 0])).all():
                                        dataTxpnSiteMult.append(TSlist[fittedSpot])
                                        fittedSpots.append(fittedSpot)
                                if len(fittedSpots) > 0:  ldNucTSPos[fn][nuc] = ldNucTSPos[fn][nuc][fittedSpots,]
                            elif params['includeMultTS'] == 1:
                                nrTS = 0
                                toRemove = []
                                for ts in range(0, len(TSlist)):
                                    if TSlist[ts] > normCyt * params['thresholdTS']:
                                        nrTS += 1
                                        if ldNucTSPos[fn][nuc].reshape(-1, 5)[ts][-1] != 0:
                                            dataTxpnSiteMult.append(TSlist[ts])
                                toRemove = (ldNucTSPos[fn][nuc].reshape(-1, 5)[:, 0] < (normCyt * params['thresholdTS'])) | (
                                            ldNucTSPos[fn][nuc].reshape(-1, 5)[:, -1] == 0)
                                toKeep = [not i for i in toRemove]
                                ldNucTSPos[fn][nuc] = ldNucTSPos[fn][nuc].reshape(-1, 5)[toKeep]
                                countTS.append(nrTS)
                dataTxpnSiteMultNorm = dataTxpnSiteMult / normCyt

            # To normalize transcription sites with a user-defined value
            if params['useNormValue'] == 1 and params['includeMultTS'] > 0:
                if channel == 0: normCyt = float(params['NormValueCy3'])
                elif channel == 1: normCyt = float(params['NormValueCy5'])
                dataTxpnSiteMultNormValue = [x / normCyt for x in dataTxpnSiteMult]
                np.savetxt(common_path + "_normalized_intensity_brightest_txpnsite_Mult_NormValue" + color + ".txt",
                           dataTxpnSiteMultNormValue, fmt='%.8e', delimiter='\t')
                np.save(common_path + "_normalized_intensity_brightest_txpnsite_Mult_NormValue" + color + ".npy",
                           dataTxpnSiteMultNormValue)

            ### convert dictionaries to np.arrays/lists

            countNucleus = []
            countCell = []
            countBoth = []

            for d in range(0, int(len(ldCountNuc))): countNucleus.extend(ldCountNuc[d].values())
            for d in range(0, int(len(ldCountCell))): countCell.extend(ldCountCell[d].values())
            for d in range(0, int(len(ldCountBoth))): countBoth.extend(ldCountBoth[d].values())

            countNucleus = [t for t in countNucleus if t != []]
            countCell = [t for t in countCell if t != []]
            countBoth = [t for t in countBoth if t != []]

            #### save results intensity spots/txpnsite, this is the filtered data. The count distribution are not filtered on whether there is a fit.
            np.savetxt(common_path + "_intensity_nuclear_spots_" + color + ".txt", dataNucleus, fmt='%.8e',
                    delimiter='\t')
            np.savetxt(common_path + "_normalized_intensity_nuclear_spots_" + color + ".txt", dataNuclNorm, fmt='%.8e',
                    delimiter='\t')
            np.savetxt(common_path + "_intensity_cytoplasmic_spots_" + color + ".txt", dataCell, fmt='%.8e',
                    delimiter='\t')
            np.savetxt(common_path + "_intensity_brightest_txpnsite_" + color + ".txt", dataTxpnSite, fmt='%.8e',
                    delimiter='\t')
            np.savetxt(common_path + "_normalized_intensity_brightest_txpnsite_" + color + ".txt", dataTxpnSiteNorm,
                    fmt='%.8e', delimiter='\t')
            np.savetxt(common_path + "_normalized_intensity_brightest_txpnsite_" + color + "_nofilter.txt",
                    dataTxpnSiteNormUnfiltered, fmt='%.8e', delimiter='\t')
            if params['includeMultTS'] > 0:
                np.savetxt(common_path + "_intensity_txpnsite_" + color + "_mult_method" + str(params['includeMultTS']) + ".txt",
                        dataTxpnSiteMult, fmt='%.8e', delimiter='\t')
                np.savetxt(common_path + "_normalized_intensity_txpnsite_" + color + "_mult_method" + str(
                    params['includeMultTS']) + ".txt", dataTxpnSiteMultNorm, fmt='%.8e', delimiter='\t')
                np.save(common_path + "_normalized_intensity_txpnsite_" + color + "_mult_method" + str(
                    params['includeMultTS']) + ".npy", dataTxpnSiteMultNorm)
                if params['includeMultTS'] == 1:
                    np.savetxt(common_path + "_count_nrTS" + color + "_threshold_" + str(params['thresholdTS']) + ".txt", countTS,
                            fmt='%.8e', delimiter='\t')
                    np.save(common_path + "_count_nrTS" + color + "_threshold_" + str(params['thresholdTS']) + ".npy", countTS)

            np.save(common_path + "_CountNuc_" + color + ".npy", countNucleus)
            np.save(common_path + "_CountCytoplasm_" + color + ".npy", countCell)
            np.save(common_path + "_CountNuc+Cytoplasm_" + color + ".npy", countBoth)
            np.save(common_path + "_normalized_intensity_txpnsite_" + color + ".npy", dataTxpnSiteNorm)

            ### calculate mask TS
            if params['CellCycle'] == 0:
                nucleusfile = get_tif(path + "_{}_max_nucleus_mask.tif".format(lfn[-1]))
            elif params['CellCycle'] == 1 and params['manualThreshold'] is None:
                nucleusfile = os.path.join(paramsexp['pathOut'], CellCycleStages[c],
                                           paramsexp['date'] + "_{}_max_nucleus_mask_threshold.tif".format(lfn[-1]))
            elif params['CellCycle'] == 1 and params['manualThreshold'] is not None:
                nucleusfile = os.path.join(paramsexp['pathOut'], 'ManualThreshold', CellCycleStages[c],
                                           paramsexp['date'] + "_{}_max_nucleus_mask_threshold.tif".format(lfn[-1]))
            with imr(nucleusfile) as im:
                size = im.shape[:2]

            for fn in range(0, len(lfn)):
                imageTS = np.zeros(size)
                squareStamp = [np.r_[
                                   -5, -5, -5, -5, -5, -5, -5, -5, -4, -3, -2, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 2, -2, -3, -4],
                               np.r_[
                                   -5, -4, -3, -2, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 4, 3, 2, -2, -3, -4, -5, -5, -5, -5, -5, -5, -5]]  # add squares around the spots in seperate channel
                if params['includeMultTS'] == 0:
                    for nuc in ldPos[fn].keys():
                        if ldPos[fn][nuc] != []:
                            if ldPos[fn][nuc][0] != 0:
                                x = int(round(ldPos[fn][nuc][1]))
                                y = int(round(ldPos[fn][nuc][2]))
                                xx = (x + squareStamp[0])[
                                    np.logical_and((y + squareStamp[1]) < 2048, (y + squareStamp[1]) >= 0)]
                                yy = (y + squareStamp[1])[
                                    np.logical_and((y + squareStamp[1]) < 2048, (y + squareStamp[1]) >= 0)]
                                xx2 = xx[np.logical_and(xx < 2048, xx >= 0)]
                                yy2 = yy[np.logical_and(xx < 2048, xx >= 0)]
                                imageTS[yy2, xx2] = 1
                #                     imageTS[ y + squareStamp[1],  x + squareStamp[0]]=1
                if params['includeMultTS'] > 0:
                    for nuc in ldNucTSPos[fn].keys():
                        if ldNucTSPos[fn][nuc] != []:
                            for ts in range(0, ldNucTSPos[fn][nuc].shape[0]):
                                x = int(round(ldNucTSPos[fn][nuc][ts][1]))
                                y = int(round(ldNucTSPos[fn][nuc][ts][2]))
                                xx = (x + squareStamp[0])[
                                    np.logical_and((y + squareStamp[1]) < 2048, (y + squareStamp[1]) >= 0)]
                                yy = (y + squareStamp[1])[
                                    np.logical_and((y + squareStamp[1]) < 2048, (y + squareStamp[1]) >= 0)]
                                xx2 = xx[np.logical_and(xx < 2048, xx >= 0)]
                                yy2 = yy[np.logical_and(xx < 2048, xx >= 0)]
                                if x != 0 and y != 0:
                                    imageTS[yy2, xx2] = 1

                if channel == 0:
                    imageTSCy3.append(imageTS)
                if channel == 1:
                    imageTSCy5.append(imageTS)

        ### add TS to mask file
        for i, fn in enumerate(tqdm(lfn, desc='Adding transcription sites to mask image')):
            path = os.path.join(paramsexp['pathOut'], paramsexp['date'])
            if not os.path.exists(path + "_{}+TS.tif".format(fn)):
                with imr(path + "_{}_max.tif".format(fn), dtype='int16') as im:
                    maskImNew = [im(channel) for channel in params['channelsToAnalyze']]
                if 0 in params['channelsToAnalyze']:
                    maskImNew.append(np.int16(imageTSCy3[i]))
                if 1 in params['channelsToAnalyze']:
                    maskImNew.append(np.int16(imageTSCy5[i]))

                tiffwrite(path + "_{}_TS.tif".format(fn), np.asarray(maskImNew), pxsize=im.pxsize, deltaz=im.deltaz)


############################################################################################################
#### makes histograms
############################################################################################################

def make_TS_count_histograms(lfn, params, paramsexp):
    print("Making histograms")

    CellCycleStages = params['CellCycleStages']
    if params['CellCycle'] == 0: Z = 1
    elif params['CellCycle'] == 1 and params['manualThreshold'] is None: Z = len(CellCycleStages)
    elif params['CellCycle'] == 1 and params['manualThreshold'] is not None: Z = 1

    for c in range(Z):
        common_path = get_common_path_out(params, paramsexp, CellCycleStages, c)
        for channel in params['channelsToAnalyze']:

            if channel == 0: color = "cy3"
            if channel == 1: color = "cy5"

            ### load data
            if params['includeMultTS'] == 0:
                try:
                    dataTS = np.load(common_path + "_normalized_intensity_txpnsite_" + color + ".npy")
                    method = "1 TS per cell"
                except Exception:
                    print("Error: No TS data found. Please run CalculateTranscriptionSites first with the same parameters")
                    continue
            elif params['includeMultTS'] == 1:
                try:
                    dataTS = np.load(common_path + "_normalized_intensity_txpnsite_" + color + "_mult_method" + str(
                    params['includeMultTS']) + ".npy")
                    method = "TS threshold " + str(params['thresholdTS'])
                    countTS = np.load(common_path + "_count_nrTS" + color + "_threshold_" + str(params['thresholdTS']) + ".npy")
                except Exception:
                    print("Error: No TS data found. Please run CalculateTranscriptionSites first with the same parameters")
                    continue
            elif params['includeMultTS'] == 2:
                try:
                    dataTS = np.load(common_path + "_normalized_intensity_txpnsite_" + color + "_mult_method" + str(
                        params['includeMultTS']) + ".npy")
                    method = str(params['nrTS2Include']) + " TSs per cell"
                except Exception:
                    print("Error: No TS data found. Please run CalculateTranscriptionSites first with the same parameters")
                    continue

            ### remove np.zeros
            dataTxpnSiteNormNonZero = [i for i in dataTS if i > (params['exclBins'] - 0.5)]  # remove TS data of excluded bins

            if channel == 0:
                if max(dataTxpnSiteNormNonZero) > params['nrbins_cy3_ts']: print("Warning! Cy3 TS fewer bins than max value")
            if channel == 1:
                if max(dataTxpnSiteNormNonZero) > params['nrbins_cy5_ts']: print("Warning! Cy5 TS fewer bins than max value")

            #### make histogram of transcription site. Negative binomial distribution as from [Raj 2006 PLoSBio] and [     ].

            if channel == 0:
                bins = np.r_[:params['nrbins_cy3_ts']]; xCenterTS = params['nrbins_cy3_ts'] / 2; freqAxisTS = params['freqAxisTS_cy3']
            elif channel == 1:
                bins = np.r_[:params['nrbins_cy5_ts']]; xCenterTS = params['nrbins_cy5_ts'] / 2; freqAxisTS = params['freqAxisTS_cy5']

            hNorm = histTSvals(dataTxpnSiteNormNonZero, bins, params['exclBins'])
            bins = bins[:-1]
            if params['CalcFracOff'] == 1 and params['exclBins'] == 0 and params['exclBinsCount'] == 0 and params['ExclEmpty'] == 1: np.savetxt(
                common_path + "_fraction_off_" + color + ".txt", [hNorm[0]], delimiter='\t')

            dataTSnonZero = []
            for i in range(0, len(dataTS)):
                if dataTS[i] > (params['exclBins'] - 1): dataTSnonZero.append(dataTS[i])

            if channel == 0: TSdistrCy3.append(dataTSnonZero)
            if channel == 1: TSdistrCy5.append(dataTSnonZero)

            #### bursting model
            if params['FitBurst'] == 1:
                [Pres, n, aT] = FitBurstModel(hNorm, bins, params['exclBins'])
            else:
                [Pres, n, aT] = [0, 0, 0]

            if params['FitPoisson'] == 1:
                [PPoissres, cT] = FitPoissonModel(hNorm, bins, params['exclBins'])
            else:
                [PPoissres, cT] = [0, 0]

            #### plot figure TS intensity distribution (and TS count)
            histTS(dataTSnonZero, hNorm[params['exclBins']:], bins[params['exclBins']:], Pres, n, aT, PPoissres, cT,
                   freqAxisTS, xCenterTS, color, params,
                   common_path + "_txpnsite_histogram_" + color + "_method=" + method.replace(" ", "_") + ".pdf")

            ### plot nr of TS per nucleus
            if params['includeMultTS'] == 1:
                if channel == 0:
                    bins = np.r_[:params['nrbins_cy3_ts']]; xCenter = params['nrbins_cy3_ts'] / 2
                elif channel == 1:
                    bins = np.r_[:params['nrbins_cy5_ts']]; xCenter = params['nrbins_cy5_ts'] / 2
                hist_TS_count = np.histogram(countTS, bins=bins - .5, density=1)[0]
                bins = bins[:-1]
                fig = plt.figure()
                plt.bar(bins, hist_TS_count)
                if freqAxisTS != "None": plt.ylim([0, freqAxisTS])
                plt.title("Number of TS in the nucleus, " + color)
                plt.xlabel('Nr TS/nucleus')
                plt.ylabel('Frequency')
                plt.text(xCenter, max(hist_TS_count),
                         "Mean TS count = " + str('{0:.2f}'.format(np.mean(countTS))) + " +/- " + str('{0:.2f}'.format(
                             (np.std(countTS) / np.sqrt(len(countTS))))) + "\nNumber of cells analyzed = " + str(
                             len(countTS)) + "\nThreshold TS = " + str(params['thresholdTS']), fontsize=12,
                         horizontalalignment='center', verticalalignment='top')
                plt.savefig(common_path + "_TS_count_" + color + "_threshold=" + str(params['thresholdTS']) + ".pdf",
                            format='pdf')
                plt.close(fig)

                # matplotlib.pyplot.violinplot(dataset, positions=None, vert=True, widths=0.5, showmeans=False, showextrema=True, showmedians=False, points=100, bw_method=None, *, data=None)[source]

            ### counting spots

            countNucleus = np.load(common_path + "_CountNuc_" + color + ".npy")
            countCell = np.load(common_path + "_CountCytoplasm_" + color + ".npy")
            countBoth = np.load(common_path + "_CountNuc+Cytoplasm_" + color + ".npy")

            if channel == 0:
                if max(countNucleus) > params['nrbins_cy3']: print("Warning! Cy3 Nuclear spots - fewer bins than max value")
                if max(countBoth) > params['nrbins_cy3']: print("Warning! Cy3 Cellular spots - fewer bins than max value")
            if channel == 1:
                if max(countNucleus) > params['nrbins_cy5']: print("Warning! Cy5 Nuclear spots - fewer bins than max value")
                if max(countBoth) > params['nrbins_cy5']: print("Warning! Cy5 Cellular spots - fewer bins than max value")

            if channel == 0:
                bins = np.r_[
                       :params['nrbins_cy3']]; xCenter = params['nrbins_cy3'] / 2; freqAxisCell = params['freqAxisCell_cy3']; freqAxisNuc = params['freqAxisNuc_cy3']
            elif channel == 1:
                bins = np.r_[
                       :params['nrbins_cy5']]; xCenter = params['nrbins_cy5'] / 2; freqAxisCell = params['freqAxisCell_cy5']; freqAxisNuc = params['freqAxisNuc_cy3']

            totalCellsCounted = len(countBoth)
            countBothNonZero = [i for i in countBoth if i > (params['exclBinsCount'] - 1)]
            totalCellsCountedNonzero = len(countBothNonZero)
            countNucleusNonZero = [i for i in countNucleus if i > (params['exclBinsCount'] - 1)]
            totalNucleusCountedNonzero = len(countNucleusNonZero)

            if channel == 0:
                CelldistrCy3.append(countBothNonZero)
                NucldistrCy3.append(countNucleusNonZero)
            if channel == 1:
                CelldistrCy5.append(countBothNonZero)
                NucldistrCy5.append(countNucleusNonZero)

            histNuc(countNucleusNonZero, bins, freqAxisNuc, xCenter, params['exclBinsCount'], color,
                    common_path + "_nuclear_count_" + color + ".pdf")
            histCell(countBothNonZero, bins, freqAxisCell, xCenter, params['exclBinsCount'], color,
                     common_path + "_cell_count_" + color + ".pdf")

            np.savetxt(common_path + "_cellular_count_" + color + ".txt", countBoth, fmt='%i', delimiter='\t')
            np.savetxt(common_path + "_nuclear_count_" + color + ".txt", countNucleus, fmt='%i', delimiter='\t')

            writestatsfile(common_path, color, countBothNonZero, countNucleusNonZero, dataTSnonZero,
                           params)

            if params['useNormValue'] == 1 and params['includeMultTS'] == 0:
                    TSNormValuedistrCy3 = []
                    TSNormValuedistrCy5 = []
                    dataTSNormValue = np.load(
                        common_path + "_normalized_intensity_brightest_txpnsite_NormValue" + color + ".npy")
                    method = "define NormValue"
                    dataTxpnSiteNormValueNonZero = [i for i in dataTSNormValue if
                                                    i > (params['exclBins'] - 0.5)]  # remove TS data of excluded bins
                    if channel == 0:
                        if max(
                            dataTxpnSiteNormValueNonZero) > params['nrbins_cy3_ts']: print("Warning! Cy3 TS fewer bins than max value")
                    if channel == 1:
                        if max(
                            dataTxpnSiteNormValueNonZero) > params['nrbins_cy5_ts']: print("Warning! Cy5 TS fewer bins than max value")
                    if channel == 0:
                        bins = np.r_[:params['nrbins_cy3_ts']]; xCenterTS = params['nrbins_cy3_ts'] / 2; freqAxisTS = params['freqAxisTS_cy3']
                    elif channel == 1:
                        bins = np.r_[:params['nrbins_cy5_ts']]; xCenterTS = params['nrbins_cy5_ts'] / 2; freqAxisTS = params['freqAxisTS_cy5']
                    hNorm = histTSvals(dataTxpnSiteNormValueNonZero, bins, params['exclBins'])
                    bins = bins[:-1]
                    dataTSNormValuenonZero = []
                    for i in range(0, len(dataTSNormValue)):
                        if dataTSNormValue[i] > (params['exclBins'] - 1): dataTSNormValuenonZero.append(dataTSNormValue[i])
                    if channel == 0: TSNormValuedistrCy3.append(dataTSNormValuenonZero)
                    if channel == 1: TSNormValuedistrCy5.append(dataTSNormValuenonZero)
                    histTS(dataTSNormValuenonZero, hNorm[params['exclBins']:], bins[params['exclBins']:], Pres, n, aT,
                           PPoissres, cT, freqAxisTS, xCenterTS, color, params,
                           common_path + "_txpnsite_histogram" + color + "_method=" + method.replace(" ", "_") + ".pdf")

            if params['useNormValue'] == 1 and params['includeMultTS'] > 0:
                    MultTSNormValuedistrCy3 = []
                    MultTSNormValuedistrCy5 = []
                    dataTSMultNormValue = np.load(common_path + "_normalized_intensity_brightest_txpnsite_Mult_NormValue" + color + ".npy")
                    method = "define NormValue"
                    dataTxpnSiteMultNormValueNonZero = [i for i in dataTSMultNormValue if i > (params['exclBins'] - 0.5)]  # remove TS data of excluded bins
                    if channel == 0:
                      if max(dataTxpnSiteMultNormValueNonZero) > params['nrbins_cy3_ts']: print("Warning! Cy3 TS fewer bins than max value")
                    if channel == 1:
                      if max(dataTxpnSiteMultNormValueNonZero) > params['nrbins_cy5_ts']: print("Warning! Cy5 TS fewer bins than max value")
                    if channel == 0:
                      bins = np.r_[:params['nrbins_cy3_ts']]; xCenterTS = params['nrbins_cy3_ts'] / 2; freqAxisTS = params['freqAxisTS_cy3']
                    elif channel == 1:
                      bins = np.r_[:params['nrbins_cy5_ts']]; xCenterTS = params['nrbins_cy5_ts'] / 2; freqAxisTS = params['freqAxisTS_cy5']
                    hNorm = histTSvals(dataTxpnSiteMultNormValueNonZero, bins, params['exclBins'])
                    bins = bins[:-1]
                    dataTSMultNormValuenonZero = []
                    for i in range(0, len(dataTSMultNormValue)):
                      if dataTSMultNormValue[i] > (params['exclBins'] - 1): dataTSMultNormValuenonZero.append(dataTSMultNormValue[i])
                    if channel == 0: MultTSNormValuedistrCy3.append(dataTSMultNormValuenonZero)
                    if channel == 1: MultTSNormValuedistrCy5.append(dataTSMultNormValuenonZero)
                    histTS(dataTSMultNormValuenonZero, hNorm[params['exclBins']:], bins[params['exclBins']:], Pres, n, aT,
                         PPoissres, cT, freqAxisTS, xCenterTS, color, params,
                         common_path + "_txpnsite_histogram" + color + "_method=" + method.replace(" ", "_") + ".pdf")

    return (CelldistrCy3, NucldistrCy3, TSdistrCy3, CelldistrCy5, NucldistrCy5, TSdistrCy5)


############################################################################################################
#### calculates correlations between channels and between cell and nascent RNA number
############################################################################################################

def calculate_singlecell_corr(lfn, params, paramsexp):
    print("Calculating correlation")
    cors = []

    CellCycleStages = params['CellCycleStages']
    if params['CellCycle'] == 0: Z = 1
    elif params['CellCycle'] == 1 and params['manualThreshold'] is None: Z = len(CellCycleStages)
    elif params['CellCycle'] == 1 and params['manualThreshold'] is not None: Z = 1

    for c in range(Z):
        common_path = get_common_path_out(params, paramsexp, CellCycleStages, c)

        ### Scatter plot between nr of nascent RNA and RNA count Cy3
        if 0 in params['channelsToAnalyze']:
            CountCy3 = np.loadtxt(common_path + "_CountNuc+Cytoplasm_cy3_nofilter.txt")
            TSCy3 = np.loadtxt(common_path + "_normalized_intensity_brightest_txpnsite_cy3_nofilter.txt")
            if params['useNormValue'] == 1: TSCy3_NormValue = np.loadtxt(common_path + "_normalized_intensity_brightest_txpnsite_NormValuecy3_nofilter.txt")
            TSCy3Pos = np.loadtxt(common_path + "_intensity_brightest_txpnsite_xyfit_cy3_nofilter.txt")
            if params['ExclEmpty'] == 1:
                EmptyCy3 = list(np.where(CountCy3 == 0)[0])  # filter on cells that are empty
            else:
                EmptyCy3 = []
            NoFitCy3 = list(np.where(TSCy3Pos[:, -1] == 0)[0])  # filter on spots np.where there was not fit in Cy3
            NoFit_EmptyCy3 = set(EmptyCy3 + NoFitCy3)  # filter on both
            NotEmptyCy3 = [x for x in range(0, len(CountCy3)) if x not in NoFit_EmptyCy3]
            corTScountCy3 = pearsonr(CountCy3[NotEmptyCy3], TSCy3[NotEmptyCy3])  # correlate nr of RNA and TS intensity
            fig = plt.figure()
            plt.scatter(CountCy3[NotEmptyCy3], TSCy3[NotEmptyCy3], s=2)
            plt.title("Cy3, correlation = " + str('{0:.3f}'.format(corTScountCy3[0])))
            plt.xlabel('Nr RNA/cell, Cy3')
            plt.ylabel('Nr of nascent RNA at TS')
            plt.savefig(common_path + "_correlation_count_with_TSintensity_cy3.pdf", format='pdf')
            plt.close(fig)
            cors.append(
                "Correlation Cy3 TS intensity and cell count = " + str(corTScountCy3[0]) + " , p-value = " + str(
                    '{:.2E}'.format(corTScountCy3[1])))
            cors.append(
                str(len(NotEmptyCy3)) + " cells were analyzed, " + str(len(NoFit_EmptyCy3)) + " cells were filtered")
            cors.append(str(len(EmptyCy3)) + " cells removed because they were empty in Cy3 channel")
            cors.append(str(len(NoFitCy3)) + " cells were removed because the Cy3 TS were not fit")

        ### Scatter plot between nr of nascent RNA and RNA count Cy5
        if 1 in params['channelsToAnalyze']:
            CountCy5 = np.loadtxt(common_path + "_CountNuc+Cytoplasm_cy5_nofilter.txt")
            TSCy5 = np.loadtxt(common_path + "_normalized_intensity_brightest_txpnsite_cy5_nofilter.txt")
            if params['useNormValue'] == 1: TSCy5_NormValue = np.loadtxt(common_path + "_normalized_intensity_brightest_txpnsite_NormValuecy5_nofilter.txt")
            TSCy5Pos = np.loadtxt(common_path + "_intensity_brightest_txpnsite_xyfit_cy5_nofilter.txt")
            if params['ExclEmpty'] == 1:
                EmptyCy5 = list(np.where(CountCy5 == 0)[0])  # filter on cells that are empty
            else:
                EmptyCy5 = []
            NoFitCy5 = list(np.where(TSCy5Pos[:, -1] == 0)[0])  # filter on spots np.where there was not fit in Cy5
            NoFit_EmptyCy5 = set(EmptyCy5 + NoFitCy5)  # filter on both
            NotEmptyCy5 = [x for x in range(0, len(CountCy5)) if
                           x not in NoFit_EmptyCy5]  # correlate nr of RNA and TS intensity
            corTScountCy5 = pearsonr(CountCy5[NotEmptyCy5], TSCy5[NotEmptyCy5])
            fig = plt.figure()
            plt.scatter(CountCy5[NotEmptyCy5], TSCy5[NotEmptyCy5], s=2)
            plt.title("Cy5, correlation = " + str('{0:.3f}'.format(corTScountCy5[0])))
            plt.xlabel('Nr RNA/cell, Cy5')
            plt.ylabel('Nr of nascent RNA at TS')
            plt.savefig(common_path + "_correlation_count_with_TSintensity_cy5.pdf", format='pdf')
            plt.close(fig)
            cors.append(
                "\nCorrelation Cy5 TS intensity and cell count = " + str(corTScountCy5[0]) + " , p-value = " + str(
                    '{:.2E}'.format(corTScountCy5[1])))
            cors.append(
                str(len(NotEmptyCy5)) + " cells were analyzed, " + str(len(NoFit_EmptyCy5)) + " cells were filtered")
            cors.append(str(len(EmptyCy5)) + " cells removed because they were empty in Cy5 channel")
            cors.append(str(len(NoFitCy5)) + " cells were removed because the Cy5 TS was not fit")

        ### Correlate between channels
        if len(params['channelsToAnalyze']) < 2:
            print("Cannot calculate Cy3-Cy5 correlations in 1 color data")
        elif len(params['channelsToAnalyze']) == 2:
            ### filter empty cells and spots that were not fit in either channel
            if params['ExclEmpty'] == 1:
                overlapEmpty = [x for x in EmptyCy3 if
                                x in EmptyCy5]  #### only trows away cells np.where no spots in both Cy3 and Cy5.
            else:
                overlapEmpty = []
            NoFitCy3Cy5 = list(set(NoFitCy3 + NoFitCy5))  # filter on cells that do not have a fit in either Cy3 or Cy5
            NoFit_EmptyCy3Cy5 = list(set(overlapEmpty + NoFitCy3Cy5))
            #           overlapNotEmpty = [x for x in range(0,len(CountCy3)) if x not in NoFit_EmptyCy3Cy5]

            ### calculate and plot distances
            if params['localize']['dimension'] == 2:
                dist = ((TSCy3Pos[:, 1] - TSCy5Pos[:, 1]) ** 2 + (TSCy3Pos[:, 2] - TSCy5Pos[:, 2]) ** 2) ** 0.5
            else:
                dist = ((TSCy3Pos[:, 1] - TSCy5Pos[:, 1]) ** 2 + (TSCy3Pos[:, 2] - TSCy5Pos[:, 2]) ** 2 +
                        (TSCy3Pos[:, 3] - TSCy5Pos[:, 3]) ** 2) ** 0.5
            np.savetxt(common_path + "_distances_TSs.txt", dist, delimiter="\n", fmt='%.8e')

            fig = plt.figure(figsize=(8, 3))
            h, x = np.histogram(np.sqrt(dist), bins=np.r_[0:20:.2])
            plt.xlabel('TS distances')
            plt.ylabel('Count')
            plt.title("Distances between Cy3 and Cy5 TS")
            plt.plot(x[:-1], h + 1)
            plt.tight_layout()
            plt.savefig(common_path + "_distances_TS.pdf", format='pdf')
            plt.close(fig)

            ### filter on distance
            if params['filterOnDistance'] == 1:
                filterDist = list(np.where(dist >= params['distThresh'])[0])
            else:
                filterDist = []

            toFilter = set(NoFit_EmptyCy3Cy5 + filterDist)
            toAnalyze = [x for x in range(0, len(CountCy3)) if x not in toFilter]

            #### only select expressing population
            fracCy3OnlyOff = sum((TSCy3[toAnalyze] < params['onoffThreshCy3']) & (TSCy5[toAnalyze] >= params['onoffThreshCy5']))
            fracCy5OnlyOff = sum((TSCy3[toAnalyze] >= params['onoffThreshCy3']) & (TSCy5[toAnalyze] < params['onoffThreshCy5']))
            fracBothOff = sum((TSCy3[toAnalyze] < params['onoffThreshCy3']) & (TSCy5[toAnalyze] < params['onoffThreshCy5']))
            fracBothOn = sum((TSCy3[toAnalyze] >= params['onoffThreshCy3']) & (TSCy5[toAnalyze] >= params['onoffThreshCy5']))
            fractionList = [["Only Cy3 off", "Only Cy5 off", "Both off", "Both on"],
                            [fracCy3OnlyOff, fracCy5OnlyOff, fracBothOff, fracBothOn],
                            [float(fracCy3OnlyOff) / float(len(toAnalyze)),
                             float(fracCy5OnlyOff) / float(len(toAnalyze)), float(fracBothOff) / float(len(toAnalyze)),
                             float(fracBothOn) / float(len(toAnalyze))]]
            np.savetxt(common_path + "_fraction_on_off.txt", fractionList, delimiter="\t", fmt="%s")

            if params['useNormValue'] == 1:
                fracCy3OnlyOff = sum((TSCy3_NormValue[toAnalyze] < params['onoffThreshCy3']) & (TSCy5_NormValue[toAnalyze] >= params['onoffThreshCy5']))
                fracCy5OnlyOff = sum((TSCy3_NormValue[toAnalyze] >= params['onoffThreshCy3']) & (TSCy5_NormValue[toAnalyze] < params['onoffThreshCy5']))
                fracBothOff = sum((TSCy3_NormValue[toAnalyze] < params['onoffThreshCy3']) & (TSCy5_NormValue[toAnalyze] < params['onoffThreshCy5']))
                fracBothOn = sum((TSCy3_NormValue[toAnalyze] >= params['onoffThreshCy3']) & (TSCy5_NormValue[toAnalyze] >= params['onoffThreshCy5']))
                fractionList = [["Only Cy3 off", "Only Cy5 off", "Both off", "Both on"],
                                [fracCy3OnlyOff, fracCy5OnlyOff, fracBothOff, fracBothOn],
                                [float(fracCy3OnlyOff) / float(len(toAnalyze)),
                                 float(fracCy5OnlyOff) / float(len(toAnalyze)), float(fracBothOff) / float(len(toAnalyze)),
                                 float(fracBothOn) / float(len(toAnalyze))]]
                np.savetxt(common_path + "_fraction_on_off_NormValue.txt", fractionList, delimiter="\t", fmt="%s")


            overlap = [value for value in np.where(TSCy3 >= params['onoffThreshCy3'])[0].tolist() if
                       value in set(np.where(TSCy5 >= params['onoffThreshCy5'])[0].tolist())]
            overlap2 = [x for x in overlap if x not in toFilter]

            corTSCy3Cy5 = pearsonr(TSCy3[toAnalyze], TSCy5[toAnalyze]) if toAnalyze else (np.nan, np.nan)
            corTSCy3Cy5bg = pearsonr(TSCy3[overlap2], TSCy5[overlap2]) if overlap2 else (np.nan, np.nan)
            corCountCy3Cy5 = pearsonr(CountCy3[toAnalyze], CountCy5[toAnalyze]) if toAnalyze else (np.nan, np.nan)
            if params['useNormValue'] == 1:
                corTSCy3Cy5NormValue = pearsonr(TSCy3_NormValue[toAnalyze], TSCy5_NormValue[toAnalyze])

            fig = plt.figure()
            plt.scatter(CountCy3[toAnalyze], CountCy5[toAnalyze], s=2)
            plt.title("Correlation = " + str('{0:.3f}'.format(corCountCy3Cy5[0])) + " , p-value = " + str(
                '{:.2E}'.format(corCountCy3Cy5[1])))
            plt.xlabel('Nr RNA/cell Cy3')
            plt.ylabel('Nr RNA/cell Cy5')
            plt.savefig(common_path + "_correlation_Count_Nuc+Cytoplasm.pdf", format='pdf')
            plt.close(fig)

            fig = plt.figure()
            plt.scatter(TSCy3[toAnalyze], TSCy5[toAnalyze], s=2)
            plt.title("Correlation = " + str('{0:.3f}'.format(corTSCy3Cy5[0])) + " , p-value = " + str(
                '{:.2E}'.format(corTSCy3Cy5[1])))
            plt.xlabel('Nr of nascent RNA at TS, Cy3')
            plt.ylabel('Nr of nascent RNA at TS, Cy5')
            plt.xlim(0, params['x_lim'])
            plt.ylim(0, params['y_lim'])
            plt.savefig(common_path + "_correlation_TSintensity.pdf", format='pdf')
            plt.close(fig)

            fig = plt.figure()
            plt.scatter(TSCy3[overlap2], TSCy5[overlap2], s=2)
            plt.title("Correlation = " + str('{0:.3f}'.format(corTSCy3Cy5bg[0])) + " , p-value = " + str(
                '{:.2E}'.format(corTSCy3Cy5bg[1])))
            plt.xlabel('Nr of nascent RNA at TS, Cy3')
            plt.ylabel('Nr of nascent RNA at TS, Cy5')
            plt.xlim(0, params['x_lim'])
            plt.ylim(0, params['y_lim'])
            plt.savefig(common_path + "_correlation_TSintensity_abovethreshold.pdf", format='pdf')
            plt.close(fig)

            if params['useNormValue'] == 1:
                fig = plt.figure()
                plt.scatter(TSCy3_NormValue[toAnalyze], TSCy5_NormValue[toAnalyze], s=2)
                plt.title("Correlation = " + str('{0:.3f}'.format(corTSCy3Cy5NormValue[0])) + " , p-value = " + str(
                    '{:.2E}'.format(corTSCy3Cy5NormValue[1])))
                plt.xlabel('Nr of nascent RNA at TS, Cy3')
                plt.ylabel('Nr of nascent RNA at TS, Cy5')
                plt.xlim(0, params['x_lim'])
                plt.ylim(0, params['y_lim'])
                plt.savefig(common_path + "_correlation_TSintensity_NormValue.pdf", format='pdf')
                plt.close(fig)
                np.savetxt(common_path + "_correlation_TSCy3_NormValue.txt", TSCy3_NormValue[toAnalyze], delimiter="\n", fmt="%s")
                np.savetxt(common_path + "_correlation_TSCy5_NormValue.txt", TSCy5_NormValue[toAnalyze], delimiter="\n", fmt="%s")

            cors.append("\nCorrelation Cy3-Cy5 TS intensity = " + str(corTSCy3Cy5[0]) + " , p-value = " + str(
                '{:.2E}'.format(corTSCy3Cy5[1])))
            cors.append("Correlation Cy3-Cy5 TS intensity (excluded below thresholds " + str(
                params['onoffThreshCy3']) + " for cy3, " + str(params['onoffThreshCy5']) + " for cy5) = " + str(
                corTSCy3Cy5bg[0]) + " , p-value = " + str('{:.2E}'.format(corTSCy3Cy5bg[1])))
            cors.append("Correlation Cy3-Cy5 cell count = " + str(corCountCy3Cy5[0]) + " , p-value = " + str(
                '{:.2E}'.format(corCountCy3Cy5[1])))
            cors.append(str(len(toAnalyze)) + "cell were analyzed, " + str(len(toFilter)) + " cells were filtered")
            cors.append(str(len(overlapEmpty)) + " cells removed because they were empty in both Cy3 and Cy5 channel")
            cors.append(str(len(
                NoFitCy3Cy5)) + " cells were removed because they were not fit in either Cy3 or Cy5 of the channels")

            TSCy3 = TSCy3[toAnalyze]
            TSCy5 = TSCy5[toAnalyze]
            np.savetxt(common_path + "_correlation_TSCy3.txt", TSCy3, delimiter="\n", fmt="%s")
            np.savetxt(common_path + "_correlation_TSCy5.txt", TSCy5, delimiter="\n", fmt="%s")
            threshCy3 = range(0, int(max(TSCy3)))
            threshCy5 = range(0, int(max(TSCy5)))
            corMatrixL = np.zeros((len(threshCy3), len(threshCy5)))
            corMatrixLS = np.zeros((len(threshCy3), len(threshCy5)))
            corMatrixS = np.zeros((len(threshCy3), len(threshCy5)))
            corMatrixSL = np.zeros((len(threshCy3), len(threshCy5)))
            for m in threshCy3:
                for n in threshCy5:
                    overlapL = [value for value in np.where(TSCy3 >= m)[0].tolist() if
                                value in set(np.where(TSCy5 >= n)[0].tolist())]
                    overlapS = [value for value in np.where(TSCy3 < m)[0].tolist() if
                                value in set(np.where(TSCy5 < n)[0].tolist())]
                    #                    non_overlapL = range(0,len(TSCy3))
                    #                    for i in overlapL: non_overlapL.remove(i)
                    non_overlapS = list(range(0, len(TSCy3)))
                    for i in overlapS: non_overlapS.remove(i)

                    if len(overlapL) > 10: corMatrixL[m, n] = pearsonr(TSCy3[overlapL], TSCy5[overlapL])[0]
                    #                    if len(overlapS)>10: corMatrixS[m,n] = pearsonr(TSCy3[overlapS], TSCy5[overlapS])[0]
                    #                    if len(non_overlapL)>10: corMatrixLS[m,n] = pearsonr(TSCy3[non_overlapL], TSCy5[non_overlapL])[0]
                    if len(non_overlapS) > 10: corMatrixSL[m, n] = pearsonr(TSCy3[non_overlapS], TSCy5[non_overlapS])[0]

            #      fig = plt.figure()
            fig, (ax1, ax2) = plt.subplots(1, 2)
            im1 = ax1.pcolor(corMatrixL, cmap="RdBu_r", vmin=-1, vmax=1)
            ax1.set_title("Corr TS intensity Cy3-Cy5, \nAND>=")
            ax1.set_xlabel("larger than threshold, Cy5")
            ax1.set_ylabel("larger than threshold, Cy3")

            im2 = ax2.pcolor(corMatrixSL, cmap="RdBu_r", vmin=-1, vmax=1)
            ax2.set_title("Corr TS intensity Cy3-Cy5, \nOR>=")
            ax2.set_xlabel("larger than threshold, Cy5")
            #            ax2.set_ylabel("larger than threshold, Cy3")

            divider = make_axes_locatable(ax2)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im2, cax=cax)

            #            fig.tight_layout()
            plt.savefig(common_path + "_correlation_TS_intensity_heatmaps_threshold.pdf", format='pdf')

            np.savetxt(common_path + "_correlations.txt", cors, delimiter="\n", fmt="%s")


############################################################################################################
#### combine replicates
############################################################################################################

def combine_reps(distributions, params):
    print("Combining datasets")
    distributionsnew = [[], [], [], [], [], []]

    for channel in params['channelsToAnalyze']:
        if channel == 0: color = "cy3"
        if channel == 1: color = "cy5"

        for k in [0, 1, 2]:
            distributionselement = []
            for i in range(0, len(params['repsToCombine'])):
                combi = []
                for j in params['repsToCombine'][i]:
                    combi.append(distributions[3 * channel + k][j - 1])
                distributionselement.append(np.concatenate(combi))
            distributionsnew[3 * channel + k] = distributionselement
    return distributionsnew

    print("Making histograms of combined datasets")

    shutil.copyfile("FISH_pipeline_parameters.yml", params['analysisfolder'] + "_FISH_pipeline_parameters.yml")

    for channel in params['channelsToAnalyze']:
        if params['CalcFracOff'] == 1 and params['exclBins'] == 0 and params['exclBinsCount'] == 0 and params['ExclEmpty'] == 1: FracOff = []
        for i in range(0, len(params['repsToCombine'])):

            ### for TS distribution
            if channel == 0: color = "cy3"; bins = np.r_[
                                                   :params['nrbins_cy3_ts']]; xCenterTS = params['nrbins_cy3_ts'] / 2; freqAxisTS = params['freqAxisTS_cy3']
            if channel == 1: color = "cy5"; bins = np.r_[
                                                   :params['nrbins_cy5_ts']]; xCenterTS = params['nrbins_cy5_ts'] / 2; freqAxisTS = params['freqAxisTS_cy5']

            hNorm = histTSvals(distributions[3 * channel + 2][i], bins, params['exclBins'])
            bins = bins[:-1]

            TSdistr = distributions[3 * channel + 2][i]
            if params['CalcFracOff'] == 1 and params['exclBins'] == 0 and params['exclBinsCount'] == 0 and params['ExclEmpty'] == 1: FracOff.append(hNorm[0])

            #### Fitting models
            if params['FitBurst'] == 1:
                [Pres, n, aT] = FitBurstModel(hNorm, bins, params['exclBins'])
            else:
                [Pres, n, aT] = [0, 0, 0]

            if params['FitPoisson'] == 1:
                [PPoissres, cT] = FitPoissonModel(hNorm, bins, params['exclBins'])
            else:
                [PPoissres, cT] = [0, 0]

            #### making TS plot

            histTS(TSdistr, hNorm[params['exclBins']:], bins[params['exclBins']:], Pres, n, aT, PPoissres,
                   cT, freqAxisTS, xCenterTS, color, params,
                   params['analysisfolder'] + "_txpnsite_histogram_" + color + "_merged_distr_" + str(i) + ".pdf")

            ## nuclear and cellular count distribution

            if channel == 0: bins = np.r_[:params['nrbins_cy3']]; freqAxisNuc = params['freqAxisNuc_cy3']; freqAxisCell = params['freqAxisCell_cy3']; xCenter = params['nrbins_cy3'] / 2
            elif channel == 1: bins = np.r_[:params['nrbins_cy5']]; freqAxisNuc = params['freqAxisNuc_cy5']; freqAxisCell = params['freqAxisCell_cy5']; xCenter = params['nrbins_cy5'] / 2

            histNuc(distributions[3 * channel + 1][i], bins, freqAxisNuc, xCenter, params['exclBinsCount'], color,
                    params['analysisfolder'] + "_nuclear_count_" + color + "_merged_distr_" + str(i) + ".pdf")

            histCell(distributions[3 * channel][i], bins, freqAxisCell, xCenter, params['exclBinsCount'], color,
                     params['analysisfolder'] + "_cell_count_" + color + "_merged_distr_" + str(i) + ".pdf")

            ## write statistics file
            writestatsfile(params['analysisfolder'] + "_distr_" + str(i), color, distributions[3 * channel][i],
                           distributions[3 * channel + 1][i], TSdistr, params)

        ##write text files
        np.savetxt(params['analysisfolder'] + "_txpnsite_" + color + "_merged_distr_" + str(i) + ".txt", TSdistr, fmt='%.8e',
                delimiter='\t')
        np.savetxt(params['analysisfolder'] + "_cell_count_" + color + "_merged_distr_" + str(i) + ".txt",
                distributions[3 * channel][i], fmt='%.8e', delimiter='\t')

        if params['CalcFracOff'] == 1 and params['exclBins'] == 0 and params['exclBinsCount'] == 0 and params['ExclEmpty'] == 1: np.savetxt(
            params['analysisfolder'] + "fraction_off_" + color + ".txt", FracOff, delimiter='\t')


############################################################################################################
#### makes violin and boxplots
############################################################################################################

def make_violin_plots(distributions, params):
    print("Making violin and boxplots")

    for channel in params['channelsToAnalyze']:
        if channel == 0: color = "cy3"; yAxisViolinCell = params['yAxisViolinCell_cy3']; yAxisViolinNuc = params['yAxisViolinNuc_cy3']; yAxisViolinTS = params['yAxisViolinTS_cy3']
        if channel == 1: color = "cy5"; yAxisViolinCell = params['yAxisViolinCell_cy5']; yAxisViolinNuc = params['yAxisViolinNuc_cy5']; yAxisViolinTS = params['yAxisViolinTS_cy5']

        countdatasets = len(distributions[3 * channel])
        if countdatasets != 1 and len(params['pvalsViolins'][0]) != 0:
            pcell = []
            pcellround = []
            pnuc = []
            pnucround = []
            pTS = []
            pTSround = []
            for ppair in range(0, len(params['pvalsViolins'])):
                dist1 = params['pvalsViolins'][ppair][0] - 1
                dist2 = params['pvalsViolins'][ppair][1] - 1
                z, p = scipy.stats.mannwhitneyu(distributions[3 * channel][dist1], distributions[3 * channel][dist2])
                p = 2 * p
                pcell.append(p)
                pcellround.append(round(p, 4))
                z, p = scipy.stats.mannwhitneyu(distributions[3 * channel + 1][dist1],
                                                distributions[3 * channel + 1][dist2])
                p = 2 * p
                pnuc.append(p)
                pnucround.append(round(p, 4))
                z, p = scipy.stats.mannwhitneyu(distributions[3 * channel + 2][dist1],
                                                distributions[3 * channel + 2][dist2])
                p = 2 * p
                pTS.append(p)
                pTSround.append(round(p, 4))

        datasets = []
        for file in range(0, len(params['fileListIn'])):
            pathIn = params['fileListIn'][file]
            date = pathIn.split("/")[-3]
            expName = pathIn.split("/")[-2]
            datasets.append([date, expName])

        if params['CombineReps'] == 0:
            datasetLabels = []
            for i in range(len(datasets)):
                datasetLabels.append("\n".join(datasets[i]))
        else:
            datasetLabels = []
            for i in range(len(params['repsToCombine'])):
                dates = []

                for j in params['repsToCombine'][i]:
                    dates.append("_".join(datasets[j - 1]))
                datasetLabels.append("\n".join(dates))


        plt.rcParams.update({'figure.autolayout': True})
        for figtype in ["Violin", "Boxplot"]:
            fig = plt.figure(figsize=(11.69 / 2., 8.27))
            if figtype == "Violin":
                plt.violinplot(distributions[3 * channel], showmeans=True, showextrema=False)
                plt.xticks(np.arange(1, countdatasets+1), np.array(datasetLabels), rotation=45,
                           horizontalalignment="right")
            elif figtype == "Boxplot":
                sns.boxplot(data=distributions[3 * channel ], color=".8")
                plt.xticks(np.arange(0, countdatasets), np.array(datasetLabels), rotation=45, horizontalalignment="right")
            plt.ylabel('Number of RNA/Cell')
            plt.xlabel('Dataset')
            y_max = max(np.concatenate(distributions[3 * channel]))
            x_max = len(distributions[3 * channel])
            if yAxisViolinCell != "None": y_max = yAxisViolinCell
            plt.ylim([0, y_max])
            if countdatasets != 1 and len(params['pvalsViolins'][0]) != 0:
                plt.text(x_max + 0.3, 0.98 * y_max,
                         "2-sided Mann-Whitney U-tests:\npvalues calculated for datasets: " + str(
                             params['pvalsViolins']) + "\npvalues: " + str(pcellround), fontsize=12, horizontalalignment='right',
                         verticalalignment='top')
            plt.savefig(params['analysisfolder']  + figtype + "_Cellular_" + color + "_combined_reps" * params['CombineReps'] +".pdf")
            plt.close(fig)

        for figtype in ["Violin", "Boxplot"]:
            fig = plt.figure(figsize=(11.69 / 2., 8.27))
            if figtype == "Violin":
                plt.violinplot(distributions[3 * channel + 1], showmeans=True, showextrema=False)
                plt.xticks(np.arange(1, countdatasets+1), np.array(datasetLabels), rotation=45,
                           horizontalalignment="right")
            elif figtype == "Boxplot":
                sns.boxplot(data=distributions[3 * channel + 1], color=".8")
                plt.xticks(np.arange(0, countdatasets), np.array(datasetLabels), rotation=45,
                           horizontalalignment="right")
            plt.ylabel('Number of RNA/Nucleus')
            plt.xlabel('Dataset')
            y_max = max(np.concatenate(distributions[3 * channel + 1]))
            x_max = len(distributions[3 * channel + 1])
            if yAxisViolinNuc != "None": y_max = yAxisViolinNuc
            plt.ylim([0, y_max])
            if countdatasets != 1 and len(params['pvalsViolins'][0]) != 0:
                plt.text(x_max + 0.3, 0.98 * y_max,
                         "2-sided Mann-Whitney U-tests:\npvalues calculated for datasets: " + str(
                             params['pvalsViolins']) + "\npvalues: " + str(pnucround), fontsize=12, horizontalalignment='right',
                         verticalalignment='top')
            plt.savefig(params['analysisfolder'] + figtype + "_Nuclear_" + color + "_combined_reps" * params['CombineReps'] +".pdf")
            plt.close(fig)

        for figtype in ["Violin", "Boxplot"]:
            fig = plt.figure(figsize=(11.69 / 2., 8.27))
            if figtype == "Violin":
                plt.violinplot(distributions[3 * channel + 2], showmeans=True, showextrema=False)
                plt.xticks(np.arange(1, countdatasets+1), np.array(datasetLabels), rotation=45,
                           horizontalalignment="right")
            elif figtype == "Boxplot":
                sns.boxplot(data=distributions[3 * channel + 2], color=".8")
                plt.xticks(np.arange(0, countdatasets), np.array(datasetLabels), rotation=45,
                           horizontalalignment="right")
            plt.ylabel('Number of nascent RNA at TS')
            plt.xlabel('Dataset')
            y_max = max(np.concatenate(distributions[3 * channel + 2]))
            x_max = len(distributions[3 * channel + 2])
            if yAxisViolinTS != "None": y_max = yAxisViolinTS
            plt.ylim([0, y_max])
            if countdatasets != 1 and len(params['pvalsViolins'][0]) != 0:
                plt.text(x_max + 0.3, 0.98 * y_max,
                         "2-sided Mann-Whitney U-tests:\npvalues calculated for datasets: " + str(
                             params['pvalsViolins']) + "\npvalues: " + str(pTSround), fontsize=12, horizontalalignment='right',
                         verticalalignment='top')
            plt.savefig(params['analysisfolder'] + figtype +"_TS_" + color + "_combined_reps" * params['CombineReps'] + ".pdf")
            plt.close(fig)

        # write pvalues to separate file
        if countdatasets != 1 and len(params['pvalsViolins'][0]) != 0:
            pvalsfile = open(params['analysisfolder'] + "pvalues_" + color + ".txt", "w+")
            pvalsfile.write("p-values cellular violin plots: " + str(pcell))
            pvalsfile.write("\np-values nuclear violin plots: " + str(pnuc))
            pvalsfile.write("\np-values TS violin plots: " + str(pTS))
            pvalsfile.close()


############################################################################################################
#### makes bargraph normalization values
############################################################################################################

def compare_norm_values(params):
    print("Comparing normalization values")
    for channel in params['channelsToAnalyze']:
        color = ('cy3', 'cy5')[channel]
        normVals = []
        datasets = []
        for file in params['fileListIn']:
            paramsexp = get_exp_info(file, params)
            fileName = os.path.join(paramsexp['pathOut'],
                                    f"{paramsexp['date']}_{paramsexp['expName']}_cyto_normalization_value_{color}.txt")
            normVal = np.loadtxt(fileName)
            normVals.append(normVal)
            datasets.append(f"{paramsexp['date']}\n{paramsexp['expName']}")

        normtable = np.c_[range(1, len(params['fileListIn']) + 1), datasets, normVals]
        np.savetxt(params['analysisfolder'] + "Cyto_normalization_values_" + color + ".txt", normtable, fmt="%s", delimiter='\t')

        plt.rcParams.update({'figure.autolayout': True})
        fig = plt.figure(figsize=(11.69 / 2., 8.27))
        plt.bar(range(1, len(params['fileListIn']) + 1), normVals, label=datasets)
        plt.ylim([0, max(normVals) * 1.2])
        plt.xticks(np.arange(1, len(params['fileListIn']) + 1), np.array(datasets), rotation=45, horizontalalignment="right")
        plt.title("Normalizing value from cytoplasmic spots, " + color)
        plt.xlabel('Dataset')
        plt.ylabel('Intensity')
        plt.savefig(params['analysisfolder'] + "Cyto_normalization_values_plot_" + color + ".pdf", format='pdf')
        plt.close(fig)


###############################################################
#### pipeline
###############################################################
TSdistrCy3 = []
TSdistrCy5 = []
CelldistrCy3 = []
CelldistrCy5 = []
NucldistrCy3 = []
NucldistrCy5 = []

def FISH_pipeline(params):
    if not isinstance(params, dict):
        parameter_file = params
        if parameter_file[-3:] == '.py':
            print('Converting py parameter file into yml format')
            misc.convertParamFile2YML(parameter_file)
            parameter_file = parameter_file[:-3]+'.yml'
        if not parameter_file[-4:] == '.yml':
            parameter_file += '.yml'
        params = misc.getParams(parameter_file, __file__.replace('.py', '_parameters_template.yml'),
                                ('folderIn', 'fileListIn', 'outputfolder', 'analysisfolder'))
    else:
        parameter_file = ''

    params = calculate_general_parameters(params)

    for fn, file in enumerate(params['fileListIn']):
        print(f"Analyzing smFISH dataset {fn + 1} out of {len(params['fileListIn'])}: {file}")

        paramsexp = get_exp_info(file, params)
        lfn = get_filelist(file, params)
        paramsexp = get_turret_filters(file, lfn, paramsexp)

        dateTimeObj = datetime.now()
        dateTime = dateTimeObj.strftime("%d-%b-%Y_%Hh%Mm%Ss")
        if os.path.exists(parameter_file):
            shutil.copyfile(os.path.abspath(parameter_file),
                            f"{paramsexp['outfilePart1']}FISH_pipeline_parameters_runtime_{dateTime}.yml")
        else:
            parameter_file = f"{paramsexp['outfilePart1']}FISH_pipeline_parameters_runtime_{dateTime}.yml"
            with open(parameter_file, 'w') as f:
                yaml.dump(params, f)
        shutil.copyfile(os.path.abspath(__file__), f"{paramsexp['outfilePart1']}_FISH_pipeline.py")
        shutil.copyfile(os.path.abspath(_version.__file__), f"{paramsexp['outfilePart1']}_smfish_version.py")

        if params['MaxProj'] == 1:
            max_projection(file, lfn, paramsexp)

        if params['RunCellprofiler'] == 1:
            run_cellprofiler(lfn, params, paramsexp)

        if params['RunFindCells'] == 1:
            find_cells(lfn, params, paramsexp)

        if params['RunOptimizeThresh'] == 1:
            run_optimize_threshold(file, lfn, params, paramsexp)

        if params['RunLocalize'] == 1:
            run_localize(file, lfn, params, paramsexp)

        if params['MakeMaskImage'] == 1:
            make_mask_image(lfn, params, paramsexp)

        if params['MakeMontage'] == 1:
            print("Making montages")
            path = os.path.join(paramsexp['pathOut'], paramsexp['date'])
            make_montage(path, lfn, paramsexp['outfilePart1'] + '_montagePython', 'max')
            make_montage(path, lfn, paramsexp['outfilePart1'] + '_montagePython', 'max_cell_mask', colormap='glasbey')
            make_montage(path, lfn, paramsexp['outfilePart1'] + '_montagePython', 'max_nucleus_mask', colormap='glasbey')

        if params['CellCycle'] == 1:
            CellCycle(lfn, params, paramsexp)

        if params['CalculateTranscriptionSites'] == 1:
            calculate_TS(lfn, params, paramsexp)

        if params['MakeHistograms'] == 1:
            rawdistributions = make_TS_count_histograms(lfn, params, paramsexp)

        if params['CalcSingleCellCorr'] == 1:
            calculate_singlecell_corr(lfn, params, paramsexp)


    if params['CombineReps'] == 1:
        if params['UsePrecalcDistr'] == 1:
            if not rawdistributions:
                print("Recalculating distributions")
                for file in params['fileListIn']:
                    paramsexp = get_exp_info(file, params)
                    lfn = get_filelist(file, params, paramsexp)
                    rawdistributions = make_TS_count_histograms(lfn, params, paramsexp['outfilePart1'])
            else:
                 print("Using pre-analyzed datasets")
        distributionscombined = combine_reps(rawdistributions, params)


    if params['Violins'] == 1:
        if params['CombineReps'] == 1:
            distributions = distributionscombined
        elif params['CombineReps'] == 0:
            if params['UsePrecalcDistr'] == 1 and not rawdistributions:
                for file in params['fileListIn']:
                    paramsexp = get_exp_info(file, params)
                    lfn = get_filelist(file, params, paramsexp)
                    rawdistributions = make_TS_count_histograms(lfn, params, paramsexp['outfilePart1'])
                    distributions = rawdistributions
            else:
                print("Using pre-analyzed datasets")
                distributions = rawdistributions
        make_violin_plots(distributions, params)

    if params['CompareNormValues'] == 1:
        compare_norm_values(params)

    try: rawdistributions
    except NameError: rawdistributions = []
    return params, rawdistributions


def main():
    ipy_debug()
    tm = time()
    if len(sys.argv) < 2:
        if os.path.exists('FISH_pipeline_parameters.yml'):
            parameter_file = 'FISH_pipeline_parameters.yml'
        elif os.path.exists('FISH_pipeline_parameters.py'):
            parameter_file = 'FISH_pipeline_parameters.py'
        else:
            raise FileNotFoundError('Could not find the parameter file.')
        print(('Using ' + parameter_file))
    else:
        parameter_file = sys.argv[1]

    params, rawdistributions = FISH_pipeline(parameter_file)

    if os.path.basename(__file__) in [os.path.basename(i) for i in psutil.Process(os.getpid()).cmdline()]:
        imr.kill_vm()  #stop java used for imread, needed to let python exit
        print('Stopped the java vm used for imread.')

    print('------------------------------------------------')
    print(misc.color('Pipeline finished, took {} seconds.'.format(time()-tm), 'g:b'))


if __name__ == '__main__':
    main()
