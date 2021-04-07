#########################################################################################################
#########################################################################################################
import os
import sys
#sys.path.insert(0, '/Users/lenstratl/Documents/dataAnalysis/src')
from misc import *
from spotDetection_functions import *
from numpy import *
from PIL import ImageFilter
from PIL import Image
from scipy import optimize
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
from tifffile import *
import matplotlib
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel
import FISH_pipeline_parameters
from FISH_pipeline_parameters import *
reload(FISH_pipeline_parameters)
from FISH_pipeline_parameters import *
import scipy.stats
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns


#########################################################################################################
#########################################################################################################

lfileListIn = len(fileListIn)
for i in range(0,lfileListIn):
    if fileListIn[i][-1] != "/": fileListIn[i] += "/"
    fileListIn[i] = folderIn+fileListIn[i]

if outputfolder[-1] != "/": outputfolder += "/"

if useControlForEmptyCells == 1 and channelsToAnalyze == [0,1]:
    if controlColor == 1:
        channelsToAnalyze = [1,0]

#########################################################################################################
### Defining functions to make nuclear, cellular and TS histograms
#########################################################################################################
def histNuc(distr, bins, freqAxis, xCenter, exclBinsCount, color, outName):
    hist_count = histogram(distr, bins = bins-.5, normed = 1)[0];
    hist_count_norm = hist_count[exclBinsCount:]/sum(hist_count[exclBinsCount:])
    bins=bins[exclBinsCount:-1]
    fig = plt.figure()
    plt.bar(bins,hist_count_norm)
    if freqAxis != "None": plt.ylim([0,freqAxis])
    plt.title("Number of RNA in the nucleus, "+color)
    plt.xlabel('Nr RNA/nucleus')
    plt.ylabel('Frequency')
    plt.text(xCenter, max(hist_count_norm), "Mean nuclear count = "+str('{0:.2f}'.format((mean(distr))))+" +/- "+str('{0:.2f}'.format((std(distr)/sqrt(len(distr)))))+"\nNumber of cells analyzed = "+str(len(distr)), fontsize = 12, horizontalalignment='center', verticalalignment='top')
    plt.savefig(outName, format = 'pdf')
    plt.close(fig)

def histCell(distr, bins, freqAxis, xCenter, exclBinsCount, color, outName):
    hist_count = histogram(distr, bins = bins-.5, normed = 1)[0];
    hist_count_norm = hist_count[exclBinsCount:]/sum(hist_count[exclBinsCount:])
    bins=bins[exclBinsCount:-1]
    fig = plt.figure()
    plt.bar(bins,hist_count_norm)
    if freqAxis != "None": plt.ylim([0,freqAxis])
    plt.title("Number of RNA in the cell, "+color)
    plt.xlabel('Nr RNA/cell')
    plt.ylabel('Frequency')
    ##
    plt.text(xCenter, max(hist_count_norm), "Mean cellular count = "+str('{0:.2f}'.format((mean(distr))))+" +/- "+str('{0:.2f}'.format((std(distr)/sqrt(len(distr)))))+"\nNumber of cells analyzed = "+str(len(distr)), fontsize = 12, horizontalalignment='center', verticalalignment='top')
    plt.savefig(outName, format = 'pdf')
    plt.close(fig)

def histTSvals(distr, bins, exclBins):
    h=histogram(distr,bins=bins-.5,normed=1)[0]
    normalize = lambda x: x/sum(x[exclBins:])
    hNorm=h/sum(h[exclBins:])
    return hNorm

def fitBurstModel(histTSvals, bins, exclBins):
    P=lambda x: x[0]**bins/(x[0]+1)**(bins+x[1])*special.binom(x[1]+bins-1,x[1]-1)
    Pnorm=lambda x: (lambda p: p/sum(p[exclBins:]))(P(x))
    [n,aT]=optimize.fmin(lambda x: sum((Pnorm(x)-histTSvals)[exclBins:]**2),[4,160./70])
    Pres=P([n,aT]); Pres*=sum(histTSvals[exclBins:])/sum(Pres[exclBins:]);
    return (Pres, n, aT)

def fitPoissonModel(histTSvals, bins, exclBins):
    PPoiss=lambda x: x**bins/special.gamma(1+bins)
    PPoissnorm=lambda x: (lambda p: p/sum(p[exclBins:]))(PPoiss(x))
    cT=optimize.fmin(lambda x: sum((PPoissnorm(x)-histTSvals)[exclBins:]**2),[.1])
    PPoissres=PPoiss(cT[0]); PPoissres*=sum(histTSvals[exclBins:])/sum(PPoissres[exclBins:]);
    return (PPoissres, cT)

def histTS(distr, hist, bins, FitBurstPar, FitPoissonPar, Pres, n, aT, PPoissres, cT, freqAxis, xCenter, exclBinsCount, color, outName):
    fig = plt.figure()
    plt.bar(bins, hist, label = "FISH")
    if FitBurstPar == 1:
        plt.plot(Pres, color = 'red', label = "Bursting model n = "+str(round(n,3))+" RNAs/burst\naT = "+str(round(aT,3)))
    if FitPoissonPar == 1:
        plt.plot(PPoissres, color = 'green', label = "Poisson model Mean = "+str(round(cT[0],3)))
    plt.legend(loc='center right')
    ymax = 1.1*max(hist)
    if freqAxis != "None": ymax = freqAxis
    plt.ylim([0,ymax])
    plt.title("Number of nascent RNA at the TS, "+color)
    plt.xlabel('Nr nascent RNA at TS')
    plt.ylabel('Frequency')
    plt.text(xCenter, 0.98*ymax, "Mean nr nascent RNA = "+str('{0:.2f}'.format((mean(distr))))+" +/- "+str('{0:.2f}'.format((std(distr)/sqrt(len(distr)))))+"\nNote: Does not include the "+str(range(0,exclBins))+" bins", fontsize = 12, horizontalalignment='center', verticalalignment='top')
    plt.savefig(outName, format = 'pdf')
    plt.close(fig)

def writestatsfile(outFilePart, color, distrCell, distrNuc, exclBinsCount, distrTS, exclBins,includeMultTS):
    statsfile = open(outFilePart+"_statistics_"+color+".txt","w+")
    statsfile.write("Number of cells analyzed: "+str(len(distrCell)))
    statsfile.write("\nCellular count "+color+": "+str('{0:.2f}'.format(mean(distrCell)))+" +/- "+str('{0:.2f}'.format(std(distrCell)/sqrt(len(distrCell)))))
    statsfile.write("\nNuclear count "+color+": "+str('{0:.2f}'.format((mean(distrNuc))))+" +/- "+str('{0:.2f}'.format(std(distrNuc)/sqrt(len(distrNuc)))))
    statsfile.write("\nNote:These counts do not include the "+str(range(0,exclBinsCount))+" bins")
    statsfile.write("\nTxpnsite average "+color+": "+str('{0:.2f}'.format(mean(distrTS)))+" +/- "+str('{0:.2f}'.format(std(distrTS)/sqrt(len(distrTS)))))
    statsfile.write("\nNote:The TS average does not include the "+str(range(0,exclBins))+" bins")
    if includeMultTS == 1:
        statsfile.write("\nThe multiple txpnsite data was calculated with TS threshold "+str(thresholdTS))
    elif includeMultTS == 2:
        statsfile.write("\nThe multiple txpnsite data was calculated with "+str(nrTS2Include)+" TSs per cell")
    statsfile.write("\nNote: errors reported above are standard errors assuming a normal distribution, but the distribution is not normal")
    statsfile.close()

#########################################################################################################
#########################################################################################################

def analysis(pathIn):
#########################################################################################################
#### make file list and output directory
#########################################################################################################

    date = pathIn.split("/")[-3]
    expName = pathIn.split("/")[-2]

    pathOut = outputfolder+date+"_"+expName+"/"
    outFilePart1 = pathOut+date+"_"+expName

    if not os.path.exists(pathOut):
        os.makedirs(pathOut)

    lfn = os.listdir(pathIn)
    if "display_and_comments.txt" in lfn: lfn.remove("display_and_comments.txt")
    if "._display_and_comments.txt" in lfn: lfn.remove("._display_and_comments.txt")
    if ".DS_Store" in lfn: lfn.remove(".DS_Store")
    if "._.DS_Store" in lfn: lfn.remove("._.DS_Store")
    lfn.sort()
    if "excluded" in lfn: lfn.remove("excluded")
    if processAll == 0:
        lfn = lfn[0:nrFilesToAnalyze]
    
    allFiles = os.listdir(pathIn+lfn[0])
    if "metadata.txt" in allFiles: allFiles.remove("metadata.txt")
    splitFiles = [x.split("_") for x in allFiles]
    Turretfilter = []
    for i in range(len(splitFiles)):
        if len(splitFiles[i])>4: Turretfilter.append("_".join(splitFiles[i][2:-1]))
        else: Turretfilter.extend(splitFiles[i][2:-1])
    ChNames = sort(uniqueFromList(Turretfilter))

    os.system("cp FISH_pipeline_parameters.py "+pathOut+date+"_"+expName+"FISH_pipeline_parameters.py")
    os.system("cp FISH_pipeline.py "+pathOut+date+"_"+expName+"FISH_pipeline.py")

###########################################################################################################
##### Does max projection of images
###########################################################################################################

    if MaxProj == 1:
        print("Making maximum intensity projections")
        for fn in range(0, len(lfn)):
            if not os.path.exists(pathOut+date+"_"+lfn[fn]+"_max.tif"):
                imagestack = [[pil.open(pathIn+lfn[fn]+'/img_%09d_%s_%03d.tif'%(0,channel,z)).convert('I') for z in r_[:zSlices]] for channel in ChNames]
                imMaxStack = []
                for c in range(len(imagestack)):
                    imtmp = asarray([array(a.convert("I").getdata()).reshape(a.size) for a in imagestack[c]])
                    imtmp=imtmp*1.
                    imMax = imtmp.max(0)
                    imMaxStack.append(imMax)
                if len(imagestack) ==2: imMaxStack = asarray([imMaxStack[0],imMaxStack[0], imMaxStack[1]]).astype('uint16') # make fake 3 color data
                imsave(pathOut+date+"_"+lfn[fn]+"_max.tif", asarray(imMaxStack).astype('uint16'), imagej=True)


    
    
    ############################################################################################################
    ##### Runs cell profiler to get cell and nucleus masks
    ############################################################################################################

    if RunCellprofiler == 1:
        print("Making nucleus and cell masks")
        os.system("cp "+PipelineCellProfiler+" "+pathOut+PipelineCellProfiler)
        for filenr in range(1, (len(lfn)+1)):
            if not os.path.exists(pathOut+date+"_"+lfn[(filenr-1)]+"_max_cell_mask.tiff"):
                os.system("cellprofiler -p "+re.escape(PipelineCellProfiler)+ " -c -i "+re.escape(pathOut)+" -o "+re.escape(pathOut)+" -f "+str(filenr)+" -l "+str(filenr))

############################################################################################################
#### Runs findCells ipython script to get cell and nucleus masks
############################################################################################################

    if RunFindCells == 1:
        print("Making nucleus and cell masks")
        if channelCellMask == "Cy3":channelCell = 0
        elif channelCellMask == "Cy5":channelCell = 1
        if channelsToAnalyze == [0]:channelNuc = 1
        elif channelsToAnalyze == [0,1]:channelNuc = 2
        for fn in range(0, (len(lfn))):
            if not os.path.exists(pathOut+date+"_"+lfn[(fn)]+"_max_cell_mask.tiff"):
                imagefile = imread(pathOut+date+"_"+lfn[fn]+"_max.tif")
                [cells,nuclei] = findcells(imagefile[channelCell,:,:], imagefile[channelNuc,:,:],cellcolormask=None, ccdist=ccdist, threshold=threshold, thresholdnuc=thresholdnuc, removeborders=removeborders)
                imsave(pathOut+date+"_"+lfn[fn]+"_max_nucleus_mask.tiff", nuclei.astype('uint16'), imagej = True)
                imsave(pathOut+date+"_"+lfn[fn]+"_max_cell_mask.tiff", cells.astype('uint16'), imagej = True)


    ############################################################################################################
    #### runs threshold optimization for Cy3
    ############################################################################################################
    if RunOptimizeThresh == 1:

        par = ["Threshold_cy3="+str(threshold_cy3),"Threshold_cy5="+str(threshold_cy5),"psfPx="+str(psfPx),"MaxDist="+str(maxDist),"minSeparation="+str(minSeparation)]

        for channel in channelsToAnalyze:
            if channel == 0: color = "cy3"
            if channel == 1: color = "cy5"
        
            print "Optimizing Threshold "+color

            tbins = 2

            threshvals = spotslist = spotslistlocalmax = zeros(50)
            threshvals = range(0,40,tbins) + range(40, 120, tbins*2) + range (120,200,tbins *4)
            locMaxValAll = []

            for fn in range(0, len(lfn)):
                if localizeDim == "2D":
                    imagefile = pathOut+date+"_"+lfn[fn]+"_max.tif"   # Input file
                    imtmp = pil.open(imagefile)
                    size = imtmp.size
                    imtmp.seek(channel)                        # 0 for 1st channel, 1 for second channel, 2 for 3rd channel
                    im = array(imtmp.convert("I").getdata()).reshape(size)
                    im=im*1.
                    imBpass=bpass(im,0.75*psfPx,psfPx) # Band-passed image
                
                elif localizeDim == "3D":
                    imagestack = [pil.open(pathIn+lfn[fn]+'/img_%09d_%s_%03d.tif'%(0,ChNames[channel],z)).convert('I') for z in r_[:zSlices]]
                    im = asarray([array(a.convert("I").getdata()).reshape(a.size) for a in imagestack])
                    im=im*1.
                    imBpass=bpass3D(im,1.,psfPx,1.,psfPxZ,zMirror=4) #*66000 # !!!
                
                locMax=rollaxis(array([roll(imBpass,i,-1)  *r_[0,ones(imBpass.shape[-1]  -2),0] for i in r_[-1:2]]).max(0),-1); # copies max values in images to left and right (except outside columns), turns axis
                locMax=rollaxis(array([roll(locMax,i,-1)*r_[0,ones(locMax.shape[-1]-2),0] for i in r_[-1:2]]).max(0),-1); # does the same along second image axis
                if localizeDim == "3D": locMax=rollaxis(array([roll(locMax,i,-1)*r_[0,ones(locMax.shape[-1]-2),0] for i in r_[-1:2]]).max(0),-1); # does the same along third image axis
                locMax=(locMax==imBpass)*1 # finds were in bandpass image the maximum value was, rest will be 0
                locMaxCoo=array(where(locMax)) # finds coordinates
                if localizeDim == "2D": locMaxCoo=tuple(locMaxCoo.T[where((locMaxCoo[0]!=0)*(locMaxCoo[0]!=imBpass.shape[0]-1)*(locMaxCoo[1]!=0)*(locMaxCoo[1]!=imBpass.shape[1]-1))].T)
                elif localizeDim == "3D": locMaxCoo=tuple(locMaxCoo.T[where((locMaxCoo[0]!=0)*(locMaxCoo[0]!=imBpass.shape[0]-1)*(locMaxCoo[1]!=0)*(locMaxCoo[1]!=imBpass.shape[1]-1)*(locMaxCoo[2]!=0)*(locMaxCoo[2]!=imBpass.shape[2]-1))].T)
                locMaxVal=imBpass[locMaxCoo];
                locMaxValAll.extend(locMaxVal)
                for i in range(0,len(threshvals)):
                    thresh = threshvals[i]
                    # nr objects without local max threshold
                    imBinary=(imBpass>thresh)*1.
                    objects=ndimage.find_objects(ndimage.label(imBinary)[0])
                    spotslist[i] = spotslist[i] + len(objects)
                    # nr objects with local max threshold
                    if localizeDim == "2D": yG,xG,valG=tuple(array(locMaxCoo+(locMaxVal,))[:,where(locMaxVal>thresh)])
                    elif localizeDim == "3D": zG,yG,xG,valG=tuple(array(locMaxCoo+(locMaxVal,))[:,where(locMaxVal>thresh)])
                    spotslistlocalmax[i] = spotslistlocalmax[i] + yG.shape[-1]

    #            diffspotslist = abs(diff(spotslist))
    #            suggthrindex = diffspotslist.argmin()
    #            suggthresh = threshvals[suggthrindex]
            locMaxValAll = asarray(locMaxValAll)
            if channel == 0:
                thresh = threshold_cy3
            elif channel == 1:
                thresh = threshold_cy5
            

            fig = plt.figure()
            plt.plot(threshvals,spotslist)
            plt.plot(threshvals,spotslistlocalmax)
#            plt.axvline(x=suggthresh, color = 'red')
#            plt.title("Suggested threshold: "+str(suggthresh))
            plt.xlabel('Threshold '+color)
            plt.ylabel('Number of spots detected')
            plt.yscale('log')
            plt.axvline(thresh,c='r',ls='-',label='Your threshold');
            plt.legend(loc='upper center')
#            plt.ylim(0,(spotslist[-1]*5))
            plt.savefig(pathOut+date+"_"+expName+"_threshold_optimization_"+color+"_"+localizeDim+".pdf", format = 'pdf')
            plt.close(fig)

            savetxt(pathOut+date+"_"+expName+"_thresholds_"+color+".txt",spotslist, delimiter = "\t")


            fig=plt.figure(figsize = (8,3))
            h,x=histogram(log10(maximum(1e-12,locMaxValAll)),bins=r_[-1:log10(locMaxValAll.max()):.01]);
            plt.xscale('log'); plt.xlabel('Pixel intensity');
            plt.yscale('log'); plt.ylabel('Count');
            plt.axvline(median(locMaxValAll),ls='--', label='median'); plt.axvline(mean(locMaxValAll),ls='-',label='mean'); plt.axvline(mean(locMaxValAll)+var(locMaxValAll)**.5,ls=':',label='mean + sd'); plt.axvline(thresh,c='r',ls='-',label='Your threshold');
            plt.text(1, max(h), "Mean + SD = "+str('{0:.2f}'.format(mean(locMaxValAll)+var(locMaxValAll)**.5))+"\nMean = "+str('{0:.2f}'.format(mean(locMaxValAll)))+"\nSD = "+str('{0:.2f}'.format(var(locMaxValAll)**.5))+"\nYour threshold = "+str(thresh), fontsize = 12, horizontalalignment='center', verticalalignment='top')
            plt.plot(10**x[:-1],h+1);
            plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.))
            plt.tight_layout();
            plt.savefig(pathOut+date+"_"+expName+"_pixel_values_local_maxima_"+color+"_"+localizeDim+".pdf", format = 'pdf')
            plt.close(fig)




    ############################################################################################################
    #### runs localize
    ############################################################################################################

    par = ["Threshold_cy3 = "+str(threshold_cy3),"Threshold_cy5 = "+str(threshold_cy5),"Localized in "+localizeDim, "psfPx = "+str(psfPx), "psfPxZ = "+str(psfPxZ),"MaxDist = "+str(maxDist),"minSeparation = "+str(minSeparation), "winSize = "+str(winSize),"winSizeZ = "+str(winSizeZ) ]

    #######################

    if RunLocalize == 1:
        print "Running Localize using "+localizeDim
        for channel in channelsToAnalyze:
            if channel == 0: color = "cy3"
            if channel == 1: color = "cy5"
            
            if localizeDim == "2D": ### localize in 2D
                for fn in range(0, len(lfn)):
                    if (color == "cy3" and not os.path.exists(pathOut+date+"_"+lfn[fn]+"_loc_results_cy3.txt")) or (color == "cy5" and not os.path.exists(pathOut+date+"_"+lfn[fn]+"_loc_results_cy5.txt")):
                        
                        # Input file
                        imagefile = pathOut+date+"_"+lfn[fn]+"_max.tif"
                        
                        # Output files
                        fnTxt=pathOut+date+"_"+lfn[fn]+"_loc_results_"+color+".txt"
                        fnTif=pathOut+date+"_"+lfn[fn]+"_loc_results_"+color+".tif"

                        # open image
                        imtmp = pil.open(imagefile)
                        size = imtmp.size                          # determine image size
                        imtmp.seek(channel)                        # 0 for 1st channel, 1 for second channel, 2 for 3rd channel
                        im = array(imtmp.convert("I").getdata()).reshape(size)
                        im=im*1.
                        imBpass=bpass(im,0.75*psfPx,psfPx)                 # bandpass image to find spots
    
                        if channel == 0:
                            thresh = threshold_cy3
                        elif channel == 1:
                            thresh = threshold_cy5
  
                        # find local maxima
                        locMax=rollaxis(array([roll(imBpass,i,-1)  *r_[0,ones(imBpass.shape[-1]  -2),0] for i in r_[-1:2]]).max(0),-1); # copies max values in images to left and right (except outside columns), turns axis
                        locMax=rollaxis(array([roll(locMax,i,-1)*r_[0,ones(locMax.shape[-1]-2),0] for i in r_[-1:2]]).max(0),-1); # does the same along second image axis
                        locMax=(locMax==imBpass)*1 # finds were in bandpass image the maximum value was, rest will be 0
                        locMaxCoo=array(where(locMax)) # finds coordinates
                        locMaxCoo=tuple(locMaxCoo.T[where((locMaxCoo[0]!=0)*(locMaxCoo[0]!=imBpass.shape[0]-1)*(locMaxCoo[1]!=0)*(locMaxCoo[1]!=imBpass.shape[1]-1))].T)
                        locMaxVal=imBpass[locMaxCoo];
                        yG,xG,valG=tuple(array(locMaxCoo+(locMaxVal,))[:,where(locMaxVal>thresh)])
                        cooGuess=r_[xG,yG].T
                        
                        imBinary=(imBpass>thresh)*1. # Binary image after fixed threshold
#                        elif localizeMethod == "bandpass":
#                            imBinary=(imBpass>thresh)*1. # Binary image after fixed threshold
#                            objects=ndimage.find_objects(ndimage.label(imBinary)[0]) # Find all the connex objects
#                            cooGuess=array([[r_[obj[1]].mean(),r_[obj[0]].mean()] for obj in objects]) # Determine their centers as initial guesses for the spot locations
#
#                        else:
#                            print "no valid localize method found, please select localmax or bandpass in paramater file"

                        # fit each spot with 2D gaussian with tilted plane. The GaussianMaskFit2 function is described in spotDetection_Functions.py. Write spots to parameter fitResults.
                        fitResults=[]
                        for i in range(cooGuess.shape[0]):
                            intensity,coo,tilt=GaussianMaskFit2(im,cooGuess[i],psfPx, winSize = winSize)
                            if intensity!=0: fit = 1
                            elif intensity == 0:
                                intensity,coo,tilt=GaussianMaskFit2(im,cooGuess[i],psfPx, optLoc = 0, winSize = winSize) # if it does not converge, fit at position of initial guess
                                fit = 0
                            # Keep only if it converged close to the inital guess
                            if intensity!=0 and sum((coo-cooGuess[i])**2)<(maxDist*psfPx)**2:
                                # Remove duplicates
                                if sum([sum((coo-a[1:3])**2)<(minSeparation*psfPx)**2 for a in fitResults])==0:
                                    fitResults.append(r_[intensity,coo,tilt,fit])
                    
                        # Save the results in a text file
                        # columns are: integrated intensity, x position, y position, offset of tilted plane, x tilt, y tilt, fit (1 for fit, 0 for forced fit)
                        savetxt(fnTxt,fitResults, delimiter = "\t")
                        
                        # Compute and save image results as a tif file
                        imSpots=im*0.                   # empty image
                        for a in fitResults: imSpots[int(a[2]+.5),int(a[1]+.5)]=1       # add spots to image

                        # for optimzing threshold
                        imGuess = im*0
                        for a in cooGuess: imGuess[int(a[1]), int(a[0])]=1       #

                        #save 6 stack image: original image, imBpass, imBinary (areas above threshold), local maxima, maxima above threshold and fitted spots. These are useful for optimizing the threshold.
                        im2save = asarray([im,imBpass,imBinary,locMax,imGuess,imSpots]).astype('uint16')
                        #im2save = asarray([im,imBpass,imBinary,imGuess,imSpots]).astype('uint16')
                        #im2save = asarray([im,imGuess,imSpots]).astype('uint16')
                        #                  im2save = asarray([((a-a.min())/(a.max()-a.min())) for a in [imBpass,locMax,imGuess,imSpots]]).astype('uint16')
                        #                      im2save2 = concatenate((int16(im).reshape(1,im.shape[0],im.shape[1]),int16(im2save)),axis = 0)
                        imsave(fnTif, im2save, imagej = True)
            
                        print "\n\n*** Found %d spots in '%s' channel. ***\nResults save in '%s' and '%s'."%(len(fitResults),color, fnTxt,fnTif)

                        if channel == 0: par.append("Image "+str(fn)+" in channel "+color+" has threshold: "+str(threshold_cy3)+", found "+str(len(fitResults))+" spots, of which " +str('{0:.0f}'.format(sum(array(fitResults)[:,-1]))) +" ("+str('{0:.2f}'.format(sum(array(fitResults)[:,-1])*100/len(fitResults))) +"%) from fit.")
                        elif channel == 1: par.append("Image "+str(fn)+" in channel "+color+" has threshold: "+str(threshold_cy5)+", found "+str(len(fitResults))+"spots, of which " +str('{0:.0f}'.format(sum(array(fitResults)[:,-1]))) +" ("+str('{0:.2f}'.format(sum(array(fitResults)[:,-1])*100/len(fitResults))) +"%) from fit.")
   
   
            elif localizeDim == "3D": ### localize in 3D
                for fn in range(0, len(lfn)):
                    if (color == "cy3" and not os.path.exists(pathOut+date+"_"+lfn[fn]+"_loc_results_cy3.txt")) or (color == "cy5" and not os.path.exists(pathOut+date+"_"+lfn[fn]+"_loc_results_cy5.txt")):
                        
                        imagestack = [pil.open(pathIn+lfn[fn]+'/img_%09d_%s_%03d.tif'%(0,ChNames[channel],z)).convert('I') for z in r_[:zSlices]]
                        im = asarray([array(a.convert("I").getdata()).reshape(a.size) for a in imagestack])
                        im=im*1.
                        imBpass=bpass3D(im,1.,psfPx,1.,psfPxZ,zMirror=4) #*66000 # !!!
                        
                        # Output files
                        fnTxt=pathOut+date+"_"+lfn[fn]+"_loc_results_"+color+".txt"
                        fnTif=pathOut+date+"_"+lfn[fn]+"_loc_results_"+color+".tif"

                        # find local maxima
                        locMax=rollaxis(array([roll(imBpass,i,-1)  *r_[0,ones(imBpass.shape[-1]  -2),0] for i in r_[-1:2]]).max(0),-1); # copies max values in images to left and right (except outside columns), turns axis
                        locMax=rollaxis(array([roll(locMax,i,-1)*r_[0,ones(locMax.shape[-1]-2),0] for i in r_[-1:2]]).max(0),-1); # does the same along second image axis
                        locMax=rollaxis(array([roll(locMax,i,-1)*r_[0,ones(locMax.shape[-1]-2),0] for i in r_[-1:2]]).max(0),-1); # does the same along third image axis
                        locMax=(locMax==imBpass)*1 # finds were in bandpass image the maximum value was, rest will be 0
                        locMaxCoo=array(where(locMax)) # finds coordinates
                        locMaxCoo=tuple(locMaxCoo.T[where((locMaxCoo[0]!=0)*(locMaxCoo[0]!=imBpass.shape[0]-1)*(locMaxCoo[1]!=0)*(locMaxCoo[1]!=imBpass.shape[1]-1)*(locMaxCoo[2]!=0)*(locMaxCoo[2]!=imBpass.shape[2]-1))].T)
                        locMaxVal=imBpass[locMaxCoo];
                      
                        if channel == 0:
                            thresh = threshold_cy3
                        elif channel == 1:
                            thresh = threshold_cy5
                        
                        zG,yG,xG,valG=tuple(array(locMaxCoo+(locMaxVal,))[:,where(locMaxVal>thresh)])
                        cooGuess=r_[zG,yG,xG].T # for 2D x and y should be swapped!
                        zDepth = (winSizeZ-1)/2
                        
                        intersect = []
                        intersect = [value for value in zG[0] if value in set(range(z)[0:zDepth]+range(z)[-zDepth:])]
     
                        # fit each spot with 3D gaussian with tilted plane. The GaussianMaskFit3D function is described in spotDetection_Functions.py. Write spots to parameter fitResults.
                        fitResults=[]
                        for i in range(cooGuess.shape[0]):
                            intensity,coo,tilt=GaussianMaskFit3D(im,cooGuess[i],psfPx, psfPxZ, winSize = winSize, winSizeZ = winSizeZ)
                            if intensity!=0: fit = 1
                            if intensity == 0:
                                intensity,coo,tilt=GaussianMaskFit3D(im,cooGuess[i],psfPx, psfPxZ, optLoc=0, winSize = winSize, winSizeZ = winSizeZ) # if it does not converge, fit at position of initial guess
                                fit = 0
                            # Keep only if it converged close to the inital guess
                            # Keep only if it converged close to the inital guess
                            if intensity!=0 and sum(((coo-cooGuess[i])/r_[psfPxZ,psfPx,psfPx])**2)<(maxDist*psfPx)**2:
                                # Remove duplicates
                                if sum([sum(((coo-a[1:4])/r_[psfPxZ,psfPx,psfPx])**2)<minSeparation**2 for a in fitResults])==0:
                                    fitResults.append(r_[intensity,coo,tilt,fit])
                    
                        # Save the results in a text file
                        # columns are: integrated intensity, x position, y position, z position, offset of tilted plane, x tilt, y tilt, z tilt, fit (1 = fit, 0 = forced fit)
                        if len(fitResults): fitResults=array(fitResults)[:,[0,3,2,1,4,7,6,5,8]]; # re-order columns
                        savetxt(fnTxt,fitResults, delimiter = "\t")
                        
                        # Compute and save image results as a tif file
                        imSpots=im*0.                   # empty image
                        for a in fitResults: imSpots[int(a[3]+0.5),int(a[2]+.5),int(a[1]+.5)]=1       # add spots to image

                        # for optimzing threshold
                        imGuess = im*0.
                        for a in cooGuess: imGuess[int(a[0]), int(a[1]), int(a[2])]=1       #

                        #save 4 stack image: original image, imBpass, local maxima, maxima above threshold and fitted spots. These are useful for optimizing the threshold.
                        #            im2save = asarray([im,imBpass,locMax,imGuess,imSpots]).astype('uint16')
                        im2save = asarray([im,imGuess,imSpots]).astype('uint16')
                        #                  im2save = asarray([((a-a.min())/(a.max()-a.min())) for a in [imBpass,locMax,imGuess,imSpots]]).astype('uint16')
                        #                  im2save2 = concatenate((int16(im).reshape(1,im.shape[0],im.shape[1]),int16(im2save)),axis = 0)
                        imsave(fnTif, im2save, imagej = True)
            
                        print "\n\n*** Found %d spots in '%s' channel. ***\nResults saved in '%s' and '%s'."%(len(fitResults),color, fnTxt,fnTif)
                        if intersect: print "\n!!!Warning: "+str(len(intersect))+" our of "+ str(len(zG[0]))+"spots were at outer frames of z focus and may not have been fit!!!"

                        if channel == 0: par.append("Image "+str(lfn[fn])+" in channel "+color+" has threshold: "+str(threshold_cy3)+", found "+str(len(fitResults))+" spots, of which " +str('{0:.0f}'.format(sum(fitResults[:,-1]))) +" ("+str('{0:.2f}'.format(sum(fitResults[:,-1])*100/len(fitResults))) +"%) from fit. "+str(len(intersect))+" spots were at outer frames of z focus")
                        elif channel == 1: par.append("Image "+str(lfn[fn])+" in channel "+color+" has threshold: "+str(threshold_cy5)+", found "+str(len(fitResults))+"spots, of which " +str('{0:.0f}'.format(sum(fitResults[:,-1]))) +" ("+str('{0:.2f}'.format(sum(fitResults[:,-1])*100/len(fitResults))) +"%) from fit. "+str(len(intersect))+" spots were at outer frames of z focus")
            

            #save parameters localize
            parameterfile = pathOut+date+"_"+"localize_parameters.txt"
            savetxt(parameterfile,par, delimiter = "\n", fmt="%s")




    ############################################################################################################
    #### makes outlines on cell and nucleus and boxes around localized spots
    ############################################################################################################

    if MakeMaskImage == 1:
        print "Combining mask with image"
        for fn in range(0, len(lfn)):
            if not os.path.exists(pathOut+date+"_"+lfn[fn]+"_mask_cell+nuc+spots.tif"):
                imagefile = pathOut+date+"_"+lfn[fn]+"_max.tif"
                if 0 in channelsToAnalyze: locfile = pathOut+date+"_"+lfn[fn]+"_loc_results_cy3.txt"
                if 1 in channelsToAnalyze: locfile2 = pathOut+date+"_"+lfn[fn]+"_loc_results_cy5.txt"
                nucleusfile = pathOut+date+"_"+lfn[fn]+"_max_nucleus_mask.tiff"
                cellfile = pathOut+date+"_"+lfn[fn]+"_max_cell_mask.tiff"
                
                
                image = pil.open(imagefile)
                im = array(image.getdata()).reshape(image.size[::-1])
                
                if 0 in channelsToAnalyze:
                    data = loadtxt(locfile)
                    if localizeDim == "3D":
                        data = data.reshape((-1,9))[:,:7]
                    Pos=im*0
                    if data.shape[0] != 0:
                        Loc_xy=(data.reshape((-1,7))[:,1:3]+.5).astype(int)
                        for spot in range(0,data.reshape((-1,7)).shape[0]): Pos[Loc_xy[spot,1],Loc_xy[spot,0]]=1
                
                if 1 in channelsToAnalyze:
                    data2 = loadtxt(locfile2)
                    if localizeDim == "3D":
                        data2 = data2.reshape((-1,9))[:,:7]
                    Pos2=im*0
                    if data2.shape[0] != 0:
                        Loc_xy2=(data2.reshape((-1,7))[:,1:3]+.5).astype(int)
                        for spot in range(0,data2.reshape((-1,7)).shape[0]): Pos2[Loc_xy2[spot,1],Loc_xy2[spot,0]]=1
                
                nucleus = pil.open(nucleusfile)
                masknuc = nucleus.convert("L").filter(ImageFilter.FIND_EDGES)
                maskarraynuc = array(masknuc) >= 1
                
                cell = pil.open(cellfile)
                maskcell = cell.convert("L").filter(ImageFilter.FIND_EDGES)
                maskarraycell = array(maskcell) >= 1
                
                maskarraytot = maskarraycell + maskarraynuc
                maskarraytot.dtype = uint8
                maskarraytot *= 60
                

                res = []
                for n in range(3):
                    image.seek(n)
                    im1 = array(image.convert('I').convert('F'))
                    res.append(int16(im1))
                
                res.append(int16(maskarraytot));
                if 0 in channelsToAnalyze: res.append(int16(Pos))
                if 1 in channelsToAnalyze: res.append(int16(Pos2))
                
                outfileMask = pathOut+date+"_"+lfn[fn]+"_mask_cell+nuc+spots.tif"
                imsave(outfileMask,asarray(res))
                
                os.system("/DATA/lenstra_lab/Fiji.app/ImageJ-linux64"+""" --no-splash --headless -eval 'open("%s"); run("Stack to Hyperstack...", "order=xyczt(default) channels=%d slices=1 frames=1 display=Composite"); run("Save", "%s")'&"""%(re.escape(outfileMask),len(res),re.escape(outfileMask)))




    ############################################################################################################
    #### calculate transcription sites
    ############################################################################################################
    if CalculateTranscriptionSites == 1:
        print "Calculating transcription sites"

        imageTSCy3= []
        imageTSCy5= []
    
        ##### create list of dictionaries (frames) with 0 for each nucleus
        ld_zero = [{} for i in range(0,len(lfn))]
        ld_zero5 = [{} for i in range(0,len(lfn))]
        
        
        for fn in range(0,len(lfn)):
            if MaskTS == "nucleus" : nucleusfile = pathOut+date+"_"+lfn[fn]+"_max_nucleus_mask.tiff"
            elif MaskTS == "cell" : nucleusfile = pathOut+date+"_"+lfn[fn]+"_max_cell_mask.tiff"
            nucleus = pil.open(nucleusfile)
            nuclist = unique(array(nucleus.getdata()).reshape(nucleus.size[::-1]))
            for nuc in nuclist[1:]:
                ld_zero[fn][nuc] = 0
                ld_zero5[fn][nuc] = array([0,0,0,0,0])

        ##### loop over frames
        for channel in channelsToAnalyze:
            if channel == 0: color = "cy3"
            if channel == 1: color = "cy5"
            
            ldCountNuc = copy2.deepcopy(ld_zero) # parameter file for number of spots in nucleus
            ldCountCell = copy2.deepcopy(ld_zero) # parameter file for number of spots in cytoplasm (excluding nucleus)
            ldCountBoth = copy2.deepcopy(ld_zero) # parameter file for number of spots in entire cell (both cytoplasm and nucleus)
            ld = copy2.deepcopy(ld_zero) # intensity for brightest TS
            ldPos = copy2.deepcopy(ld_zero5) # for brightest TS with position information
            ldNucTS = copy2.deepcopy(ld_zero) # for all nucl spots (filter later for TS)
            ldNucTSPos = copy2.deepcopy(ld_zero5) # for all nucl spots (filter later for TS) wiht position information

            dataNucleus = []    # parameter file for spot intensities in nucleus
            dataCell = []       # parameter file for spot intensities in cytoplasm (excluding nucleus)

            for fn in range(0, len(lfn)):
                sys.stdout.write("\rComputing TS %d/%d ..."%(fn+1,int(len(lfn))+1)); sys.stdout.flush()
                locfile = pathOut+date+"_"+lfn[fn]+"_loc_results_"+color+".txt"
                nucleusfile = pathOut+date+"_"+lfn[fn]+"_max_nucleus_mask.tiff"
                cellfile = pathOut+date+"_"+lfn[fn]+"_max_cell_mask.tiff"

                ##### load data
                data = loadtxt(locfile)
                if localizeDim == "3D":
                    data = data.reshape((-1,9))[:,r_[:5,-1]]
                if localizeDim == "2D":
                    data = data.reshape((-1,7))[:,r_[:5,-1]]
                nucleus = pil.open(nucleusfile)
                cell = pil.open(cellfile)
                
                nucleusArray= array(nucleus.getdata()).reshape(nucleus.size[::-1])
                cellArray = array(cell.getdata()).reshape(cell.size[::-1])
                
                # select for cell and nucleus spots
                for row in data.reshape((-1,6)):
                    x = int(round(row[1]))
                    y = int(round(row[2]))
                    nucnr = nucleusArray[y, x]
                    cellnr = cellArray[y, x]
                    if nucnr > 0 and nucnr != 32768:
                        ldCountNuc[fn][nucnr] +=1                               ### add count to ldCountNuc
                        ldCountBoth[fn][nucnr] +=1                              ### add count to ldCountBoth
                        if row[-1] > 0:                                         ### check if fit was not forced
                            dataNucleus.append(row[0])                          ### add intensity to dataNucleus
                        if includeMultTS > 0:
                            ldNucTS[fn][nucnr] = vstack((ldNucTS[fn][nucnr],row[0]))                    ### add intensity nuclear spot
                            ldNucTSPos[fn][nucnr] = vstack((ldNucTSPos[fn][nucnr], row[[0,1,2,3,5]]))   ### add intensity and x,y location, fit
                        if row[0] > ld[fn][nucnr]:
                            ld[fn][nucnr] = row[0]                          ### replace intensity in ld if spot has higher intensity
                            ldPos[fn][nucnr] =  row[[0,1,2,3,5]]              ### replace intensity and x,y location, fit
                    elif cellnr  > 0 and cellnr != 32768:                       ### if spot is in cytoplasm but not in nucleus
                        ldCountCell[fn][cellnr] +=1                             ### add count to ldCountCell
                        ldCountBoth[fn][cellnr] +=1                             ### add count to ldCountBoth
                        if row[-1] > 0:                                         ### check if fit is not forced
                            dataCell.append(row[0])                             ### add intensity to dataCell
                        if MaskTS == "cell":
                            if includeMultTS > 0:
                                ldNucTS[fn][cellnr] = vstack((ldNucTS[fn][cellnr],row[0]))                      ### add intensity cell spot
                                ldNucTSPos[fn][cellnr] = vstack((ldNucTSPos[fn][cellnr],row[[0,1,2,3,5]]))     ### add intensity and x,y location, fit
                            if row[0] > ld[fn][cellnr]:
                                ld[fn][cellnr] = row[0]                         ### replace intensity in ld if spot has higher intensit
                                ldPos[fn][cellnr] = row[[0,1,2,3,5]]              ### replace intensity and x,y location
                    else:
                        continue

            countUnfiltered = []
            for d in range(0,int(len(ldCountBoth))): countUnfiltered.extend(ldCountBoth[d].values())
            
            dataTxpnSiteUnfiltered = []
            for d in range(0,int(len(ld))): dataTxpnSiteUnfiltered.extend(ld[d].values())
            savetxt(outFilePart1+"_CountNuc+Cytoplasm_"+color+"_nofilter.txt",countUnfiltered, fmt='%.8e', delimiter='\t')
            savetxt(outFilePart1+"_intensity_brightest_txpnsite_"+color+"_nofilter.txt", dataTxpnSiteUnfiltered, fmt='%.8e', delimiter='\t')
            
            dataTxpnSiteUnfilteredPos = []
            for d in range(0,int(len(ldPos))): dataTxpnSiteUnfilteredPos.extend(ldPos[d].values())
            savetxt(outFilePart1+"_intensity_brightest_txpnsite_xyfit_"+color+"_nofilter.txt", dataTxpnSiteUnfilteredPos, fmt='%.8e', delimiter='\t')
     
            ### filter spots
            
            ### remove unfitted spots from TS distribution, (only for single TS, for multiple TS see code below)
            for fn in range(0, len(lfn)):
                for nuc in ldPos[fn].keys():
                    if ldPos[fn][nuc][-1] == 0:
                        ld[fn][nuc] = []
                        ldPos[fn][nuc] = []
#                    if includeMultTS > 0:
#                        for TS in range(0,r_[ldNucTS[fn][nuc]].shape[0]):
#                            if ldNucTSPos[fn][nuc][TS][3] == 0:
#                                ldNucTS[fn][nuc][TS] = NaN
#                                ldNucTSPos[fn][nuc][TS] = [NaN,NaN,NaN,NaN]
     
     
            ### remove empty cells from TS list and count list (use either own channel or other channel)
            if ExclEmpty == 1:
                for fn in range(0, len(lfn)):
                    for nuc in ld[fn].keys():
                        if ldCountBoth[fn][nuc] == 0:
                            ld[fn][nuc] = []
                            ldPos[fn][nuc] = []
                            if includeMultTS > 0:
                                ldNucTS[fn][nuc] = []
                                ldNucTSPos[fn][nuc] = []
                            ldCountNuc[fn][nuc] = []
                            ldCountBoth[fn][nuc] = []
                            ldCountCell[fn][nuc] = []
            elif useControlForEmptyCells == 1:
                if channel == controlColor:
                    CountControl = ldCountBoth
                for fn in range(0, len(lfn)):
                    for nuc in ldNucTS[fn].keys():
                        if CountControl[fn][nuc] == 0:
                            ld[fn][nuc] = []
                            ldPos[fn][nuc] = []
                            if includeMultTS > 0:
                                ldNucTS[fn][nuc] = []
                                ldNucTSPos[fn][nuc] = []
                            ldCountNuc[fn][nuc] = []
                            ldCountBoth[fn][nuc] = []
                            ldCountCell[fn][nuc] = []




            #### define normaling spots
            normSpots = []
            if AreaNorm == "cyto":
                normSpots = dataCell
            elif AreaNorm == "nucleus":
                normSpots = dataNucleus
            elif AreaNorm == "cell":
                for d in range(0,len(dataCell)): normSpots.append(d)
                for d in range(0,len(dataNucleus)): normSpots.append(d)
            
            ### save cytoplasmic intensity distribution
            dataCellMedian = median(array(normSpots))
            maxbin = 8 * dataCellMedian
            nobins = 100
            binvalues = zeros(nobins)
            for i in range(0,nobins):
                binvalues[i] = maxbin * i/nobins
            histo = histogram(normSpots,binvalues,normed = True)
            bincenters = zeros(len(binvalues)-1)
            for i in range(0,nobins-1):
                bincenters[i] = 0.5*(binvalues[i]+binvalues[i+1])

            mod = GaussianModel()
            pars = mod.guess(histo[0],x=bincenters)
            out = mod.fit(histo[0],pars,x=bincenters)
            normvalue = out.params['center'].value
            
            if useFit == 0: cal = "median"
            elif useFit == 1: cal = "fit"

            fig = plt.figure()
            plt.bar(bincenters,histo[0],width = 0.5*binvalues[1])
            plt.plot(bincenters,out.best_fit)
            plt.axvline(x=normvalue, color = 'red')
            plt.title("Fit="+str('{0:.1f}'.format(normvalue))+", median="+str('{0:.1f}'.format(dataCellMedian))+", calibration used: "+cal)
            plt.xlabel('Cytoplasmic spot intensity')
            plt.ylabel('Frequency')
            plt.xlim([0,8*dataCellMedian])
            plt.savefig(outFilePart1+"_cyto_intensity_distribution_"+color+".pdf", format = 'pdf')
            plt.close(fig)

            if useFit == 0:
                normCyt = dataCellMedian
            else:
                normCyt = normvalue
            
            savetxt(outFilePart1+"_cyto_normalization_value_"+color+".txt", array([normCyt]), fmt='%.8e', delimiter='\t')

            ### normalize Txpn site
            dataNucleus = array(dataNucleus)

            dataTxpnSite = []

            for d in range(0,int(len(ld))): dataTxpnSite.extend(ld[d].values())
            dataTxpnSite = [x for x in dataTxpnSite if isinstance(x,(int,float))]

            dataNuclNorm = dataNucleus / normCyt
            dataTxpnSiteNorm = dataTxpnSite / normCyt
            dataTxpnSiteNormUnfiltered = dataTxpnSiteUnfiltered / normCyt

            ### filter and define TS for multiple TSs
            if includeMultTS > 0:
                dataTxpnSiteMult = []
                countTS = []
                dataTxpnSiteMultNorm = []
                for fn in range(0,len(ldNucTS)):
                    for nuc in ldNucTS[fn].keys():
                        TSlist = array(ldNucTS[fn][nuc]).reshape(-1,1).ravel().tolist()
                        if len(TSlist)> 1: TSlist.sort()
                        if ldNucTSPos[fn][nuc] != []:
                            ldNucTSPos[fn][nuc].view('i8,i8,i8,i8,i8').sort(order=['f0'], axis=0)
                            if includeMultTS == 2:
                                if nrTS2Include > len(TSlist):
                                    for addzero in range(len(TSlist),nrTS2Include): TSlist.append(0.)
                                    for addzero2 in range(ldNucTSPos[fn][nuc].reshape(-1,5).shape[0],nrTS2Include): ldNucTSPos[fn][nuc] = vstack((array([0,0,0,0,0]),ldNucTSPos[fn][nuc]))
                                fittedSpots = []
                                for fittedSpot in range(len(TSlist)-nrTS2Include,len(TSlist)):
                                    if ldNucTSPos[fn][nuc][fittedSpot][-1] != 0:
                                        dataTxpnSiteMult.append(TSlist[fittedSpot])
                                        fittedSpots.append(fittedSpot)
                                    elif (ldNucTSPos[fn][nuc][fittedSpot][-1] == 0) & (ldNucTSPos[fn][nuc][fittedSpot,] == array([0,0,0,0,0])).all():
                                        dataTxpnSiteMult.append(TSlist[fittedSpot])
                                        fittedSpots.append(fittedSpot)
                                if len(fittedSpots) > 0:  ldNucTSPos[fn][nuc] = ldNucTSPos[fn][nuc][fittedSpots,]
                            elif includeMultTS == 1:
                                nrTS = 0
                                toRemove = []
                                for ts in range(0, len(TSlist)):
                                    if TSlist[ts] > normCyt * thresholdTS:
                                        nrTS += 1
                                        if ldNucTSPos[fn][nuc].reshape(-1,5)[ts][-1] != 0:
                                            dataTxpnSiteMult.append(TSlist[ts])
                                toRemove =(ldNucTSPos[fn][nuc].reshape(-1,5)[:,0]< (normCyt * thresholdTS)) | (ldNucTSPos[fn][nuc].reshape(-1,5)[:,-1] == 0)
                                toKeep = [not i for i in toRemove]
                                ldNucTSPos[fn][nuc] = ldNucTSPos[fn][nuc].reshape(-1,5)[toKeep]
                                countTS.append(nrTS)
                dataTxpnSiteMultNorm =  dataTxpnSiteMult/ normCyt
            

            ### convert dictionaries to arrays/lists
    
            countNucleus = []
            countCell = []
            countBoth =[]
        
            for d in range(0,int(len(ldCountNuc))): countNucleus.extend(ldCountNuc[d].values())
            for d in range(0,int(len(ldCountCell))): countCell.extend(ldCountCell[d].values())
            for d in range(0,int(len(ldCountBoth))): countBoth.extend(ldCountBoth[d].values())

            countNucleus = filter(lambda t: t!=[],countNucleus)
            countCell = filter(lambda t: t!=[],countCell)
            countBoth = filter(lambda t: t!=[],countBoth)

            #### save results intensity spots/txpnsite
            savetxt(outFilePart1+"_intensity_nuclear_spots_"+color+".txt", dataNucleus, fmt='%.8e', delimiter='\t')
            savetxt(outFilePart1+"_normalized_intensity_nuclear_spots_"+color+".txt", dataNuclNorm, fmt='%.8e', delimiter='\t')
            savetxt(outFilePart1+"_intensity_cytoplasmic_spots_"+color+".txt", dataCell, fmt='%.8e', delimiter='\t')
            savetxt(outFilePart1+"_intensity_brightest_txpnsite_"+color+".txt", dataTxpnSite, fmt='%.8e', delimiter='\t')
            savetxt(outFilePart1+"_normalized_intensity_brightest_txpnsite_"+color+".txt", dataTxpnSiteNorm, fmt='%.8e', delimiter='\t')
            savetxt(outFilePart1+"_normalized_intensity_brightest_txpnsite_"+color+"_nofilter.txt", dataTxpnSiteNormUnfiltered, fmt='%.8e', delimiter='\t')
            if includeMultTS > 0:
                savetxt(outFilePart1+"_intensity_txpnsite_"+color+"_mult_method"+str(includeMultTS)+".txt", dataTxpnSiteMult, fmt='%.8e', delimiter='\t')
                savetxt(outFilePart1+"_normalized_intensity_txpnsite_"+color+"_mult_method"+str(includeMultTS)+".txt", dataTxpnSiteMultNorm, fmt='%.8e', delimiter='\t')
                save(outFilePart1+"_normalized_intensity_txpnsite_"+color+"_mult_method"+str(includeMultTS)+".npy", dataTxpnSiteMultNorm)
                if includeMultTS == 1:
                    savetxt(outFilePart1+"_count_nrTS"+color+"_threshold_"+str(thresholdTS)+".txt", countTS, fmt='%.8e', delimiter='\t')
                    save(outFilePart1+"_count_nrTS"+color+"_threshold_"+str(thresholdTS)+".npy", countTS)

            save(outFilePart1+"_CountNuc_"+color+".npy",countNucleus)
            save(outFilePart1+"_CountCytoplasm_"+color+".npy",countCell)
            save(outFilePart1+"_CountNuc+Cytoplasm_"+color+".npy",countBoth)
            save(outFilePart1+"_normalized_intensity_txpnsite_"+color+".npy", dataTxpnSiteNorm)


            ### calculate mask TS
            nucleusfile = pathOut+date+"_"+lfn[fn]+"_max_nucleus_mask.tiff"
            size = pil.open(nucleusfile).size

            for fn in range(0,len(lfn)):
                imageTS = zeros(size)
                squareStamp=[r_[-5,-5,-5,-5,-5,-5,-5,-5,-4,-3,-2,2,3,4,5,5,5,5,5,5,5,5,4,3,2,-2,-3,-4],r_[-5,-4,-3,-2,2,3,4,5,5,5,5,5,5,5,5,4,3,2,-2,-3,-4,-5,-5,-5,-5,-5,-5,-5]] # add squares around the spots in seperate channel
                if includeMultTS == 0:
                    for nuc in ldPos[fn].keys():
                        if ldPos[fn][nuc] != []:
                            if ldPos[fn][nuc][0] != 0:
                                x = int(round(ldPos[fn][nuc][1]))
                                y = int(round(ldPos[fn][nuc][2]))
                                xx = (x+squareStamp[0])[logical_and((y+squareStamp[1])<2048, (y+squareStamp[1])>=0)]
                                yy = (y+squareStamp[1])[logical_and((y+squareStamp[1])<2048, (y+squareStamp[1])>=0)]
                                xx2 = xx[logical_and(xx<2048, xx>=0)]
                                yy2 = yy[logical_and(xx<2048, xx>=0)]
                                imageTS[ yy2, xx2]=1
                #                     imageTS[ y + squareStamp[1],  x + squareStamp[0]]=1
                if includeMultTS > 0 :
                    for nuc in ldNucTSPos[fn].keys():
                        if ldNucTSPos[fn][nuc] != []:
                            for ts in range(0,ldNucTSPos[fn][nuc].shape[0]):
                                x = int(round(ldNucTSPos[fn][nuc][ts][1]))
                                y = int(round(ldNucTSPos[fn][nuc][ts][2]))
                                xx = (x+squareStamp[0])[logical_and((y+squareStamp[1])<2048, (y+squareStamp[1])>=0)]
                                yy = (y+squareStamp[1])[logical_and((y+squareStamp[1])<2048, (y+squareStamp[1])>=0)]
                                xx2 = xx[logical_and(xx<2048, xx>=0)]
                                yy2 = yy[logical_and(xx<2048, xx>=0)]
                                if x != 0 and y != 0:
                                    imageTS[ yy2,  xx2]=1

                if channel == 0:
                   imageTSCy3.append(imageTS)
                if channel == 1:
                   imageTSCy5.append(imageTS)


        ### add TS to mask file
        for fn in range(0,len(lfn)):
            if not os.path.exists(pathOut+date+"_"+lfn[fn]+"+TS.tif"):
                sys.stdout.write("\rAdding transcription sites to mask image %d/%d ..."%(fn+1,len(lfn))); sys.stdout.flush()
                file = pathOut+date+"_"+lfn[fn]+"_max.tif"
                maskIm = pil.open(file)
        
                maskImNew = []
                for stack in channelsToAnalyze:
                    maskIm.seek(stack)
                    im = array(maskIm.convert('I').convert('F'))
                    maskImNew.append(int16(im))
            
                if 0 in channelsToAnalyze:
                    maskImNew.append(int16(imageTSCy3[fn]))
                if 1 in channelsToAnalyze:
                    maskImNew.append(int16(imageTSCy5[fn]))

                outfileMask = pathOut+date+"_"+lfn[fn]+"+TS.tif"
                imsave(outfileMask,asarray(maskImNew))

                os.system("/DATA/lenstra_lab/Fiji.app/ImageJ-linux64"+""" --no-splash --headless -eval 'open("%s"); run("Stack to Hyperstack...", "order=xyczt(default) channels=%d slices=1 frames=1 display=Composite"); run("Save", "%s")'&"""%(re.escape(outfileMask),len(maskImNew),re.escape(outfileMask)))


    ############################################################################################################
    #### makes histograms
    ############################################################################################################

    if MakeHistograms == 1:
        print "Making histograms"
        for channel in channelsToAnalyze:

            if channel == 0: color = "cy3"
            if channel == 1: color = "cy5"
        
            ### load data
            if includeMultTS == 0:
                dataTS = load(outFilePart1+"_normalized_intensity_txpnsite_"+color+".npy")
                method = "1 TS per cell"
            elif includeMultTS == 1:
                dataTS = load(outFilePart1+"_normalized_intensity_txpnsite_"+color+"_mult_method"+str(includeMultTS)+".npy")
                method = "TS threshold "+str(thresholdTS)
                countTS = load(outFilePart1+"_count_nrTS"+color+"_threshold_"+str(thresholdTS)+".npy")
            elif includeMultTS == 2:
                dataTS = load(outFilePart1+"_normalized_intensity_txpnsite_"+color+"_mult_method"+str(includeMultTS)+".npy")
                method = str(nrTS2Include) + " TSs per cell"

            ### remove zeros
            dataTxpnSiteNormNonZero = [i for i in dataTS if i > (exclBins-0.5)] #remove TS data of excluded bins

            if channel == 0:
                if max(dataTxpnSiteNormNonZero)>nrbins_cy3_ts: print "Warning! Cy3 TS fewer bins than max value"
            if channel == 1:
                if max(dataTxpnSiteNormNonZero)>nrbins_cy5_ts: print "Warning! Cy5 TS fewer bins than max value"
        
            #### make histogram of transcription site. Negative binomial distribution as from [Raj 2006 PLoSBio] and [     ].
        
        
            if channel == 0: bins=r_[:nrbins_cy3_ts]; xCenterTS = nrbins_cy3_ts/2; freqAxisTS = freqAxisTS_cy3
            elif channel == 1: bins=r_[:nrbins_cy5_ts]; xCenterTS = nrbins_cy5_ts/2; freqAxisTS = freqAxisTS_cy5
           
            hNorm = histTSvals(dataTxpnSiteNormNonZero,bins,exclBins)
            bins=bins[:-1]
            if CalcFracOff == 1 and exclBins == 0 and exclBinsCount == 0 and ExclEmpty == 1: savetxt(outFilePart1+"_fraction_off_"+color+".txt", [hNorm[0]], delimiter='\t')

            
            dataTSnonZero = []
            for i in range(0,len(dataTS)):
                if dataTS[i]>(exclBins-1): dataTSnonZero.append(dataTS[i])
            
            if channel == 0: TSdistrCy3.append(dataTSnonZero)
            if channel == 1: TSdistrCy5.append(dataTSnonZero)
            
            #### bursting model
            if FitBurst == 1: [Pres, n, aT] = fitBurstModel(hNorm, bins, exclBins)
            else: [Pres, n, aT] = [0, 0, 0]

            if FitPoisson == 1: [PPoissres, cT] = fitPoissonModel(hNorm, bins, exclBins)
            else: [PPoissres, cT] = [0, 0]

            #### plot figure TS intensity distribution (and TS count)
            histTS(dataTSnonZero, hNorm[exclBins:], bins[exclBins:], FitBurst, FitPoisson, Pres, n, aT, PPoissres, cT, freqAxisTS, xCenterTS, exclBins, color, outFilePart1+"_txpnsite_histogram_"+color+"_method="+method.replace(" ", "_")+".pdf")

            ### plot nr of TS per nucleus
            if includeMultTS == 1:
                if channel == 0: bins = r_[:nrbins_cy3_ts]; xCenter = nrbins_cy3_ts/2
                elif channel == 1: bins = r_[:nrbins_cy5_ts]; xCenter = nrbins_cy5_ts/2
                hist_TS_count = histogram(countTS, bins = bins-.5, normed = 1)[0];
                bins=bins[:-1]
                fig = plt.figure()
                plt.bar(bins,hist_TS_count)
                if freqAxisTS != "None": plt.ylim([0,freqAxisTS])
                plt.title("Number of TS in the nucleus, "+color)
                plt.xlabel('Nr TS/nucleus')
                plt.ylabel('Frequency')
                plt.text(xCenter, max(hist_TS_count), "Mean TS count = "+str('{0:.2f}'.format(mean(countTS)))+" +/- "+str('{0:.2f}'.format((std(countTS)/sqrt(len(countTS)))))+"\nNumber of cells analyzed = "+str(len(countTS))+"\nThreshold TS = "+str(thresholdTS), fontsize = 12, horizontalalignment='center', verticalalignment='top')
                plt.savefig(outFilePart1+"_TS_count_"+color+"_threshold="+str(thresholdTS)+".pdf", format = 'pdf')
                plt.close(fig)
            
     
                #matplotlib.pyplot.violinplot(dataset, positions=None, vert=True, widths=0.5, showmeans=False, showextrema=True, showmedians=False, points=100, bw_method=None, *, data=None)[source]

                
            ### counting spots
  
            countNucleus = load(outFilePart1+"_CountNuc_"+color+".npy")
            countCell = load(outFilePart1+"_CountCytoplasm_"+color+".npy")
            countBoth = load(outFilePart1+"_CountNuc+Cytoplasm_"+color+".npy")


            if channel == 0:
                if max(countNucleus)>nrbins_cy3: print "Warning! Cy3 Nuclear spots - fewer bins than max value"
                if max(countBoth)>nrbins_cy3: print "Warning! Cy3 Cellular spots - fewer bins than max value"
            if channel == 1:
                if max(countNucleus)>nrbins_cy5: print "Warning! Cy5 Nuclear spots - fewer bins than max value"
                if max(countBoth)>nrbins_cy5: print "Warning! Cy5 Cellular spots - fewer bins than max value"

            if channel == 0: bins = r_[:nrbins_cy3]; xCenter = nrbins_cy3/2; freqAxisCell = freqAxisCell_cy3; freqAxisNuc = freqAxisNuc_cy3
            elif channel == 1: bins = r_[:nrbins_cy5]; xCenter = nrbins_cy5/2; freqAxisCell = freqAxisCell_cy5; freqAxisNuc = freqAxisNuc_cy3
    
            totalCellsCounted = len(countBoth)
            countBothNonZero = [i for i in countBoth if i > (exclBinsCount-1)]
            totalCellsCountedNonzero = len(countBothNonZero)
            countNucleusNonZero = [i for i in countNucleus if i > (exclBinsCount-1)]
            totalNucleusCountedNonzero = len(countNucleusNonZero)

            if channel == 0:
                CelldistrCy3.append(countBothNonZero)
                NucldistrCy3.append(countNucleusNonZero)
            if channel == 1:
                CelldistrCy5.append(countBothNonZero)
                NucldistrCy5.append(countNucleusNonZero)

            histNuc(countNucleusNonZero, bins, freqAxisNuc, xCenter, exclBinsCount, color, outFilePart1+"_nuclear_count_"+color+".pdf")
            histCell(countBothNonZero, bins, freqAxisCell, xCenter, exclBinsCount, color, outFilePart1+"_cell_count_"+color+".pdf")

            savetxt(outFilePart1+"_cellular_count_"+color+".txt", countBoth, fmt='%i', delimiter='\t')
            savetxt(outFilePart1+"_nuclear_count_"+color+".txt", countNucleus, fmt='%i', delimiter='\t')

            writestatsfile(outFilePart1, color, countBothNonZero, countNucleusNonZero, exclBinsCount, dataTSnonZero, exclBins, includeMultTS)



    ############################################################################################################
    #### calculates correlations between channels and between cell and nascent RNA number
    ############################################################################################################

    if CalcSingleCellCorr == 1:
        print "Calculating correlation"
        cors = []
        
        ### Scatter plot between nr of nascent RNA and RNA count Cy3
        if 0 in channelsToAnalyze:
            CountCy3 = loadtxt(outFilePart1+"_CountNuc+Cytoplasm_cy3_nofilter.txt")
            TSCy3 = loadtxt(outFilePart1+"_normalized_intensity_brightest_txpnsite_cy3_nofilter.txt")
            TSCy3Pos = loadtxt(outFilePart1+"_intensity_brightest_txpnsite_xyfit_cy3_nofilter.txt")
            if ExclEmpty == 1:
                EmptyCy3 = list(where(CountCy3 ==0 )[0])
            else:
                EmptyCy3 = []
            NoFitCy3 = list(where(TSCy3Pos[:,-1] == 0)[0])
            NoFit_EmptyCy3 = uniqueFromList(EmptyCy3 + NoFitCy3)
            NotEmptyCy3 = [x for x in range(0,len(CountCy3)) if x not in NoFit_EmptyCy3]
            corTScountCy3  = pearsonr(CountCy3[NotEmptyCy3], TSCy3[NotEmptyCy3])
            fig = plt.figure()
            plt.scatter(CountCy3[NotEmptyCy3], TSCy3[NotEmptyCy3], s = 2)
            plt.title("Cy3, correlation = "+str('{0:.3f}'.format(corTScountCy3[0])))
            plt.xlabel('Nr RNA/cell, Cy3')
            plt.ylabel('Nr of nascent RNA at TS')
            plt.savefig(outFilePart1+"_correlation_count_with_TSintensity_cy3.pdf", format = 'pdf')
            plt.close(fig)
            cors.append("Correlation Cy3 TS intensity and cell count = "+str(corTScountCy3[0])+" , p-value = "+ str('{:.2E}'.format(corTScountCy3[1])))
            cors.append(str(len(NotEmptyCy3)) +" cells were analyzed, "+str(len(NoFit_EmptyCy3))+" cells were filtered")
            cors.append(str(len(EmptyCy3))+" cells removed because they were empty in Cy3 channel")
            cors.append(str(len(NoFitCy3))+" cells were removed because the Cy3 TS were not fit")
            
        ### Scatter plot between nr of nascent RNA and RNA count Cy5
        if 1 in channelsToAnalyze:
            CountCy5 = loadtxt(outFilePart1+"_CountNuc+Cytoplasm_cy5_nofilter.txt")
            TSCy5 = loadtxt(outFilePart1+"_normalized_intensity_brightest_txpnsite_cy5_nofilter.txt")
            TSCy5Pos = loadtxt(outFilePart1+"_intensity_brightest_txpnsite_xyfit_cy5_nofilter.txt")
            if ExclEmpty == 1:
                EmptyCy5 = list(where(CountCy5 == 0)[0])
            else:
                EmptyCy5 = []
            NoFitCy5 = list(where(TSCy5Pos[:,-1] == 0)[0])
            NoFit_EmptyCy5 = uniqueFromList(EmptyCy5 + NoFitCy5)
            NotEmptyCy5 = [x for x in range(0,len(CountCy5)) if x not in NoFit_EmptyCy5]
            corTScountCy5  = pearsonr(CountCy5[NotEmptyCy5], TSCy5[NotEmptyCy5])
            fig = plt.figure()
            plt.scatter(CountCy5[NotEmptyCy5], TSCy5[NotEmptyCy5], s = 2)
            plt.title("Cy5, correlation = "+str('{0:.3f}'.format(corTScountCy5[0])))
            plt.xlabel('Nr RNA/cell, Cy5')
            plt.ylabel('Nr of nascent RNA at TS')
            plt.savefig(outFilePart1+"_correlation_count_with_TSintensity_cy5.pdf", format = 'pdf')
            plt.close(fig)
            cors.append("\nCorrelation Cy5 TS intensity and cell count = "+str(corTScountCy5[0])+" , p-value = "+ str('{:.2E}'.format(corTScountCy5[1])))
            cors.append(str(len(NotEmptyCy5)) +" cells were analyzed, "+str(len(NoFit_EmptyCy5))+" cells were filtered")
            cors.append(str(len(EmptyCy5))+" cells removed because they were empty in Cy5 channel")
            cors.append(str(len(NoFitCy5))+" cells were removed because the Cy5 TS was not fit")

        ### Correlate between channels
        if len(channelsToAnalyze) < 2:
            print "Cannot calculate Cy3-Cy5 correlations in 1 color data"
        elif len(channelsToAnalyze) == 2:
            ### filter empty cells and spots that were not fit in either channel
            if ExclEmpty ==1:
                overlapEmpty = [x for x in EmptyCy3 if x in EmptyCy5]  #### only trows away cells where no spots in both Cy3 and Cy5.
            else:
                overlapEmpty = []
            NoFitCy3Cy5 = uniqueFromList(NoFitCy3 + NoFitCy5)
            NoFit_EmptyCy3Cy5 = uniqueFromList(overlapEmpty + NoFitCy3Cy5)
 #           overlapNotEmpty = [x for x in range(0,len(CountCy3)) if x not in NoFit_EmptyCy3Cy5]
            
            ### calculate and plot distances
            if localizeDim == "2D": dist = ((TSCy3Pos[:,1]-TSCy5Pos[:,1])**2 +(TSCy3Pos[:,2]-TSCy5Pos[:,2])**2)**0.5
            if localizeDim == "3D": dist = ((TSCy3Pos[:,1]-TSCy5Pos[:,1])**2 +(TSCy3Pos[:,2]-TSCy5Pos[:,2])**2+(TSCy3Pos[:,3]-TSCy5Pos[:,3])**2)**0.5
            savetxt(outFilePart1+"_distances_TSs.txt",dist, delimiter = "\n", fmt='%.8e')
            
            fig=plt.figure(figsize = (8,3))
            h,x=histogram(sqrt(dist),bins=r_[0:20:.2]);
            plt.xlabel('TS distances');
            plt.ylabel('Count');
            plt.title("Distances between Cy3 and Cy5 TS")
            plt.plot(x[:-1],h+1);
            plt.tight_layout();
            plt.savefig(outFilePart1+"_distances_TS.pdf", format = 'pdf')
            plt.close(fig)
            
            ### filter on distance
            if filterOnDistance == 1:
                filterDist = list(where(dist >= distThresh)[0])
            else:
                filterDist = []
            
            toFilter = uniqueFromList(NoFit_EmptyCy3Cy5 + filterDist)
            toAnalyze = [x for x in range(0,len(CountCy3)) if x not in toFilter]
       
            #### only select expressing population
            fracCy3OnlyOff = sum((TSCy3[toAnalyze] < onoffThreshCy3) & (TSCy5[toAnalyze] >= onoffThreshCy5))
            fracCy5OnlyOff = sum((TSCy3[toAnalyze] >= onoffThreshCy3) & (TSCy5[toAnalyze] < onoffThreshCy5))
            fracBothOff = sum((TSCy3[toAnalyze] < onoffThreshCy3) & (TSCy5[toAnalyze] < onoffThreshCy5))
            fracBothOn = sum((TSCy3[toAnalyze] >= onoffThreshCy3) & (TSCy5[toAnalyze] >= onoffThreshCy5))
            fractionList =[["Only Cy3 off", "Only Cy5 off", "Both off", "Both on"],[fracCy3OnlyOff, fracCy5OnlyOff, fracBothOff, fracBothOn],[float(fracCy3OnlyOff)/float(len(toAnalyze)), float(fracCy5OnlyOff)/float(len(toAnalyze)), float(fracBothOff)/float(len(toAnalyze)), float(fracBothOn)/float(len(toAnalyze))]]
            savetxt(outFilePart1+"_fraction_on_off.txt", fractionList, delimiter = "\t", fmt= "%s")
            
            
            overlap  = [value for value in where(TSCy3>=onoffThreshCy3)[0].tolist() if value in set(where(TSCy5>=onoffThreshCy5)[0].tolist())]
            overlap2 = [x for x in overlap if x not in toFilter]

            corTSCy3Cy5 = pearsonr(TSCy3[toAnalyze], TSCy5[toAnalyze])
            corTSCy3Cy5bg = pearsonr(TSCy3[overlap2], TSCy5[overlap2])
            corCountCy3Cy5 = pearsonr(CountCy3[toAnalyze], CountCy5[toAnalyze])
            
            fig = plt.figure()
            plt.scatter(CountCy3[toAnalyze], CountCy5[toAnalyze], s = 2)
            plt.title("Correlation = "+str('{0:.3f}'.format(corCountCy3Cy5[0]))+" , p-value = "+ str('{:.2E}'.format(corCountCy3Cy5[1])))
            plt.xlabel('Nr RNA/cell Cy3')
            plt.ylabel('Nr RNA/cell Cy5')
            plt.savefig(outFilePart1+"_correlation_Count_Nuc+Cytoplasm.pdf", format = 'pdf')
            plt.close(fig)
            
            fig = plt.figure()
            plt.scatter(TSCy3[toAnalyze], TSCy5[toAnalyze], s = 2)
            plt.title("Correlation = "+str('{0:.3f}'.format(corTSCy3Cy5[0])) +" , p-value = "+ str('{:.2E}'.format(corTSCy3Cy5[1])))
            plt.xlabel('Nr of nascent RNA at TS, Cy3')
            plt.ylabel('Nr of nascent RNA at TS, Cy5')
            plt.savefig(outFilePart1+"_correlation_TSintensity.pdf", format = 'pdf')
            plt.close(fig)
            
            fig = plt.figure()
            plt.scatter(TSCy3[overlap2], TSCy5[overlap2], s = 2)
            plt.title("Correlation = "+str('{0:.3f}'.format(corTSCy3Cy5bg[0])) +" , p-value = "+ str('{:.2E}'.format(corTSCy3Cy5bg[1])))
            plt.xlabel('Nr of nascent RNA at TS, Cy3')
            plt.ylabel('Nr of nascent RNA at TS, Cy5')
            plt.savefig(outFilePart1+"_correlation_TSintensity_abovethreshold.pdf", format = 'pdf')
            plt.close(fig)
            
            cors.append("\nCorrelation Cy3-Cy5 TS intensity = "+str(corTSCy3Cy5[0]) +" , p-value = "+ str('{:.2E}'.format(corTSCy3Cy5[1])))
            cors.append("Correlation Cy3-Cy5 TS intensity (excluded below thresholds "+str(onoffThreshCy3)+" for cy3, "+str(onoffThreshCy5)+" for cy5) = "+str(corTSCy3Cy5bg[0]) +" , p-value = "+ str('{:.2E}'.format(corTSCy3Cy5bg[1])))
            cors.append("Correlation Cy3-Cy5 cell count = "+str(corCountCy3Cy5[0])+" , p-value = "+ str('{:.2E}'.format(corCountCy3Cy5[1])))
            cors.append(str(len(toAnalyze)) +"cell were analyzed, "+str(len(toFilter))+" cells were filtered")
            cors.append(str(len(overlapEmpty))+" cells removed because they were empty in both Cy3 and Cy5 channel")
            cors.append(str(len(NoFitCy3Cy5))+" cells were removed because they were not fit in either Cy3 or Cy5 of the channels")
 
            TSCy3 = TSCy3 [toAnalyze]
            TSCy5 = TSCy5 [toAnalyze]
            threshCy3 = range(0,int(max(TSCy3)))
            threshCy5 = range(0,int(max(TSCy5)))
            corMatrixL = zeros((len(threshCy3),len(threshCy5)))
            corMatrixLS = zeros((len(threshCy3),len(threshCy5)))
            corMatrixS = zeros((len(threshCy3),len(threshCy5)))
            corMatrixSL = zeros((len(threshCy3),len(threshCy5)))
            for m in threshCy3:
                for n in threshCy5:
                    overlapL  = [value for value in where(TSCy3>=m)[0].tolist() if value in set(where(TSCy5>=n)[0].tolist())]
                    overlapS  = [value for value in where(TSCy3<m)[0].tolist() if value in set(where(TSCy5<n)[0].tolist())]
#                    non_overlapL = range(0,len(TSCy3))
#                    for i in overlapL: non_overlapL.remove(i)
                    non_overlapS = range(0,len(TSCy3))
                    for i in overlapS: non_overlapS.remove(i)
                    
                    if len(overlapL)>10: corMatrixL[m,n] = pearsonr(TSCy3[overlapL], TSCy5[overlapL])[0]
#                    if len(overlapS)>10: corMatrixS[m,n] = pearsonr(TSCy3[overlapS], TSCy5[overlapS])[0]
#                    if len(non_overlapL)>10: corMatrixLS[m,n] = pearsonr(TSCy3[non_overlapL], TSCy5[non_overlapL])[0]
                    if len(non_overlapS)>10: corMatrixSL[m,n] = pearsonr(TSCy3[non_overlapS], TSCy5[non_overlapS])[0]
                    
      #      fig = plt.figure()
            fig,(ax1,ax2) = plt.subplots(1,2)
            im1 = ax1.pcolor(corMatrixL, cmap = "RdBu_r", vmin = -1, vmax = 1)
            ax1.set_title("Corr TS intensity Cy3-Cy5, \nAND>=")
            ax1.set_xlabel("larger than threshold, Cy5")
            ax1.set_ylabel("larger than threshold, Cy3")
              
            im2 = ax2.pcolor(corMatrixSL, cmap = "RdBu_r", vmin = -1, vmax = 1 )
            ax2.set_title("Corr TS intensity Cy3-Cy5, \nOR>=")
            ax2.set_xlabel("larger than threshold, Cy5")
#            ax2.set_ylabel("larger than threshold, Cy3")
            
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im2, cax=cax)
            
#            fig.tight_layout()
            plt.savefig(outFilePart1+"_correlation_TS_intensity_heatmaps_threshold.pdf", format = 'pdf')



        savetxt(outFilePart1+"_correlations.txt",cors, delimiter = "\n", fmt="%s")
        
    return (CelldistrCy3, NucldistrCy3, TSdistrCy3,CelldistrCy5, NucldistrCy5, TSdistrCy5)

############################################################################################################
#### running function on paths
############################################################################################################
CelldistrCy3 = []
NucldistrCy3 = []
TSdistrCy3 = []
CelldistrCy5 = []
NucldistrCy5 = []
TSdistrCy5 = []

# This part tests whether the distributions have already been calculated and are stored in the memory. If so, they are not recaculcated. If they are not stored in the memory, they will be recalculated.
if UsePrecalcDistr == 1:
    try: rawdistributions
    except NameError:
        print "Analyzing datsets"
        for file in range(0,lfileListIn):
            print 'Analyzing smFISH dataset '+str(file+1)+' out of '+str(lfileListIn)+': '+fileListIn[file]
            rawdistributions = analysis(fileListIn[file])
    else: print "Using pre-analyzed datasets"
else:
    print "Analyzing datsets"
    for file in range(0,lfileListIn):
        print 'Analyzing smFISH dataset '+str(file+1)+' out of '+str(lfileListIn)+': '+fileListIn[file]
        rawdistributions = analysis(fileListIn[file])

distributions = rawdistributions

if CombineReps == 1:
    print "Combining datasets"
    distributionsnew = [[],[],[],[],[],[]]
    
    for channel in channelsToAnalyze:
        if channel == 0: color = "cy3"
        if channel == 1: color = "cy5"

        for k in [0,1,2]:
            distributionselement = []
            for i in range(0,len(repsToCombine)):
                combi = []
                for j in repsToCombine[i]:
                    combi.append(distributions[3*channel+k][j-1])
                distributionselement.append(concatenate(combi))
            distributionsnew[3*channel+k] = distributionselement
    distributions = distributionsnew

if CombineReps == 1:
    print "Making histograms of combined datasets"
    
    os.system("cp FISH_pipeline_parameters.py "+analysisfolder+"FISH_pipeline_parameters.py")

    for channel in channelsToAnalyze:
        if CalcFracOff == 1 and exclBins == 0 and exclBinsCount == 0 and ExclEmpty == 1: FracOff = []
        for i in range(0,len(repsToCombine)):
            
             ### for TS distribution
            if channel == 0: color = "cy3"; bins=r_[:nrbins_cy3_ts]; xCenterTS = nrbins_cy3_ts/2; freqAxisTS = freqAxisTS_cy3
            if channel == 1: color = "cy5"; bins=r_[:nrbins_cy5_ts]; xCenterTS = nrbins_cy5_ts/2; freqAxisTS = freqAxisTS_cy5

            hNorm = histTSvals(distributions[3*channel+2][i],bins,exclBins)
            bins=bins[:-1]

            TSdistr = distributions[3*channel+2][i]
            if CalcFracOff == 1 and exclBins == 0 and exclBinsCount == 0 and ExclEmpty == 1: FracOff.append(hNorm[0])

            #### Fitting models
            if FitBurstCombi == 1: [Pres, n, aT] = fitBurstModel(hNorm, bins, exclBins)
            else: [Pres, n, aT] = [0, 0, 0]

            if FitPoissonCombi == 1: [PPoissres, cT] = fitPoissonModel(hNorm, bins, exclBins)
            else: [PPoissres, cT] = [0, 0]

            #### making TS plot

            histTS(TSdistr, hNorm[exclBins:], bins[exclBins:], FitBurstCombi, FitPoissonCombi, Pres, n, aT, PPoissres, cT, freqAxisTS, xCenterTS, exclBins, color, analysisfolder+"txpnsite_histogram_"+color+"_merged_distr_"+str(i)+".pdf")

            ## nuclear and cellular count distribution
        
            if channel == 0: bins=r_[:nrbins_cy3]; freqAxisNuc = freqAxisNuc_cy3; freqAxisCell = freqAxisCell_cy3; xCenter = nrbins_cy3/2
            if channel == 1: bins=r_[:nrbins_cy5]; freqAxisNuc = freqAxisNuc_cy5; freqAxisCell = freqAxisCell_cy5; xCenter = nrbins_cy5/2

            histNuc(distributions[3*channel+1][i], bins, freqAxisNuc, xCenter, exclBinsCount, color, analysisfolder+"nuclear_count_"+color+"_merged_distr_"+str(i)+".pdf")

            histCell(distributions[3*channel][i], bins, freqAxisCell, xCenter, exclBinsCount, color, analysisfolder+"cell_count_"+color+"_merged_distr_"+str(i)+".pdf")

        ## write statistics file
            writestatsfile(analysisfolder+"distr_"+str(i), color, distributions[3*channel][i], distributions[3*channel+1][i], exclBinsCount, TSdistr, exclBins)
        
        if CalcFracOff == 1 and exclBins == 0 and exclBinsCount == 0 and ExclEmpty == 1: savetxt(analysisfolder+"fraction_off_"+color+".txt", FracOff, delimiter='\t')

if Violins == 1:
    print "Making violin plots"

    for channel in channelsToAnalyze:
        if channel == 0: color = "cy3"; yAxisViolinCell = yAxisViolinCell_cy3; yAxisViolinNuc = yAxisViolinNuc_cy3; yAxisViolinTS = yAxisViolinTS_cy3
        if channel == 1: color = "cy5"; yAxisViolinCell = yAxisViolinCell_cy5; yAxisViolinNuc = yAxisViolinNuc_cy5; yAxisViolinTS = yAxisViolinTS_cy5

        countdatasets = len(distributions[3*channel])
        if countdatasets != 1 and len(pvalsViolins[0]) != 0:
            pcell = [];
            pcellround = [];
            pnuc = [];
            pnucround = [];
            pTS = [];
            pTSround = [];
            for ppair in range(0,len(pvalsViolins)):
                dist1 = pvalsViolins[ppair][0]-1
                dist2 = pvalsViolins[ppair][1]-1
                z, p = scipy.stats.mannwhitneyu(distributions[3*channel][dist1],distributions[3*channel][dist2])
                p = 2*p
                pcell.append(p)
                pcellround.append(round(p,4))
                z, p = scipy.stats.mannwhitneyu(distributions[3*channel+1][dist1],distributions[3*channel+1][dist2])
                p = 2*p
                pnuc.append(p)
                pnucround.append(round(p,4))
                z, p = scipy.stats.mannwhitneyu(distributions[3*channel+2][dist1],distributions[3*channel+2][dist2])
                p = 2*p
                pTS.append(p)
                pTSround.append(round(p,4))

        fig = plt.figure()
        plt.violinplot(distributions[3*channel], showmeans = True, showextrema = False)
        plt.ylabel('Number of RNA/Cell')
        plt.xlabel('Dataset number')
        y_max = max(concatenate(distributions[3*channel]))
        x_max = len(distributions[3*channel])
        if yAxisViolinCell != "None": y_max = yAxisViolinCell
        plt.ylim([0,y_max])
        if countdatasets != 1 and len(pvalsViolins[0]) != 0:
            plt.text(x_max+0.3, 0.98*y_max, "2-sided Mann-Whitney U-tests:\npvalues calculated for datasets: "+str(pvalsViolins)+"\npvalues: "+str(pcellround), fontsize = 12, horizontalalignment='right', verticalalignment='top')
        plt.savefig(analysisfolder+"Violin_Cellular_"+color+".pdf")
        plt.close(fig)

        fig = plt.figure()
        plt.violinplot(distributions[3*channel+1], showmeans = True, showextrema = False)
        plt.ylabel('Number of RNA/Nucleus')
        plt.xlabel('Dataset number')
        y_max = max(concatenate(distributions[3*channel+1]))
        x_max = len(distributions[3*channel+1])
        if yAxisViolinNuc != "None": y_max = yAxisViolinNuc
        plt.ylim([0,y_max])
        if countdatasets != 1 and len(pvalsViolins[0]) != 0:
            plt.text(x_max+0.3, 0.98*y_max, "2-sided Mann-Whitney U-tests:\npvalues calculated for datasets: "+str(pvalsViolins)+"\npvalues: "+str(pnucround), fontsize = 12, horizontalalignment='right', verticalalignment='top')
        plt.savefig(analysisfolder+"Violin_Nuclear_"+color+".pdf")
        plt.close(fig)

        fig = plt.figure()
        plt.violinplot(distributions[3*channel+2], showmeans = True, showextrema = False)
        plt.ylabel('Number of nascent RNA at TS')
        plt.xlabel('Dataset number')
        y_max = max(concatenate(distributions[3*channel+2]))
        x_max = len(distributions[3*channel+2])
        if yAxisViolinTS != "None": y_max = yAxisViolinTS
        plt.ylim([0,y_max])
        if countdatasets != 1 and len(pvalsViolins[0]) != 0:
            plt.text(x_max+0.3, 0.98*y_max, "2-sided Mann-Whitney U-tests:\npvalues calculated for datasets: "+str(pvalsViolins)+"\npvalues: "+str(pTSround), fontsize = 12, horizontalalignment='right', verticalalignment='top')
        plt.savefig(analysisfolder+"Violin_TS_"+color+".pdf")
        plt.close(fig)

        #write pvalues to separate file
        if len(pvalsViolins[0]) != 0:
            pvalsfile = open(analysisfolder+"pvalues_"+color+".txt", "w+")
            pvalsfile.write("p-values cellular violin plots: "+str(pcell))
            pvalsfile.write("\np-values nuclear violin plots: "+str(pnuc))
            pvalsfile.write("\np-values TS violin plots: "+str(pTS))
            pvalsfile.close()


if CompareNormValues == 1:
    print "Comparing normalization values"
    for channel in channelsToAnalyze:
        if channel == 0: color = "cy3"
        if channel == 1: color = "cy5"
        normVals = []
        datasets = []
        for file in range(0,len(fileListIn)):
            pathIn = fileListIn[file]
            date = pathIn.split("/")[-3]
            expName = pathIn.split("/")[-2]
            fileName = outputfolder+date+"_"+expName+"/"+date+"_"+expName+"_cyto_normalization_value_"+color+".txt"
            normVal = loadtxt(fileName)
            normVals.append(normVal)
            datasets.append(date+"_"+expName)
        
        normtable = c_[range(1,len(fileListIn)+1),datasets, normVals]
        savetxt(analysisfolder+"Cyto_normalization_values_"+color+".txt", normtable, fmt="%s", delimiter='\t')
    
        fig = plt.figure()
        plt.bar(range(1,len(fileListIn)+1),normVals)
        plt.ylim([0,max(normVals)*1.2])
        plt.title("Normalizing value from cytoplasmic spots, "+color)
        plt.xlabel('Dataset number')
        plt.ylabel('Intensity')
        plt.savefig(analysisfolder+"Cyto_normalization_values_plot_"+color+".pdf", format = 'pdf')
        plt.close(fig)









