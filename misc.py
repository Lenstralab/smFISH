"Some of this code was adapted from our collaborator Antoine Coulon (2020) <antoine.coulon@curie.fr>"

from scipy import *; from scipy import special, linalg, fftpack
from PIL import Image as pil
import copy as copy2, sys, os, time, cPickle, commands, re

cmdExists=lambda cmd: (len(commands.getoutput('command -v '+cmd))>0)*cmd
if sys.platform=='linux2':
  viewerImg=cmdExists('eog') or ('false',sys.stdout.write('Warining: No image viewer found.\n'))[0]
  viewerPdf=cmdExists('evince') or ('false',sys.stdout.write('Warining: No PDF viewer found.\n'))[0]
  fijiCmd=cmdExists('fiji') or cmdExists('imagej')\
    or ('false',sys.stdout.write('Warining: Fiji or ImageJ not found.\n'))[0]
elif sys.platform=='darwin':
  viewerImg=cmdExists('open -a /Applications/Preview.app/Contents/MacOS/Preview')\
    or ('false',sys.stdout.write('Warining: No image viewer found.\n'))[0]
  viewerPdf=cmdExists('open -a /Applications/Preview.app/Contents/MacOS/Preview')\
    or ('false',sys.stdout.write('Warining: No PDF viewer found.\n'))[0]
  fijiCmd=cmdExists('/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx')\
    or cmdExists('/Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub')\
    or ('false',sys.stdout.write('Warining: Fiji or ImageJ not found.\n'))[0]


class objFromDict:
  "Usage: objFromDict(**dictionary).\nPrint give the list of attributes."
  def __init__(self, **entries):
    self.__dict__.update(entries)
    for a in self.__dict__.keys():
      self.__dict__[a.replace('-','_').replace('*','_').replace('+','_').replace('/','_')]=self.__dict__.pop(a)
  def __repr__(self): return '** objFromDict attr. --> '+', '.join(filter((lambda s: (s[:2]+s[-2:])!='____'),dir(self)))



#----------------------------
#import Gnuplot;
#class GpDisplays:
#  def __init__(self,n=1):
#    self.gs=[]; self.addDisp(n)
#  def __getitem__(self,i):
#    if len(self.gs)>i: return self.gs[i]
#    else:
#      print "Not enough displays. Use g.addDisp(n)."
#      return (lambda a: None)
#  def addDisp(self,n=1):
#    for i in range(n): self.gs.append(Gnuplot.Gnuplot())
#  def __call__(self,a): self.gs[0](a)
#  def __repr__(self): return "%d display(s)."%len(self.gs)
#if not (globals().has_key('GpRenew') and GpRenew==False):
#  g=GpDisplays(6); g1=Gnuplot.Gnuplot(); g2=Gnuplot.Gnuplot(); g3=Gnuplot.Gnuplot(); g4=Gnuplot.Gnuplot(); g5=Gnuplot.Gnuplot(); g6=Gnuplot.Gnuplot();
#
#
#class GpDisplaysPDF:
#  def __init__(self,n=1):
#    self.gs=[]; self.formatSpec=[]; self.isNew=[]; self.addDisp(n); self.fn='g%02d.pdf'
#  def __getitem__(self,i):
#    if self.fn.find('%')==-1: fn=self.fn
#    else: fn=self.fn%i
#    def gPDFloc(s='',disp=1,extra=None):
#      if (s.find('reset')==-1) and self.isNew[i] and len(s): s='reset; '+s
#      self.gs[i](s.replace('reset','reset; '+self.formatSpec[i]+'set output "%s";'%fn));
#      if disp:
#        self.gs[i]=Gnuplot.Gnuplot(); self.isNew[i]=1
#        if disp==1: os.system(viewerPdf+' '+re.escape(fn)+'&')
#        elif disp==2:
#          if extra==None: extra=(.05,0)
#          os.system('convert -delay %f -loop %d -density 130 %s -bordercolor white -border 0 -alpha off -layers Optimize %s; firefox %s'%((extra[0]*100,extra[1],fn,)+(fn[:-4]+'.gif',)*2))
#      else: self.isNew[i]=0
#    if len(self.gs)>i: return gPDFloc
#    else:
#      print "Not enough displays. Use g.addDisp(n)."
#      return (lambda a: None)
#  def addDisp(self,n=1):
#    for i in range(n):
#      self.gs.append(Gnuplot.Gnuplot())
#      self.formatSpec.append('set term pdf size 6,4.5 lw 2 font "Helvetica,9"; ')
#      self.isNew.append(1)
#  def __call__(self,*a): self[0](*a)
#  def __repr__(self): return "%d display(s).\n"%len(self.gs)+"  .fn='%s'\n"%self.fn+'\n'.join(["  .formatSpec[%d]='%s'"%(i,self.formatSpec[i]) for i in range(len(self.formatSpec))])
#if not (globals().has_key('GpRenew') and GpRenew==False): gPDF=GpDisplaysPDF(6);
### Example:
##gPDF[0]('reset; unset key; set zeroa;',0); [gPDF[0]('plot [-1:1] [-1:1] sin((x+%e)*2*pi)'%a,0) for a in r_[0:1:.05]]*0; gPDF[0]()
##gPDF[0]('',2)
#
#
#
#gTmp=(lambda s: Gnuplot.Gnuplot(persist=1)(s))
#
#if not globals().has_key('GpRenew'): GpRenew=True;
#
#def GpD(e,**args):
#  if e.shape[0]: d=Gnuplot.Data(e)
#  else: return '"<(echo \'\')"'
#  if args.has_key('lines'):
#    cont=d.content.split('\n'); lines=args['lines']
#    d.content='\n'.join(['"%s" %s'%(lines[i].replace('\n','\\n'),cont[i]) for i in range(len(cont)-2)]+['',''])
#  if args.has_key('cols'): d.content=args.has_key('cols')*'"" '+' '.join(['"%s"'%a for a in args['cols']])+'\n'+d.content
#  return d.get_base_command_string()
#def GpD3d(data,n):
#  name=GpD(data)[1:-1]; name2=name.replace('gnuplot/','gnuplot_')+"_3d"
#  os.system("cat "+name+" | sed '"+"n;"*(n-1)+"G;' > "+name.replace('gnuplot/','gnuplot_')+"_3d")
#  return "'%s'"%name2
#def GpD3d2(data,n):
#  name=GpD(data)[1:-1]; name2=name.replace('gnuplot/','gnuplot_')+"_3d"
#  os.system("cat "+name+" | sed '"+"n;"*(n-1)+"G;G;' > "+name.replace('gnuplot/','gnuplot_')+"_3d")
#  return "'%s'"%name2


#----------------------------
from scipy import *;

def image2array(f_im):
  im=pil.open(f_im)
  if im.mode == "RGB": return fromstring(im.tostring(), uint8).reshape(im.size[1],im.size[0],3)
  elif im.mode == "RGBA": return fromstring(im.tostring(), uint8).reshape(im.size[1],im.size[0],4)
  elif im.mode == "L": return fromstring(im.tostring(), uint8).reshape(im.size[1],im.size[0])
  raise ValueError, "Not able to convert. Mode is "+im.mode

def array2image(a,f_im):
  if len(a.shape)==3:
    if a.shape[2]==3: return pil.fromstring("RGB", (a.shape[1], a.shape[0]), a.astype(uint8).tostring()).save(f_im)
    elif a.shape[2]==4: return pil.fromstring("RGBA", (a.shape[1], a.shape[0]), a.astype(uint8).tostring()).save(f_im)
  elif len(a.shape)==2: return pil.fromstring("L", (a.shape[1], a.shape[0]), a.astype(uint8).tostring()).save(f_im)
  raise ValueError, "Incorrect data format. Shape is "+str(a.shape)

def tiff2array(fn,shape=(-1,)):
  im=pil.open(fn); res=[]
  try:
    while 1:
      res.append(array(im.getdata())); im.seek(im.tell()+1)
  except EOFError: pass
  return array(res).reshape(shape+im.size[::-1])

def array2tiff_Old(a,fn):
  if not (type(a.flatten()[0]) in [uint8,uint16,uint32]): print "Error: Data type has to be either uint8, uint16 or uint32. e.g. array2tiff(uint16(data),'file.tif')"
  else:
    array2D2tiff=(lambda a,fn: pil.fromstring({uint8:'L', uint16:'I;16', uint32:'I;32'}[type(a.flatten()[0])],a.shape[::-1],a.tostring()).save(fn))
    if len(a.shape)==2: array2D2tiff(a,fn)
    elif len(a.shape)>2:
      a=a.reshape((-1,)+a.shape[-2:])
      for i in r_[:a.shape[0]]: array2D2tiff(a[i],'.tmp-%09d-%s'%(i,fn))
      os.system('convert -compress None .tmp-?????????-%s %s; rm .tmp-?????????-%s'%(fn,fn,fn))
    else: print "!! Wrong dimensions !!"

def array2tiff(a,fn):
  if not (type(a.flatten()[0]) in [uint8,uint16,uint32]): print "Error: Data type has to be either uint8, uint16 or uint32. e.g. array2tiff(uint16(data),'file.tif')"
  else:
    array2D2tiff=(lambda a,fn: pil.fromstring({uint8:'L', uint16:'I;16', uint32:'I;32'}[type(a.flatten()[0])],a.shape[::-1],a.tostring()).save(fn))
    if len(a.shape)==2: array2D2tiff(a,fn)
    elif len(a.shape)>2:
      a=a.reshape((-1,)+a.shape[-2:])
      for i in r_[:a.shape[0]]: array2D2tiff(a[i],'%s.tmp-%09d.tif'%(fn,i))
      os.system('convert -compress None %s.tmp-?????????.tif %s; rm %s.tmp-?????????.tif'%(fn,fn,fn))
    else: print "!! Wrong dimensions !!"

def array2RGBtiff(a,fn):
  if type(a.flatten()[0])!=uint8: print "Error: Data type has to be uint8. e.g. array2RGBtiff(uint8(data),'file.tif')"
  else:
    array2D2tiff=(lambda a,fn: pil.fromstring('RGB',a.shape[:0:-1],a.T.swapaxes(0,1).astype(uint8).tostring()).save(fn))
    if len(a.shape)<3: print "!! Wrong dimensions !!"
    else:
      if a.shape[-3]==2: a=a.swapaxes(-1,-3); a=c_[a,zeros(a.shape[:-1]+(1,))].swapaxes(-1,-3)
      if len(a.shape)==3: array2D2tiff(a,fn)
      else:
        a=a.reshape((-1,)+a.shape[-3:])
        for i in r_[:a.shape[0]]: array2D2tiff(a[i],'.tmp-%09d-%s'%(i,fn))
        os.system('convert -compress None .tmp-?????????-%s %s; rm .tmp-?????????-%s'%(fn,fn,fn))


def saveTiffMultipage(lim,fn, rescaleTo8bit=False, rgbOrder="rgba", **params):
    """Adapted from priithon (http://code.google.com/p/priithon/)
    lim=list of images."""
    fp = open(fn, 'w+b')
    ifd_offsets=[]
    params["_debug_multipage"] = True
    for z in range(len(lim)):
        ii = lim[z]
        fp.seek(0,2)
        if z==0: ifdOffset = 8
        else: ifdOffset = fp.tell()
        ii.save(fp, format="TIFF", **params)
        if z>0: # correct "next" entry of previous ifd -- connect !
            ifdo = ifd_offsets[-1]
            fp.seek(ifdo)
            ifdLength = ii._debug_multipage.i16(fp.read(2))
            fp.seek(ifdLength*12,1) # go to "next" field near end of ifd
            fp.write(ii._debug_multipage.o32( ifdOffset ))
        ifd_offsets.append(ifdOffset)
    fp.close()
  #e.g. RGB-color: saveTiffMultipage([    pil.merge('RGB',[pil.fromarray(uint8(rand(100,100)*255)) for a in r_[:3]])    for i in range(10)], 'tmpPy.tif'); os.system(fijiCmd+' --no-splash tmpPy.tif')
  #e.g 8-bit grayscale: saveTiffMultipage([    pil.fromarray(uint8(rand(100,100)*255))    for i in range(10)], 'tmpPy.tif'); os.system(fijiCmd+' --no-splash tmpPy.tif')
  #e.g 16-bit grayscale: saveTiffMultipage([    pil.fromarray(int16(meshgrid(r_[:1:.01],ones(100))[0]*65535))    for i in range(10)], 'tmpPy.tif'); os.system(fijiCmd+' --no-splash tmpPy.tif')
  #e.g 16-bit with multipl channels and slices: saveTiffMultipage([    pil.fromarray(int16(rand(100,100)*65535))    for i in range(3*4*10)], 'tmpPy.tif'); os.system(fijiCmd+""" --no-splash -eval 'open("tmpPy.tif"); run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=4 frames=10 display=Composite");'&""") #run("Save", "tmpPy.tif")


rasterList=['0110000100011001111010000111100110011110011000110000000000000000011110',
            '1001001100100100001010100100001000000010100101001000000000000110010010',
            '1001000100001000010011110111001110000100011000111011110000001111010010',
            '1001000100010001001000100000101001001000100100001000000011001000010010',
            '0110000100111100110000100111000110001000011000110000000011000110011110',]
rasterList=array([[map(int,list(a[i*5:(i+1)*5])) for a in rasterList] for i in r_[:13]])
def rasterTxt(s):
  if not type(s)==str: s=str(s)
  res=[]
  for c in s:
    try:
      res.append(rasterList['0123456789-.e'.index(c)])
    except: res.append(rasterList[-1])
  return array(res).swapaxes(0,1).reshape(5,-1)


#----------------------------
lut_gray=lambda vv: (lambda x: c_[x,x,x].clip(0,1))(vv.reshape(-1)).reshape(vv.shape+(3,))
defaultLut=lut_gray
norm01=(lambda a: (lambda b: b/max(b.flatten()))(a*1.-min(a.flatten())))

def showMat(m,lut=defaultLut,fname='tmpImg.png',newWindow=True,detach=True,disp=True): array2image((lut(m)*255),fname); os.system((viewerImg+' '+'-n '*newWindow+fname+'&'*detach)*(disp!=0))
def showImg(img,fname='tmpImg.png',newWindow=True,detach=True,disp=True): array2image(img*255,fname); os.system((viewerImg+' '+'-n '*newWindow+fname+'&'*detach)*(disp!=0))


#----------------------------
histoF=(lambda hr: c_[hr[1][1:],hr[0]])
histoFlog10=(lambda hr: c_[10**((hr[1][1:]+hr[1][:-1])/2.),hr[0]])
cumDistF=(lambda hr: c_[hr[1][1:],cumsum(hr[0])/sum(hr[0])])
cumDistFlog10=(lambda hr: c_[10**((hr[1][1:]+hr[1][:-1])/2.),cumsum(hr[0])/sum(hr[0])])


#----------------------------
def leastsqWrap(func, x0, force=None, args=(), Dfun=None, col_deriv=0, ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag2=None):
  if force==None: force=zeros(x0.shape[0])
  wif=where(force==0)[0]
  func2=lambda x: func(array([x[where(wif==i)[0][0]] if i in wif else x0[i] for i in r_[:x0.shape[0]]]))
  resFit=optimize.leastsq(func2, x0[wif], args=args, Dfun=Dfun, full_output=1, col_deriv=col_deriv, ftol=ftol, xtol=xtol, gtol=gtol, maxfev=maxfev, epsfcn=epsfcn, factor=factor, diag=diag2)
  x,err=resFit[0],diag(resFit[1])**.5
  x=array([x[where(wif==i)[0][0]] if i in wif else x0[i] for i in r_[:x0.shape[0]]])
  err=array([err[where(wif==i)[0][0]] if i in wif else 0. for i in r_[:x0.shape[0]]])
  return x,err,resFit

#----------------------------
def uniqueFromList(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list
#----------------------------
#findCells

import skimage.feature
import scipy, skimage
import numpy as np
from tqdm.auto import tqdm

def findcells(im, imnuc=None, cellcolormask=None, ccdist=None, threshold=None, thresholdnuc=None, removeborders=True):
    #""" segement cells and nuclei from an image (nxm array)
    #    wp@tl20190710
    #
    #    im: 2d numpy array with an image
    #
    #    optional:
    #    imnuc: 2d numpy array with an image containing the nuclei, for example with DAPI
    #    cellcolormask: use a cellmask from another frame to assign the same colors
    #                   to the cells found by findcells
    #    ccdist:        (approximate) minimum distance between the centers of cells, good values: yeast: 25, mamalian: 150
    #    threshold:     value used to threshold the image, default: use Otsu
    #    thresholdnuc:  value used to threshold the image with nuclei, default: use Otsu
    #    removeborders: remove any cells (and their nuclei) touching any borders (case with nan's unhandled atm), default: True
    #"""
    
    # ----------------------- Helper functions ---------------------------------
    def collectr(cf):
        #""" Makes a 1d corrfun out of a 2d corrfun
        #    The result of the maximum of cf at each distance from the center
        #    wp@tl20191102
        #"""
        x, y = np.meshgrid(range(cf.shape[0]), range(cf.shape[1]))
        c = [(i-1)/2 for i in cf.shape]
        r = np.sqrt((x-c[0])**2 + (y-c[1])**2).flatten()
        idx = np.argsort(r)
        r = np.round(r[idx]).astype('int')
        cf = np.round(cf.flatten()[idx]).astype('int')
        rr = np.unique(r)
        cf = [np.max(cf[r==i]) for i in rr]
        return rr, cf
    
    def maskpk(pk, mask):
        #    """ remove points in nx2 array which are located outside mask
        #        wp@tl20190709
        #        """
        pk = np.round(pk)
        idx = []
        for i in range(pk.shape[0]):
            if mask[pk[i,0], pk[i,1]]:
                idx.append(i)
        return pk[idx,:]
    
    def disk(s):
        #    """ make a disk shaped structural element to be used with
        #        morphological functions
        #        wp@tl20190709
        #        """
        d = zeros((s,s))
        c = (s-1)/2
        x, y = meshgrid(range(s), range(s))
        d2 = (x-c)**2 + (y-c)**2
        d[d2<s**2/4] = 1
        return d
    
    def fill_nan(im):
        #""" Assigns the value of the nearest finite valued pixel to infinite valued pixels
        #    wp@tl20190910
        #"""
        im = im.copy()
        a = np.where(np.isfinite(im))
        b = np.where(~np.isfinite(im))
        v = scipy.interpolate.griddata(a, im[a[0], a[1]], b, 'nearest')
        for i, j, w in zip(b[0], b[1], v):
            im[i, j] = w
        return im
    
    def collectr(cf):
        #""" Makes a 1d corrfun out of a 2d corrfun
        #    The result of the maximum of cf at each distance from the center
        #    wp@tl20191102
        #"""
        x, y = np.meshgrid(range(cf.shape[0]), range(cf.shape[1]))
        c = [(i-1)/2 for i in cf.shape]
        r = np.sqrt((x-c[0])**2 + (y-c[1])**2).flatten()
        idx = np.argsort(r)
        r = np.round(r[idx]).astype('int')
        cf = np.round(cf.flatten()[idx]).astype('int')
        rr = np.unique(r)
        cf = [np.max(cf[r==i]) for i in rr]
        return rr, cf
    
    def corrfft(im,jm):
        #'''
        #% usage: d, cfunc = corrfft(images)
        #%
        #% input:
        #%   im, jm: images to be correlated
        #%
        #% output:
        #%   d:      offset (x,y) in px
        #%   cfunc:  correlation function
        #'''
        
        im = im.astype(float)
        jm = jm.astype(float)
        
        im -= np.nanmean(im)
        im /= np.nanstd(im)
        jm -= np.nanmean(jm)
        jm /= np.nanstd(jm)
        
        im[np.isnan(im)] = 0
        jm[np.isnan(jm)] = 0
        
        nY = np.shape(im)[0]
        nX = np.shape(im)[1]
        
        cfunc = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fft2(im)*np.conj(np.fft.fft2(jm)))))
        y, x = np.unravel_index(np.nanargmax(cfunc), cfunc.shape)
        
        d = [x-np.floor((nX)/2), y-np.floor(nY/2)]
        
        #peak at x=nX-1 means xoffset=-1
        if d[0]>nX/2:
            d[0] -= nX
        if d[1]>nY/2:
            d[1] -= nY
        return d, cfunc
    
    def otsu_local(im, res=16):
        th = np.zeros((res, res))
        s = im.shape[0]/res
        for i in range(res):
            for j in range(res):
                th[i, j] = skimage.filters.threshold_otsu(im[i*s:(i+1)*s, j*s:(j+1)*s])
        th = skimage.transform.resize(th, im.shape, anti_aliasing=False, mode='edge')
        return scipy.ndimage.gaussian_filter(th, 256, mode='reflect')
    
    # ----------------------- Finally the function body ---------------------------------
    
    im = im.astype('float')
    if not imnuc is None:
        imnuc = imnuc.astype('float')
    
    if ccdist is None:
        if imnuc is None:
            cf = corrfft(im.copy(), im.copy())[1]
        else:
            cf = corrfft(imnuc.copy(), imnuc.copy())[1]
        r, cf = collectr(cf)
        ccdist = r[scipy.signal.find_peaks(-np.array(cf))[0][0]]
        ccdist = np.round(np.clip(ccdist, 10, np.sqrt(np.prod(im.shape))/5)).astype('int')

    if imnuc is None:
        p = np.where(np.isfinite(im))
        q = np.where(~np.isfinite(im))
        LB = fill_nan(im)
        LA = scipy.ndimage.gaussian_filter(LB, 2.5, mode='nearest')
        th = threshold or skimage.filters.threshold_otsu(LA[p[0], p[1]])
        mask = np.zeros(im.shape)
        mask[LA>th] = 1
        mask[q[0], q[1]] = 0
        mask = skimage.morphology.binary_opening(mask, disk(5)) #remove small features
#        mask = skimage.morphology.binary_dilation(mask, disk(5))
        mask = skimage.morphology.binary_erosion(mask, disk(5))
        pk = skimage.feature.peak_local_max(fill_nan(LA), footprint=disk(ccdist), exclude_border=False)
        pk = maskpk(pk, mask)
        markers = np.zeros(im.shape)
        for i in range(pk.shape[0]):
            if cellcolormask is None:
                markers[pk[i,0], pk[i,1]] = i+1
            else:
                markers[pk[i,0], pk[i,1]] = cellcolormask[pk[i,0], pk[i,1]]
        WS = skimage.segmentation.watershed(-LA, markers=markers, mask=mask)
        lbl = list(np.unique(WS))
        lbl.remove(0)
        cells = np.zeros(im.shape)
        nuclei = np.zeros(im.shape)
        for i in tqdm(lbl, disable=len(lbl)<25):
            cellmask = scipy.ndimage.morphology.binary_fill_holes((WS==i))
            cell = (cellmask*LA).flatten()
            cell = cell[cell>0]
            cell = cell[cell>np.percentile(cell, 25)]
            if len(np.unique(cell))==0:
                th = 0
            elif len(np.unique(cell))==1:
                th = cell[0]
            else:
                th = skimage.filters.threshold_otsu(cell)
            nucleus = (cellmask*LA)>th
            cells += float(i)*cellmask
            nuclei += float(i)*nucleus
    else:
        nuclei, _ = findcells(imnuc, cellcolormask=cellcolormask, ccdist=ccdist, threshold=thresholdnuc, removeborders=False)
        p = np.where(np.isfinite(im))
        q = np.where(~np.isfinite(im))
        LB = fill_nan(im)
        #LB -= scipy.ndimage.gaussian_filter(LB, 256, mode='nearest')
        LA = scipy.ndimage.gaussian_filter(LB, 2.5, mode='nearest')
        #th = thresholdnuc or skimage.filters.threshold_otsu(LA[p[0], p[1]])
        th = thresholdnuc or otsu_local(LA)
        #mask = np.zeros(im.shape)
        #mask[LA>th] = 1
        #mask[q[0], q[1]] = 0
        mask = skimage.morphology.binary_closing(LA>th, disk(5))
#        mask = skimage.morphology.binary_dilation(mask, disk(5))
        mask = skimage.morphology.binary_erosion(mask, disk(6))
        WS = skimage.segmentation.watershed(-LA, markers=nuclei, mask=mask)
        
        lbl = list(np.unique(WS))
        if 0 in lbl:
            lbl.remove(0)
        cells = np.zeros(im.shape)
        for i in tqdm(lbl, disable=len(lbl)<25):
            cellmask = scipy.ndimage.morphology.binary_fill_holes((WS==i))
            cells += float(i)*cellmask

    if removeborders:
        for lbl in np.unique(np.hstack((cells[:,0], cells[:,-1], cells[0,:], cells[-1,:]))):
            cells[cells==lbl] = 0
            nuclei[nuclei==lbl] = 0
    return cells, nuclei
