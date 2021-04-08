##################################################################
##
## Spot-fitting and bandpass routines
##
## Antoine Coulon, 2020
## Contact: <antoine.coulon@curie.fr>
##
## Code under GPL v3 license

from scipy import *
from scipy import ndimage
from skimage import filters
from misc import *

def GaussianMaskFit2(im,coo,s,optLoc=1,bgSub=2,winSize=13,convDelta=.01,nbIter=20):
  """Applies the algorithm from [Thompson et al. (2002) PNAS, 82:2775].
Parameters:
- im: a numpy array with the image
- coo: approximate coordinates (in pixels) of the spot to localize and measure. Note, the coordinates are x,y!
- s: width of the PSF in pixels
- s: width of the PSF in pixels
- optLoc: If 1, applied the iterative localization refinement algorithm, starting with the coordinates provided in coo. If 0, only measures the spot intensity at the coordinates provided in coo.
- bgSub: 0 -> no background subtraction. 1 -> constant background subtraction. 2 -> tilted plane background subtraction.
- winSize: Size of the window (in pixels) around the position in coo, used for the iterative localization and for the background subtraction.
- convDelta: cutoff to determine convergence, i.e. the distance (in pixels) between two iterations
- nbIter: the maximal number of iterations.

Returns
- the intensity value of the spot.
- the corrdinates of the spot.

If convergence is not found after nbIter iterations, return 0 for both intensity value and coordinates.
"""
  from scipy import optimize; coo=array(coo);
  for i in range(nbIter):
    if not prod(coo-winSize/2.>=0)*prod(coo+winSize/2.<=im.shape[::-1]): return 0.,r_[0.,0.], 0.
    winOrig=(coo-int(winSize)/2).astype(int)
    i,j=meshgrid(winOrig[0]+r_[:winSize],winOrig[1]+r_[:winSize]);
    N=exp(-(i-coo[0])**2/(2*s**2)-(j-coo[1])**2/(2*s**2))/(2*pi*s**2)
    S=im[:,winOrig[0]:winOrig[0]+winSize][winOrig[1]:winOrig[1]+winSize]*1.
    if bgSub==2:
      xy=r_[:2*winSize]%winSize-(winSize-1)/2.
      bgx=polyfit(xy,r_[S[0],S[-1]],1); S=(S-xy[:winSize]*bgx[0]).T;
      bgy=polyfit(xy,r_[S[0],S[-1]],1); S=(S-xy[:winSize]*bgy[0]).T;
      bg=mean([S[0],S[-1],S[:,0],S[:,-1],]); S-=bg
      bg=r_[bg,bgx[0],bgy[0]]
    if bgSub==1:
      bg=mean([S[0],S[-1],S[:,0],S[:,-1],]); S-=bg
    S=S.clip(0) # Prevent negative values !!!!
    if optLoc:
      SN=S*N; ncoo=r_[sum(i*SN),sum(j*SN)]/sum(SN)
      #ncoo=ncoo+ncoo-coo # Extrapolation
      if abs(coo-ncoo).max()<convDelta: return sum(SN)/sum(N**2),coo,bg
      else: coo=ncoo
    else: return sum(S*N)/sum(N**2),coo,bg
  return 0.,r_[0.,0.], 0.


def GaussianMaskFit3D(im,coo,s,sZ=None,optLoc=1,bgSub=2,winSize=None,winSizeZ=None,convDelta=.05,nbIter=20):
  """Applies the algorithm from [Thompson et al. (2002) PNAS, 82:2775] adapted to 3D images.
Parameters:
- im: a numpy array with the image
- coo: approximate z,y,x coordinates (in pixels) of the spot to localize and measure.
- s: width of the PSF in x,y in pixels
- sZ: width of the PSF in z in pixels. Defaults to the same value as s.
- optLoc: If 1, applied the iterative localization refinement algorithm, starting with the coordinates provided in coo. If 0, only measures the spot intensity at the coordinates provided in coo.
- bgSub: 0 -> no background subtraction. 1 -> constant background subtraction. 2 -> tilted plane background subtraction.
- winSize: Size of the x,y window (in pixels) around the position in coo, used for the iterative localization and for the background subtraction.
- winSizeZ: Same as winSize, in the z dimension.
- convDelta: cutoff to determine convergence, i.e. the distance (in pixels) between two iterations
- nbIter: the maximal number of iterations.

Returns
- the intensity value of the spot.
- the coordinates of the spot (z,y,x).
- the background level:
   - either a constant value if bgSub=1
   - or [offset, tilt in z, tilt in y, tilt in x] if bgSub=2

If convergence is not found after nbIter iterations, return 0 for both intensity value and coordinates.
"""
  coo=array(coo);
  if sZ==None: sZ=s
  if winSize ==None: winSize =int(ceil(s*8./2))*2+1
  if winSizeZ==None: winSizeZ=int(ceil(sZ*4./2))*2+1
  for i in range(nbIter):
    if not (winSizeZ/2.<=coo[0]<=im.shape[0]-winSizeZ/2.)*prod([winSize/2.<=coo[j]<=im.shape[j]-winSize/2. for j in [1,2]]):
      return 0.,r_[0.,0.,0.], 0.
    winOrig=r_[coo[0]-int(winSizeZ/2),coo[1:]-int(winSize/2)].astype(int)
    iy,iz,ix=meshgrid(winOrig[1]+r_[:winSize],winOrig[0]+r_[:winSizeZ],winOrig[2]+r_[:winSize]);
    N=exp(-(iz-coo[0])**2/(2*sZ**2)-(iy-coo[1])**2/(2*s**2)-(ix-coo[2])**2/(2*s**2))/((2*pi)**1.5*s*s*sZ)
    S=im[winOrig[0]:winOrig[0]+winSizeZ][:,winOrig[1]:winOrig[1]+winSize][:,:,winOrig[2]:winOrig[2]+winSize]*1.
    if bgSub==2:
      cxy=r_[:winSize]-(winSize-1)/2.
      cz=r_[:winSizeZ]-(winSizeZ-1)/2.
      bgx=polyfit(cxy,mean(r_[S[:,0],S[:,-1]],0),1)[0];
      bgy=polyfit(cxy,mean(r_[S[:,:,0],S[:,:,-1]],0),1)[0];
      bgz=polyfit(cz,mean(c_[S[:,0],S[:,-1],S[:,1:-1,0],S[:,1:-1,-1]],1),1)[0];
      S=rollaxis(rollaxis(rollaxis(S-cxy*bgx,2)-cxy*bgy,2)-cz*bgz,2)
      bg=mean([S[:,0],S[:,-1],S[:,:,0],S[:,:,-1],]); S-=bg
      bg=r_[bg,bgz,bgy,bgx]
    if bgSub==1:
      bg=mean([S[:,0],S[:,-1],S[:,:,0],S[:,:,-1],]); S-=bg
    #S=S.clip(0) # Prevent negative values !!!!
    if optLoc:
      SN=S*N; ncoo=r_[sum(iz*SN),sum(iy*SN),sum(ix*SN)]/sum(SN)
      #ncoo+=ncoo-coo # Extrapolation of localization step !!!!
      #ncoo+=(ncoo-coo)*.7 # Extrapolation of localization step !!!!
      #print(i,ncoo,abs(coo-ncoo).max())
      if abs(coo-ncoo).max()<convDelta: return sum(SN)/sum(N**2),coo,bg
      else: coo=ncoo
    else: return sum(S*N)/sum(N**2),coo,bg
  return 0.,r_[0.,0.,0.], 0.

sHS=fftpack.fftshift # Swap half-spaces. sHS(matrix[, axes]). axes=all by default
def hS(m,axes=None):
    if axes==None: axes=range(rank(m))
    elif type(axes)==int: axes=[axes]
    elif axes==[]: return m
    return hS(m.swapaxes(0,axes[-1])[:m.shape[axes[-1]]/2].swapaxes(0,axes[-1]),axes[:-1])

def sHSM(m,axes=None):
    if axes==None: axes=range(rank(m))
    elif type(axes)==int: axes=[axes]
    m=m.swapaxes(0,axes[0]); max=m[1]+m[-1]; m=(m+max/2)%max-max/2; m=m.swapaxes(0,axes[0])
    return sHS(m,axes)


def bpass(im,r1=1.,r2=1.7):
  ker1x=exp(-(sHS(sHSM(r_[:im.shape[1]]))/r1)**2/2); ker1x/=sum(ker1x); fker1x=fft(ker1x);
  ker1y=exp(-(sHS(sHSM(r_[:im.shape[0]]))/r1)**2/2); ker1y/=sum(ker1y); fker1y=fft(ker1y);
  ker2x=exp(-(sHS(sHSM(r_[:im.shape[1]]))/r2)**2/2); ker2x/=sum(ker2x); fker2x=fft(ker2x);
  ker2y=exp(-(sHS(sHSM(r_[:im.shape[0]]))/r2)**2/2); ker2y/=sum(ker2y); fker2y=fft(ker2y);
  fim=fftpack.fftn(im)
  return fftpack.ifftn((fim*fker1x).T*fker1y-(fim*fker2x).T*fker2y).real.T

def bpass3D(im,r1=1.,r2=1.7,rz1=1.,rz2=1.7,zMirror=False):
    psfPx = r2
    psfPxZ = rz2
    return filters.gaussian(im,1.,mode='mirror')-filters.gaussian(im,r_[psfPxZ,psfPx,psfPx],mode='mirror')


