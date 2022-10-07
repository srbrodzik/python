from __future__ import division     #For python2 users only.
import numpy as np
import math

def ZtoDBZ(z):
  return 10*np.log10(z)

def DBZtoZ(dbz):
  return 10**(0.1*dbz)

def radialDistanceMask(minRadius,maxRadius,xdim,dx):
  mask = np.zeros((xdim, xdim), dtype='float')
  centerPixelX = int(np.floor(xdim/2))
  centerPixelY = int(np.floor(xdim/2))
  for i in range(0,xdim):
    for j in range(0,xdim):
      x_range_sq=((centerPixelX-i)*dx) * ((centerPixelX-i)*dx)
      y_range_sq=((centerPixelY-j)*dx) * ((centerPixelY-j)*dx)
      r = np.sqrt(x_range_sq+y_range_sq)
      if np.logical_and(r<=maxRadius,r>=minRadius):
        mask[i,j] = 1
  return mask

def makebgmask(backgrndradius,dx):

  global bgmask

  bgrange = int( np.floor((float(backgrndradius)/dx)*2) )
  if np.floor(bgrange/2) == bgrange/2.0:
    bgrange = bgrange + 1
  print bgrange

  bgmask = radialDistanceMask(0,backgrndradius,bgrange,dx)
  
  return bgmask

def makeconvmask(maxConvRadius,dx):

  print 'New version'
  
  global maskcell

  maxConvDiameter = int( np.floor((float(maxConvRadius)/dx)*2) )
  if np.floor(maxConvDiameter/2) == maxConvDiameter/2.0:
    maxConvDiameter = maxConvDiameter + 1
  print "maxConvDiameter = ",maxConvDiameter
    
  maskcell = np.zeros([maxConvDiameter,maxConvDiameter,maxConvRadius])
  
  for iradius in range(0,maxConvRadius):
    mask = np.zeros([maxConvDiameter,maxConvDiameter])
    mask = radialDistanceMask(0,iradius+1,maxConvDiameter,dx)
    maskcell[:,:,iradius] = mask

  return maskcell


def get_background_refl(Z,bgmask):

  from scipy import signal as sg

  #Create background variable.
  background = np.zeros(Z.shape)
  Z[np.isnan(Z)] = 0

  #First convolve the mask with the reflectivity data. The result will be
  #a matrix that reports background reflectivity even where there is zero
  #reflectivity.
  bg1 = sg.convolve2d(Z,bgmask,mode='same',boundary='symm',fillvalue=0)
  Z[Z != 0] = 1
  #Next, convolve the mask with ones and zeros, where the ones are where
  #reflectivity is nonzero.
  bg2 = sg.convolve2d(Z,bgmask,mode='same',boundary='symm',fillvalue=0)
  #Normalize to determine the actual background reflectivity.
  background = bg1/bg2
  #Make sure non-existent values in background are NaN
  background[Z==0] = np.nan

  return background

def chopmask(inmask,topchop,rightchop,btmchop,leftchop):

  #Simply trims inmask so it fits data on edge of domain. Returns newmask. 
  newmask = 0
  masksize = inmask.shape[1]

  if topchop != 0:
      newmask = inmask[:,0:masksize-topchop] 
  if rightchop != 0:
      newmask = inmask[0:masksize-rightchop,:]
  if btmchop != 0:
      newmask = inmask[:,btmchop:]
  if leftchop != 0:
      newmask = inmask[leftchop:,:]

  return newmask

def makedBZcluster(refl,isCore,convsfmat,weakechothres,minsize,maxsize,startslope,
                   shallowconvmin,truncZconvthres,types,dx):

  from scipy import ndimage as nd

  #Allocate matrix indicating whether rain is occurring.
  rain = np.zeros((refl.shape),dtype=np.int)

  #If echo is strong enough, rain = 1.
  rain[refl>=weakechothres] = 1

  #Create truncvalue, which has same shape as reflectivity data. It indicates the 
  #reflectivity over which an echo is automatically classified as some sort of 
  #ISOLATED CONVECTIVE echo. See details below.
  truncvalue = np.ones((refl.shape),dtype=np.float64)*truncZconvthres

  #This is a blob detector. Detects contiguous areas of raining pixels. Diagonally
  #touching pixels that share a corner don't count. Edges must touch.
  #echoes contains the blob objects, numechoes is just a count of them.
  (echoes,numechoes) = nd.label(rain)

  for i in range(0,numechoes):
    #Find 2D indices of echo object.
    (I,J) = (rain*(echoes==i+1)==1).nonzero()
    
    #Compute the total areal coverage of the echo object.
    clusterarea = (dx**2)*len(I)    #In km^2

    #Any echo object with a size between minsize and maxsize is considered 
    #ISOLATED CONVECTION. First, make all of it FRINGE.
    if clusterarea >= minsize and clusterarea <= maxsize:
      convsfmat[I,J] = types['ISO_CONV_FRINGE'] 

    #Very small echo objects are dismissed as WEAK ECHO.  
    if clusterarea < minsize:
      isCore[I,J] = 0
      convsfmat[I,J] = types['WEAK_ECHO']
    #Echo objects with size between minsize and startslope get a small truncvalue
    #equal to shallowconvmin.
    elif clusterarea >= minsize and clusterarea < startslope:
      truncvalue[I,J] = shallowconvmin
    #Echo objects with size between startslope and maxsize get a truncvalue that 
    #is linearly interpolated between shallowconvmin and truncZconvthres depending
    #on the size relative to startslope and maxsize.
    elif clusterarea >= startslope and clusterarea <= maxsize:
      truncvalue[I,J] = shallowconvmin + float((clusterarea-startslope)/(maxsize-startslope))*(truncZconvthres-shallowconvmin)

    #Evaluate isCore with size of echo object accounted for.
    #First, if reflectivity exceeds truncvalue, classify it as ISOLATED CONVECTIVE CORE.
    isCore[refl >= truncvalue] = types['ISO_CS_CORE']

    #But if reflectivity exceeds original reflectivity threshold, classify as CONVECTIVE core.
    isCore[refl >= truncZconvthres] = types['CS_CORE']

  return (convsfmat,isCore)


def radial_distance_mask(min_radius, max_radius, xdim, ydim, x_spacing, y_spacing):

    """
    Description: Creates a radar coverage mask for a square grid

    Inputs:

    Outputs:

    """
    mask = np.zeros(shape=(ydim,xdim))

    center_y = int(math.floor(ydim/2.)) - 0.5
    center_x = int(math.floor(xdim/2.)) - 0.5
      
    for j in range(0,ydim):
        for i in range(0,xdim):
            y_range_sq = math.pow( ((center_y-j)*y_spacing), 2 )
            x_range_sq = math.pow( ((center_x-i)*x_spacing), 2 )
            dist = math.sqrt(x_range_sq+y_range_sq)
            if dist < max_radius and dist > min_radius:
                mask[j,i] = 1

    return mask
