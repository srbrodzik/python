"""
Rain-type Classification code of Powell et al. (2016, JTECH)
Author: Scott Powell and Stacy Brodzik, University of Washington
Date: April 11, 2016

Amendments by
Stacy Brodzik (brodzik@uw.edu)
rev 3 5/2/2016 - added missing value mask to match input reflectivity file

""" 

from __future__ import division  #For python2 only. Alternatively, run interpreter with -Q flag. 
import numpy as np
import os
from rel_raintype_srb import algorithm as alg
from rel_raintype_srb import rtfunctions as rtf
import math
#from rel_raintype_srb import netcdf_io as net
import logging as log
import sys
import matplotlib.pyplot as plt

def raintype(fname, fileDir, refl=None, refl_missing_val=-9999, refl_xdim=620, refl_dx=0.75,
             minZdiff=10, deepcoszero=64, truncZconvthres=40, dBZformaxconvradius=40,
             weakechothres=5, backgrndradius=11, maxConvRadius=5, minradius=12, mindBZuse=-50):
  #Values above are default values. Code will use what is designated in runraintype.py
  #first and will use these only if no value(s) is/are given.

  """
  Description: This is the driver code for the Steiner et al. (1995) convective/stratiform 
  classification code for use with Cartesian gridded datasets. For details, see
  Steiner, M., R. A. Houze, Jr., and S. E. Yuter, 1995: Climatological characterization of 
  three-dimensional storm structure from operational radar and rain gauge data. J. Appl. Meteor., 
  34, 1978-2007.  (http://www.atmos.washington.edu/MG/PDFs/JAM95_stei_climatological.pdf).

  Inputs:
  refl = Reflectivity
  refl_missing_val = missing value in reflectivity data
  refl_dx (km) = horizontal spacing of grid
  minZdiff = factor for comparing echo to background reflectivity; see equation (1) in journal 
     article referenced above
  deepcoszero = see equation (1) in journal article referenced above
  truncZconvthres = reflectivity threshold at or above which echos are classified as convective;  
     The value in Powell et al. (2016) was used for an S-band radar with a beam width of 0.91 degrees. 
     For C-band and/or larger beam width, this value will probably need to be decreased.  Rain
     type classification is most sensitive to this input.
  weakechothres = minimum dBZ for classification as not weak echo; don't change this without a good 
     reason.  7 is about as low as we can go without getting into Bragg scatter territory.
  backgrndradius (km) = radius within which background reflectivity is computed
  maxConvRadius (km) = maximum radius around convective core for possible uncertain classification; 
     Powell et al. (2016) tested 5, and showed that too much convection was included 
     in stratiform region.  Don't lower this number without a good reason.
  minradius = minimum radius for consf output (in km)
  mindBZuse = minimum dBZ for echo

  Outputs:
  rain_type = rain type classification
  types = dict of rain types and their values
  
  """
  global bgmask
  global maskcell

  ## *****************  BEGIN OUTPUT CONSTANTS *****************
  
  # Output constants: Do not change these without a good reason!
  types = {'NO_SFC_ECHO':0,'STRATIFORM':1,'CONVECTIVE':2,'WEAK_ECHO':3,'CS_CORE':4}
  
  ## ***************** END OUTPUT CONSTANTS   ******************

  log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

  #Check inputs
  if refl is None:
    log.info( "No reflectivity provided to raintype module.  Exiting ..." )
    rtout = None
    return rtout, types
  
  #Create missing value mask and turn refl missing values into NaN's
  mask = np.zeros(refl.shape)
  mask[(refl == refl_missing_val)] = 1
  mask[(np.isnan(refl))] = 1
  refl[(mask == 1)] = np.nan

  #Convert dBZ to Z
  Z = rtf.DBZtoZ(refl)

  #If this is the first file in a batch, create a background reflectivity mask and
  #a mask for identifying regions around convective cores.
  if fname == os.listdir(fileDir)[0]:
    bgmask = rtf.makebgmask(backgrndradius,refl_dx)
    maskcell = rtf.makeconvmask(maxConvRadius,refl_dx)
    print maskcell.shape
    # plot array
    #fig = plt.figure()
    #plt.imshow(maskcell[:,:,4])
    #plt.show()
    sys.exit()

  #Now determine the background reflectivity at each grid point.
  background = rtf.get_background_refl(Z,bgmask)

  #Convert Z to dBZ
  background = rtf.ZtoDBZ(background)

  #Run convectivecore.
  rtout = alg.convectivecore(background,
                             refl,
                             types,
                             minZdiff,
                             deepcoszero,
                             truncZconvthres,
                             dBZformaxconvradius,
                             weakechothres,
                             backgrndradius,
                             maxConvRadius,
                             minradius,
                             mindBZuse,
                             refl_dx,
                             maskcell)

  #Apply missing value mask to raintype array
  rtout[mask == 1] = refl_missing_val

  return rtout, types

