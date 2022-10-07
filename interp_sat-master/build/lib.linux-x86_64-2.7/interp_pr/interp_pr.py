# This program will interpolate the reflectivity files that getTRMM found to be part of
# a storm system into a regular grid defined by gridLatitude-gridLongitude and height.
# The interpolated data will be output in zebra netcdf format.

from __future__ import absolute_import
from __future__ import division
import netCDF4 as nc4
import sys
import os
import numpy as np
import logging as log
import math
import glob
import ascii_netcdf_io as net
#import subprocess

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
binDir = '/home/disk/stellar2/brodzik/interp/src'
#binDir = '/data/keeling/a/atmos-brodzik/bin'
interpExe = binDir + '/' + 'PRInterp_varRays'
#dataBaseDir = '/home/disk/bob/trmm_v7/TEST/ascii'
#dataBaseDir = '/data/keeling/a/atmos-brodzik/data/trmm'
dataBaseDir = '/home/disk/bob/gpm/m2_cmed_ku/classify/ex_data'

## years and months to process
#months=['01','02','03','04','05','06','07','08','09','10','11','12']
months=['05']
years=['2014']

## set interpolation parameters
radius = 4.25
missingVal = -99

## set satellite - must be either 'trmm' or 'gpm'
satellite = 'gpm'

## info for netcdf output file
title = 'Interpolated GPM 2Ku Data'
institution = 'University of Washington'
# ------------------------------------- END INPUTS -------------------------------------

# set up logging
log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# set params based on satellite type
if satellite == 'trmm':
    levels = 80
    delta_z = 0.250
    nrays = 49
elif satellite == 'gpm':
    levels = 176
    delta_z = 0.125
    nrays = 49
else:
    log.info('satellite = {} - must be either trmm or gpm.  Exiting.'.format(satellite) )
    sys.exit()

for i,iyear in enumerate(years):
    for j,imonth in enumerate(months):
        dataDir = dataBaseDir + '/' + iyear + '/' + imonth
        log.info('dataDir = {}'.format(dataDir) )
        os.chdir(dataDir)
        for file in os.listdir(dataDir):
            if file.endswith('gridColumns.txt'):
                if satellite == 'trmm':
                    prefix = file[0:27]
                else:
                    prefix = file[0:31]
                    
                log.info('prefix = {}'.format(prefix) )
                
                ## run interpolation
                cmd = interpExe + ' ' + prefix + ' ' + str(radius) + ' ' + str(levels) + ' ' + str(nrays)
                os.system(cmd)
                
                ## read ascii files
                if satellite == 'trmm':
                    (baseTime,lats,lons,alts,raintype,swath,refl) = net.readTRMMascii(dataDir,prefix,levels,
                                                                                         delta_z,missingVal)
                else:
                    (baseTime,lats,lons,alts,rainType,rainTypeRaw,phaseType,phaseTypeRaw,shallowRainType,
                     nearSurfRain,bbWidth,bbHeight,swath,refl,hdfFileName) = net.readGPMascii(dataDir,prefix,levels,
                                                                                              delta_z,missingVal)

                log.info('Done readying ascii files.')

                ## output netcdf file
                if satellite == 'trmm':
                    net.writeZeb_TRMMnetcdf(dataDir,prefix,baseTime,lats,lons,alts,rainType,swath,refl,
                                            title,institution,missingVal)
                else:
                    net.writeZeb_GPMnetcdf(dataDir,prefix,baseTime,lats,lons,alts,rainType,rainTypeRaw,
                                           phaseType,phaseTypeRaw,shallowRainType,nearSurfRain,bbWidth,
                                           bbHeight,swath,refl,hdfFileName,title,institution,missingVal)
                                            
                log.info('   Done writing zeb netcdf file')

                ## remove all ascii files
                #for textFile in glob.glob(prefix+'*.txt'):
                #    os.remove(textFile)
                              
