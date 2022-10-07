# This program will interpolate the reflectivity files that getTRMM found to be part of
# a storm system into a regular grid defined by gridLatitude-gridLongitude and height.
# The interpolated data will be output in zebra netcdf format.

from __future__ import absolute_import
from __future__ import division
import netCDF4 as nc4
import os
import numpy as np
import logging as log
import math
import glob
import trmm_pr_netcdf_io as net
#import subprocess

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
binDir = '/home/disk/stellar2/brodzik/interp/src'
interpExe = binDir + '/' + 'PRInterp'
dataBaseDir = '/home/disk/bob/trmm_v7/TEST/ascii'

## years and months to process
#months=['01','02','03','04','05','06','07','08','09','10','11','12']
months=['11']
years=['2007']

## set interpolation parameters
radius = 4.25
levels = 80
delta_z = 0.250
missingVal = -99

## info for netcdf output file
title = 'Interpolated GPM 2Ku Data'
institution = 'University of Washington'
# ------------------------------------- END INPUTS -------------------------------------

# set up logging
log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for i,iyear in enumerate(years):
    for j,imonth in enumerate(months):
        dataDir = dataBaseDir + '/' + iyear + '/' + imonth
        log.info('dataDir = {}'.format(dataDir) )
        os.chdir(dataDir)
        for file in os.listdir(dataDir):
            if file.endswith('gridColumns.txt'):
                prefix = file[0:27]
                log.info('prefix = {}'.format(prefix) )
                
                ## run interpolation
                cmd = interpExe + ' ' + prefix + ' ' + str(radius) + ' ' + str(levels)
                os.system(cmd)
                
                ## read ascii files
                (baseTime,lats,lons,alts,raintype,swath,refl) = net.readTRMMascii(dataDir,prefix,levels,
                                                                                  delta_z,missingVal)

                ## output netcdf file
                net.writeZeb_netcdf(dataDir,prefix,baseTime,lats,lons,alts,raintype,swath,refl,
                                    title,institution,missingVal)
                log.info('   Done writing zeb netcdf file')

                ## remove all ascii files
                for textFile in glob.glob(prefix+'*.txt'):
                    os.remove(textFile)
                              
