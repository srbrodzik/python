#Author: Stacy Brodzik
#Date: 20161101
#Description: This code separates the D3R rhi files into east and west facing
#sectors. We need this for two reasons: 1. To get east-sector rhi tick marks
#to show in CIDD and 2. To create east and west facing cfads

import netCDF4 as nc4
import numpy as np
import os
import logging as log

inDir = '/home/disk/bob/olympex/cfradial/moments/d3r_ku_qc/rhi/20151110'

#Make east and west subdirs if they don't exist
eastDir = inDir+'/east'
westDir = inDir+'/west'
if not os.path.exists(eastDir):
    os.makedirs(eastDir)
if not os.path.exists(westDir):
    os.makedirs(westDir)

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for fname in os.listdir(inDir):

    if os.path.isfile(inDir+'/'+fname):

        log.info( "file = {}".format(fname) )

        #Read file
        dset = nc4.Dataset(inDir+'/'+fname,'a')
        
        #Get min azimuth
        #numSweeps = dset.dimensions['sweep']
        angles = dset.variables['azimuth'][:]
        minAngle = min(angles)

        #Separate into sectors
        if minAngle < 180.0:
            os.rename(inDir+'/'+fname, eastDir+'/'+fname)
        else:
            os.rename(inDir+'/'+fname, westDir+'/'+fname)
        
        #Close file
        dset.close()

