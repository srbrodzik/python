#Author: Stacy Brodzik
#Date: 20160729
#Description: This code makes correctons to the D3R cfradial files created by
#   RadxConvert.  1. The fixed_angle values for all scans may need correcting
#   and 2. The sweep_mode for the rhis needs to be changed from
#   'azimuth_surveillance' to 'rhi'

import netCDF4 as nc4
import numpy as np
import os
import logging as log

inDir = '/home/disk/bob/olympex/cfradial/moments/d3r_ka_qc/sur/20160115'
swpType = 'sur'

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for fname in os.listdir(inDir):

    if os.path.isfile(inDir+'/'+fname):

        log.info( "file = {}".format(fname) )

        #Open file
        dset = nc4.Dataset(inDir+'/'+fname,'a')

        #Rename sur files with numSwps != 3
        numSwps =  len(dset.dimensions['sweep'])
        log.info( "   numSwps = {}".format(numSwps) )
        if swpType == 'sur' and numSwps != 3:
            dset.close()
            if not os.path.exists(inDir+'/rhiScans'):
                os.makedirs(inDir+'/rhiScans')
            os.rename(inDir+'/'+fname, inDir+'/rhiScans/'+fname)
        #Assign new values for 'fixed_angle' and 'sweep_mode'
        else:
            if swpType == 'sur':
                angles = dset.variables['elevation'][:]
            elif swpType == 'rhi':
                angles = dset.variables['azimuth'][:]
                dset.variables['sweep_mode'][:, 0] = 'r'
                dset.variables['sweep_mode'][:, 1] = 'h'
                dset.variables['sweep_mode'][:, 2] = 'i'
                dset.variables['sweep_mode'][:, 3:] = ''

            #Correct fixed_angle
            dset.variables['fixed_angle'][:] = np.unique(angles)

            #Close file
            dset.close()

