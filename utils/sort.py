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

inDir = '/home/disk/bob/olympex/cfradial/moments/d3r_ka_qc/sur/20151113'

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for fname in os.listdir(inDir):

    if os.path.isfile(inDir+'/'+fname):

        log.info( "file = {}".format(fname) )

        #Open file
        dset = nc4.Dataset(inDir+'/'+fname,'a')

        numSwps =  len(dset.dimensions['sweep'])
        log.info( "   numSwps = {}".format(numSwps) )
        if numSwps != 3:
            dset.close()
            if not os.path.exists(inDir+'/rhiScans'):
                os.makedirs(inDir+'/rhiScans')
            os.rename(inDir+'/'+fname, inDir+'/rhiScans/'+fname)
