#Author: Stacy Brodzik
#Date: 20160729
#Description: This code makes correctons to the D3R cfradial files created by
#   RadxConvert.  1. The fixed_angle values for all scans may need correcting
#   and 2. The sweep_mode for the rhis needs to be changed from
#   'azimuth_surveillance' to 'rhi'

#import netCDF4 as nc4
import pyart.io.cfradial as pyio
import numpy as np
#import numpy.ma as ma
import os
import logging as log

inDir = '/home/disk/bob/olympex/cfradial/moments/d3r_ka_qc/sur/20151108'
swpType = 'sur'

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for fname in os.listdir(inDir):

    if os.path.isfile(inDir+'/'+fname):

        log.info( "file = {}".format(fname) )

        #Read file
        radar = pyio.read_cfradial(inDir+'/'+fname)
        #radar = pyio.read_cfradial(inDir+'/'+fname,field_names=["fixed_angle","azimuth","elevation","sweep_mode"])
    
        #Correct fixed_angle
        if swpType == 'rhi':
            angles = radar.azimuth["data"]
        elif swpType == 'sur':
            angles = radar.elevation["data"]
            
        radar.fixed_angle["data"] = np.unique(angles)

        #Change sweep_mode
        if swpType == 'rhi':
            numMask = len(radar.sweep_mode["data"].mask)
            radar.sweep_mode["data"][:,0] = 'r'
            radar.sweep_mode["data"][:,1] = 'h'
            radar.sweep_mode["data"][:,2] = 'i'
            radar.sweep_mode["data"].mask[:,3:numMask-1] = True
    
        #Write fixed_angle and sweep_mode back to original file
        pyio.write_cfradial(inDir+'/'+fname,radar)

