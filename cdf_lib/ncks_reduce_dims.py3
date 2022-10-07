#!/usr/bin/python3

import os
import numpy as np
import netCDF4 as nc4
import shutil

indir = '/home/disk/archive3/trmm/uw2/WPO/interp_data'
years = ['1998']
months = ['02','03','04','05','06','07','08','09','10','11','12']

# Test file - ROV lons = 162.8 / 165.45
#file = 'TPR7_uw2_00943.19980126.164449_WPO.nc4'
# Test file - out of range; lons = 146.9 / 154.7
#file = 'TPR7_uw2_00944.19980126.181925_WPO.nc4'
# Test file - crosses 180; ROV lons = [170.8 / -177.55]
#file = 'TPR7_uw2_00969.19980128.075258_WPO.nc4'
#file = 'TPR7_uw2_00564.19980102.150852_WPO.nc4'

# WPO region
minLat = -40.
maxLat = 40.
minLon = 145.
maxLon = 178.

for iyear in years:
    for imonth in months:

        outOfRangeDir = indir+'/'+iyear+'/'+imonth+'/outOfRange'
        if not os.path.exists(outOfRangeDir):
            os.makedirs(outOfRangeDir)

        for file in os.listdir(indir+'/'+iyear+'/'+imonth):
            if file.endswith('nc4'):
                
                print(file)
                file_fullpath = indir+'/'+iyear+'/'+imonth+'/'+file

                # Open file and read lat and lon
                ncid = nc4.Dataset(file_fullpath,'r')
                lats = np.array(ncid.variables['latitude'])
                lons = np.array(ncid.variables['longitude'])
                ncid.close()

                # If necessary, adjust lons so they are in [0,360] range
                if np.min(lons) < 0:
                    cmd = "ncap2 -s 'where(longitude<0) longitude=longitude+360.' "+file_fullpath+" "+file_fullpath+".tmp"
                    os.system(cmd)
                    shutil.move(file_fullpath+'.tmp',file_fullpath)
                
                # check to see if file in WPO area
                indices = np.where(np.logical_and(lons>=minLon, lons<=maxLon))
                if np.size(indices) > 0:
                    print('   in range')
                    cmd = 'ncks -d latitude,'+str(minLat)+','+str(maxLat)+' -d longitude,'+str(minLon)+','+str(maxLon)+' '+file_fullpath+' '+file_fullpath+'.wpo'
                    os.system(cmd)
                else:
                    print('   out of range')
                    shutil.move(indir+'/'+iyear+'/'+imonth+'/'+file,
                                outOfRangeDir+'/'+file)
                    

    
