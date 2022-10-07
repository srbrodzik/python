# This function is used to remove files from ex_data_v05 for H03 that do not fall in proper domain

# Imports and functions
import os
import netCDF4 as nc4
import numpy as np
import shutil

# User inputs
indir = '/home/disk/bob/gpm/h03_ku/classify/ex_data_v05'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
years = ['2014','2015','2016','2017']

# H03 domain
minLatitude = -67
maxLatitude = -35
minLongitude = -30
maxLongitude = 75

for iyear in years:
    for imonth in months:
        # Create outOfBounds directory if it doesn't exist
        if not os.path.exists(indir+'/'+iyear+'/'+imonth+'/outOfBounds'):
            os.makedirs(indir+'/'+iyear+'/'+imonth+'/outOfBounds')
            
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
            if fname.endswith('nc'):

                print fname

                # open and read lat and lon from input file
                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'r')
                lat = np.array(ncid.variables['lat'])
                #lon = np.array(ncid.variables['lon'])
                ncid.close()

                # find min and max lat and lon
                minLat = np.amin(lat)
                #minLon = np.amin(lon)
                maxLat = np.amax(lat)
                #maxLon = np.amax(lon)
                
                if minLat > maxLatitude:
                    fullpath_infile  = indir+'/'+iyear+'/'+imonth+'/'+fname
                    fullpath_outfile = indir+'/'+iyear+'/'+imonth+'/outOfBounds/'+fname
                    shutil.move(fullpath_infile,fullpath_outfile)
                    
