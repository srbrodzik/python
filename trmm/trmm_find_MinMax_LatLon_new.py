#!/usr/bin/python

# Imports and functions
import os
import netCDF4 as nc4
import numpy as np

# User inputs
indir = '/home/disk/bob/trmm/afc_pr/classify/ex_data_v07'
months = ['01']
years = ['1998']

minLatitude = 99
maxLatitude = -99
minLongitude = 999
maxLongitude = -999

for iyear in years:
    for imonth in months:
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
        #for fname in os.listdir(indir+'/'+iyear+'/'+imonth+'/outOfBounds'):
            if fname.endswith('nc'):

                print fname

                # open and read lat and lon from input file
                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'r')
                #ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/outOfBounds/'+fname,'r')
                lat = np.array(ncid.variables['lat'])
                lon = np.array(ncid.variables['lon'])
                ncid.close()

                # find min and max lat and lon
                minLat = np.amin(lat)
                minLon = np.amin(lon)
                maxLat = np.amax(lat)
                maxLon = np.amax(lon)
                
                if minLat < minLatitude:
                    minLatitude = minLat
                if minLon < minLongitude:
                    minLongitude = minLon
                if maxLat > maxLatitude:
                    maxLatitude = maxLat
                if maxLon > maxLongitude:
                    maxLongitude = maxLon

print 'min/max lat = ',minLatitude,'/',maxLatitude
print 'min/max lon = ',minLongitude,'/',maxLongitude

