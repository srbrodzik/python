import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt
import os

indir = '/home/disk/bob/gpm/aka_ku/classify/ex_data_v05'
years = ['2014']
months = ['03']

absolute_max_lat = 0.
absolute_min_lat = 90.
absolute_max_lon = -180.
absolute_min_lon = 180.

for yr in years:
    for mo in months:
        for file in os.listdir(indir+'/'+yr+'/'+mo):
            if file.endswith('nc'):

                print 'file = ',file
                
                # open compressed netcdf4 file and read lats
                ncid = nc4.Dataset(indir+'/'+yr+'/'+mo+'/'+file,'r')
                lats = ncid.variables['lat'][:]
                lons = ncid.variables['lon'][:]
                ncid.close()

                # find min/max lat/lon and print
                min_lat = np.amin(lats)
                if min_lat < absolute_min_lat:
                    absolute_min_lat = min_lat
                max_lat = np.amax(lats)
                if max_lat > absolute_max_lat:
                    absolute_max_lat = max_lat

                min_lon = np.amin(lons)
                if min_lon < absolute_min_lon:
                    absolute_min_lon = min_lon
                max_lon = np.amax(lons)
                if max_lon > absolute_max_lon:
                    absolute_max_lon = max_lon

print 'min lat = ',absolute_min_lat
print 'max lat = ',absolute_max_lat

print 'min lon = ',absolute_min_lon
print 'max lon = ',absolute_max_lon
