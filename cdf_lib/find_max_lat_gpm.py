import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt
import os

indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
years = ['2014']
months = ['03','04']

absolute_max = 0.

for yr in years:
    for mo in months:
        for file in os.listdir(indir+'/'+yr+'/'+mo):
            if file.endswith('nc'):
                
                # open compressed netcdf4 file and read lats
                ncid = nc4.Dataset(indir+'/'+yr+'/'+mo+'/'+file,'r')
                lats = ncid.variables['lat'][:]
                ncid.close()

                # find max lat and print
                max_lat = np.amax(lats)
                if max_lat > absolute_max:
                    absolute_max = max_lat
                if max_lat > 60.0:
                    print file, max_lat

print 'absolute_max=',absolute_max

