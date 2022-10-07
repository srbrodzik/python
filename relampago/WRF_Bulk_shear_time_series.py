#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Feb 16 16:05:05 2022

@author: crs326

Calculate average bulk shear over SDC region and saves to file

"""

import matplotlib
matplotlib.use('Agg') 

import csv
import math
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.dates as mdates
#import metpy.calc
#from metpy.units import units
from netCDF4 import Dataset
import numpy as np
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import numpy.ma as ma
import os
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, vinterp, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)
import xarray

plt.rcParams.update({'font.size': 40})

# ----------------- reading in file/plotting terrain -----------------

# Open the NetCDF file
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

#change the start and end days to be plotted
start_date = 20181101
end_date = 20181116

########### for finding LLJs ###########

## 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
#crit = [2 ,3200, 5700, 23.326, 11.663]
#
## min and max pressure to plot
#min_pres = 450
#max_pres = 975
#
#crit_num = crit[0]
#max_search_hgt = crit[1]
#min_search_hgt = crit[2]
#max_wind_threshold = crit[3]
#decrease_to_min_threshold = crit[4]

######################################

# get dates from times to only run for correct day folders
input_start_day = pd.to_datetime(str(start_date), format='%Y%m%d', errors='ignore')
input_end_day = pd.to_datetime(str(end_date), format='%Y%m%d', errors='ignore')

# Set the regions for SALLJ 

area_names = ['main_SDC']

areas_lat_bottom_left = [-32.0, -32.0, -30.0, -34.0]
areas_lon_bottom_left = [-64.0, -67.0, -64.0, -64.0]

areas_lat_top_right = [-30.0, -30.0, -28.0, -32.0]
areas_lon_top_right = [-62.0, -65.0, -62.0, -62.0]

area_names = ['COR']
lat_s = -31.298
lon_s = -64.212

#area_names = ['north_SDC', 'east', 'south_SDC', 'west']
#
#areas_lat_bottom_left = [-30.0, -32.0, -34.0, -32.0]
#areas_lon_bottom_left = [-64.0, -61.0, -64.0, -67.0]
#
#areas_lat_top_right = [-28.0,  -30.0, -32.0, -30.0]
#areas_lon_top_right = [-62.0, -59.0, -62.0, -65.0]

# loop through regions to find SALLJ coverage
for i in range(len(area_names)):
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(100,20))
    
    print('area', area_names[i])

    time_plot = []
    bulk_shear_0_1km_list = []

    # loop through files (times)
    for folder_date in sorted(os.listdir(path_wrf)):

        # convert the folder name to datetime object
        folder_date_dt = pd.to_datetime(folder_date, format='%Y%m%d', errors='ignore')

        # only check folder if it is named with a date
        if folder_date.startswith('20') and (folder_date_dt >= input_start_day) and (folder_date_dt <= input_end_day):

            folder_path = os.path.join(path_wrf, folder_date)

            sorted_hour_files = sorted(os.listdir(folder_path))

            # go through hourly files for the chosen date
            for hourly_file in sorted_hour_files:
                
                if hourly_file.startswith('wrfout'):

                    # for each new hour set count of SALLJ points to 0
                    count_points_SALLJ_hour = 0
                    count_total_hour = 0

                    file_path = os.path.join(folder_path, hourly_file)

                    file_hour = int(hourly_file[22:24])

                    # get the netCDF
                    ncfile = Dataset(file_path,'r')
                    file_name = hourly_file

                    # get xy for regions
#                    bottom_left_xy = ll_to_xy(ncfile, areas_lat_bottom_left[i], areas_lon_bottom_left[i])
#                    top_right_xy = ll_to_xy(ncfile, areas_lat_top_right[i], areas_lon_top_right[i])  

                    print('file time', str(file_name[16:24]))

                    # get xy for station
                    Station_xy = ll_to_xy(ncfile, lat_s, lon_s)

                    # read in the variables 
                    pres = getvar(ncfile, 'pressure')[:,Station_xy[1],Station_xy[0]]
                    u = getvar(ncfile, 'ua', units='kt')[:,Station_xy[1],Station_xy[0]] # in kts
                    v = getvar(ncfile, 'va', units='kt')[:,Station_xy[1],Station_xy[0]] # in kts
                    hght = getvar(ncfile, 'height', msl=False, units='m')[:,Station_xy[1],Station_xy[0]] # in m NOTE: this is AGL not MSL!!
                    speed, drct = getvar(ncfile, 'wspd_wdir', units='kt') # in kts
                    
                    speed = speed[:,Station_xy[1],Station_xy[0]]
                    drct = drct[:,Station_xy[1],Station_xy[0]]
                    
                    # calculate bulk shear
                    #bulk_shear_0_1km = metpy.calc.bulk_shear(pres, u, v, depth=1000 * units.meter)
                    
                    # now some shear calculations
                    wspd1km = speed[np.argmin(abs(hght-1000))]
                    #print('wspd1km', wspd1km)
                    wdir1km = drct[np.argmin(abs(hght-1000))]
                    #print('wdir1km', wdir1km)

                    udiff = wspd1km*np.cos(np.radians(270-wdir1km)) - speed[3]*np.cos(np.radians(270-drct[3]))
                    vdiff = wspd1km*np.sin(np.radians(270-wdir1km)) - speed[3]*np.sin(np.radians(270-drct[3]))
                    #print('udiff', udiff)
                    #print('vdiff', vdiff)
                    # print udiff, vdiff
                    shear1km = np.sqrt(udiff**2 + vdiff**2)
                    if (math.isnan(shear1km)):
                        shear1km = 0

                    print('shear1km', shear1km)
                    
                    # append bulk shear to list
                    bulk_shear_0_1km_list.append(shear1km)
                    #Average_MUCAPE_subset_list.append(np.nanmean(MCAPE_subset))
                    time_plot.append(pd.to_datetime(hourly_file[11:24], format='%Y-%m-%d_%H', errors='ignore'))

                    #print('np.nanmean(MCAPE_subset)', np.nanmean(MCAPE_subset))

                    ncfile.close()
                
    # write data to .csv file
    with open('PNNL_WRF_SDC_bulk_shear_0_1km_time_series_%s_%s_%s.csv' % (area_names[i], str(start_date)[4:], str(end_date)[4:]), 'w') as f:
        writer = csv.writer(f)

        # write the times to a row
        writer.writerow(time_plot)

        # write the data to the next rows
        writer.writerow(bulk_shear_0_1km_list)

print('done')
