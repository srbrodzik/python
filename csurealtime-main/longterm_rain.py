# Brody Fuchs, Sept 2017
# brfuchs@atmos.colostate.edu

# This will read in the saved rainfall data and make longer term accumulation plots

import numpy as np 
import matplotlib.pyplot as plt 
from netCDF4 import Dataset
import glob
import matplotlib as mpl
from matplotlib.colors import LogNorm

from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, 
                            csu_kdp, csu_misc, csu_blended_rain_tropical)

from mpl_toolkits.basemap import Basemap 

import argparse
from copy import deepcopy
import datetime
import time
import raintype
from csuram import RadarData, Cell, general_tools, Case, RadarConfig
from skewt import SkewT
import cPickle as pickle
import os
import sys
from matplotlib.ticker import LogFormatter 
plt.switch_backend('agg')
import Config
import logging
import traceback


def file_time_location(time1, time_arr):
    out = None  
    for i in range(len(time_arr)-1):
    	if time1>time_arr[i] and time1<time_arr[i+1]:
    	    out = i

    return out



def plot_point(point, angle, length):
     '''
     point - Tuple (x, y)
     angle - Angle you want your end point at in degrees.
     length - Length of the line you want to plot.

     Will plot the line on a 10 x 10 plot.
     '''

     # unpack the first point
     x, y = point

     # find the end point
     endy = length * np.sin(np.radians(angle))
     endx = length * np.cos(np.radians(angle))

     return x+endx, y+endy


def plot_range_rings(a, rlat, rlon, rads, minor_rads=None):


    for icr in major_circle_rads:


        icr_conv = icr*nm2km/111.0
	dummycircle = plt.Circle((rlon, rlat), icr_conv, facecolor='none', edgecolor='black', 
	    				linewidth=0.5, alpha=0.5, linestyle='dashed', zorder=9)
	a.add_artist(dummycircle)
	if icr_conv < lat_width/2.0:
	    a.text(rlon, rlat-icr_conv, '%d nm'%(icr), ha='center', va='center', alpha=0.5, zorder=10)



    for imcr in minor_rads:

	imcr_conv = imcr*nm2km/111.0
	dummycircle_minor = plt.Circle((rlon, rlat), imcr_conv, facecolor='none', edgecolor='black', 
	    				linewidth=0.5, alpha=0.25, linestyle='dashed', zorder=9)
	a.add_artist(dummycircle_minor)
	if imcr_conv < lat_width/2.0:
	    a.text(rlon, rlat-imcr_conv, '%d nm'%(imcr), ha='center', va='center', alpha=0.25, zorder=10)



def plot_azimuths(a, rlat, rlon, azs, minor_azs=None):
    for az in azs:
	fx, fy = plot_point((rlon, rlat), az, 200.0/111.0)
	a.plot([rlon, fx], [rlat, fy], color='black', linewidth=0.5, alpha=0.5, linestyle='dashed')

    if minor_azs is not None:
    	for maz in minor_azs:
	    mfx, mfy = plot_point((rlon, rlat), maz, 200.0/111.0)
	    a.plot([rlon, mfx], [rlat, mfy], color='black', linewidth=0.5, alpha=0.13, linestyle='dashed')


def convert_decimal_to_degree(val):
    val_string = str(val)
    vs_split = val_string.split('.')
    decimal_string = vs_split[1]
    if len(decimal_string) == 1:
    	decimal_string += '0'
    decimal_float = float(decimal_string)*60.0/100.0
    out_string = "%d$^{\circ}$ %d'"%(int(vs_split[0]), int(decimal_float))
    return out_string


start = time.time()

base_path = '/home/rmet/SPURS2/realtime_radar'
save_path = '%s/saved_data'%(base_path)
file_type = 'npz'

parser = argparse.ArgumentParser(description='Put in a file to be processed')

#parser.add_argument('--noarg', action="store_true", default=False)
parser.add_argument('--file', action="store", dest="file")
#parser.add_argument('--realtime', action="store", dest="realtime", type=bool, default=False)
parser.add_argument('--config', action="store", dest="config", default=None)


pargs = parser.parse_args()

#print pargs


if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v



    #print 'Processing radar files in realtime'
# Check for the accumulated radar times pickled list
load_check = glob.glob('%s/%s/radar_times_accum.npz'%(base_path, ycfg['saved_data_path']))
if len(load_check): # if it's there
	load_data = np.load('%s/%s/radar_times_accum.npz'%(base_path, ycfg['saved_data_path']))
	rainfall = load_data['rainfall']
	radar_times = load_data['times']
	radar_lats = load_data['radar_lats']
	radar_lons = load_data['radar_lons']
	#print 'load data: {}'.format(load_data)


# ******* Thresholds and settings
accum_threshs = [2*60.0, 4*60.0, 6*60.0, 12*60.0]


rainfall_age_thresh = 5.0 # 10 mm/hr
default_dtime = ycfg['default_dtime'] # time between scans to calculate accumulated rainfall
accum_lag_thresh = ycfg['accum_lag_thresh']
DPI = 150 # I chose this because it's the optimal size for PNGs on the monitor on the ship

image_type = 'png' # of tk looper only handles pngs
km2nm = 0.54 # 1 km is 0.54 nautical miles
nm2km = 1.0/km2nm # Just to have both conversions

# ****** Colorbar options
cb_pad = 0.05 # colorbar options
cb_frac = 0.046


major_circle_rads = np.arange(20, 100, 20)
minor_circle_rads = np.arange(10, 110, 20)
azimuths = np.arange(0, 360, 45)
minor_azimuths = np.arange(22.5, 360, 45)


# This says cs_colors here, but it's actually the rainfall
cs_colors = ['#FFFFFF', '#ccd1d1', '#99a3a4', '#707b7c', '#f4ec3f', '#f3bd4a', '#f39c12', '#ef776a', '#C85046',
			'#8B1412', '#600000']
cs_cmap = mpl.colors.ListedColormap(cs_colors)
#cs_bounds = np.arange(0,5)
#cs_bounds = np.array([0.0, 0.1, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 150.0])
cs_bounds = np.array([0, 1, 2, 3, 5, 8, 12, 18, 25, 40, 75, 150])
cs_norm = mpl.colors.BoundaryNorm(cs_bounds, cs_cmap.N)
rrage_colors = ['#505050', '#AE33DF', '#7E00BF', '#1D00C6', '#00BCD4', '#00DB81',
			'#E7F000', '#F78000', '#FF0000', '#777777', '#AAAAAA', '#CECECE', '#EDEDED']

rrage_cmap = mpl.colors.ListedColormap(rrage_colors)
rrage_bounds = np.array([0, 5, 10, 20, 30, 40, 60, 80, 100, 120, 180, 240, 360, 480])
rr_labs = ['0', '5', '10', '20', '30', '40', '60', '80', '100', '2h', '3h', '4h', '6h', '8h']
rrage_norm = mpl.colors.BoundaryNorm(rrage_bounds, rrage_cmap.N)


rraccum_colors = ['#ffffff', '#cecece', '#7effa1', '#27DA55', '#228e32', '#4Dc1FF', '#007cff', '#3333ef',
            '#E05BF3', '#A629DD', '#5E009F']

rraccum_cmap = mpl.colors.ListedColormap(rraccum_colors)
rraccum_bounds = np.array([0, 1, 2, 3, 5, 10, 20, 30, 50, 75, 150, 300])
rraccum_norm = mpl.colors.BoundaryNorm(rraccum_bounds, rraccum_cmap.N)

# ******** Some map options

map_limits = ycfg['map_limits']


# this is how wide the plot will be
lat_width = ycfg['lat_width'] # degrees
lon_width = ycfg['lon_width'] # degrees

# ********* Read in file(s) *************

map_res = 1.0 # km

dlat = map_res/111.0
dlon = map_res/(111.0*np.cos(np.radians(map_limits['lat'][0])))

nlat = int(np.round((map_limits['lat'][1]-map_limits['lat'][0])/dlat))
nlon = int(np.round((map_limits['lon'][1]-map_limits['lon'][0])/dlon))

big_lat_array = np.arange(map_limits['lat'][0], map_limits['lat'][1], dlat)[:nlat]
big_lon_array = np.arange(map_limits['lon'][0], map_limits['lon'][1], dlon)[:nlon]

tdiffs = np.array([(radar_times[-1] - _).total_seconds() for _ in radar_times])

#print 'before tdiffs: {}'.format(before_tdiffs)
if (len(tdiffs) > 0) and (tdiffs[-1] <= accum_lag_thresh*60.0) and (tdiffs[-1] > 0.):
    print 'before tdiff: {}'.format(before_tdiffs[-1])
    dtime = tdiffs[-1]/60.0
else:
    dtime = default_dtime

print 'dtime: {} minutes'.format(dtime)



# These are the x (lon) and y (lat) limits for plotting later on
xl = [radar_lons[-1]-lon_width/2.0, radar_lons[-1]+lon_width/2.0]
yl = [radar_lats[-1]-lat_width/2.0, radar_lats[-1]+lat_width/2.0]

print 'xlim: {}'.format(xl)
print 'ylim: {}'.format(yl)


# Need to add in the rain rate as a field to the RadarData object, also need to add the convective/stratiform



accums = []

for at in accum_threshs:
    good_file_ind = np.where( (tdiffs <= at*60.0) & (tdiffs >= 0) )
    accums.append(np.ma.masked_invalid(np.sum(rainfall[good_file_ind], axis=0)))


# start by getting the time differences between the current file and previous ones

#print 'after tdiffs: {}'.format(after_tdiffs)
# Figure out which files are within the accum_minute_thresh1 (30 minutes)
#good_files1 = np.where( (after_tdiffs <= accum_minute_thresh1*60.0) & (after_tdiffs >= 0) )
# Same with the 60 minutes. We could change these values if needed
#good_files2 = np.where( (after_tdiffs <= accum_minute_thresh2*60.0) & (after_tdiffs >= 0) )

# Here we sum up all the files that fall within the certain timeframes
#accum_rainfall_good_times1 = np.ma.masked_invalid(np.sum(rainfall[good_files1], axis=0))
#accum_rainfall_good_times2 = np.ma.masked_invalid(np.sum(rainfall[good_files2], axis=0))



#before_ind = np.where(after_tdiffs >= 0.)[0]
#print 'before ind: {}'.format(before_ind)

# ***** Alright, now how do I do the rainfall age thing??


# Commenting out the rainfall age for now, I think this is just going to be accumulations here
#rainfall_age = np.zeros_like(big_rain_array)
#rainfall_age = np.ma.masked_where(rainfall_age == 0.0, rainfall_age)


# set to 0's here
#accrain_array = np.array(rainfall)
#enough_rain = np.where(rainfall[before_ind]/(default_dtime/60.0) >= rainfall_age_thresh)

# now need to loop thru each entry in enough_rain
#xy_pairs = zip(enough_rain[1], enough_rain[2])

# Only go thru unique x,y pairs
#uniq_xy_pairs = list(set(xy_pairs))


#for pt in uniq_xy_pairs:
## pt is an x,y pair
## find all the times when it rained more than the threshold value at the given point
#    all_pts = np.where( (enough_rain[1] == pt[0]) & (enough_rain[2] == pt[1]) )
#    t_inds_thresh = enough_rain[0][all_pts[0]]
## grab the last time
#    last_t = t_inds_thresh[-1] # figure out the time index where it was last raining
##print t_inds_thresh
## figure out what file we're looking at here. This logic is done to allow you to run a file
## that is not the most current and it will only look backwards relative to that file
#    this_file_ind = np.where(radar_times == file_time)[0][0]
## Here we just multiply the difference in index by the dtime, which will be set to 4 minutes
#    this_age = (radar_times[this_file_ind] - radar_times[last_t]).total_seconds()/60.0
#this_age = (this_file_ind - last_t)*default_dtime
#    if this_age <= ycfg['rainfall_age_thresh']: # only assign a value if it's under the threshold
#		# if we didn't have any thresholds, the map would be covered and unintelligible
#	rainfall_age[pt] = this_age




#print 'Making rainfall figure'
#print 'in realtime figure making: {}'.format(pargs.realtime)
# Now plot the rainfall products here. Since we've made them part of the radar object, we can use the plotting framework
fig, ax = plt.subplots(2,2, figsize=(10.5, 8.5))
axf = ax.flatten()

#rr_pc = radar.cappi('RR', z=z_val, ax=axf[0], cmap=cs_cmap, norm=cs_norm, coords='ll', xlim=xl, ylim=yl)
#acc_pc1 = axf[0].pcolormesh(map_lons, map_lats, accums[0], cmap=cs_cmap, norm=cs_norm)
#rr_pc = radar.cappi('RR', z=z_val, ax=axf[0], cmap=plt.cm.GnBu, coords='ll', xlim=xl, ylim=yl, norm=LogNorm(vmin=1e-1, vmax=100))

arr_pc1 = axf[0].pcolormesh(big_lon_array, big_lat_array, accums[0], cmap=rraccum_cmap, norm=rraccum_norm)
arr_cb1 = plt.colorbar(arr_pc1, orientation='vertical', ax=axf[0], fraction=cb_frac, pad=cb_pad)

arr_cb1.set_ticks(rraccum_bounds)
arr_cb1.set_ticklabels(rraccum_bounds)

axf[0].set_title('Last %d hours'%(accum_threshs[0]/60.0))
#axf[0].set_xlabel('') # override the default here
#axf[0].set_ylabel('')


arr_pc2 = axf[1].pcolormesh(big_lon_array, big_lat_array, accums[1], cmap=rraccum_cmap, norm=rraccum_norm)

arr_cb2 = plt.colorbar(arr_pc2, orientation='vertical', ax=axf[1], fraction=cb_frac, pad=cb_pad)
arr_cb2.set_ticks(rraccum_bounds)
arr_cb2.set_ticklabels(rraccum_bounds)

axf[1].set_title('Last %d hours'%(accum_threshs[1]/60.0))


arr_pc3 = axf[2].pcolormesh(big_lon_array, big_lat_array, accums[2], cmap=rraccum_cmap, norm=rraccum_norm)

arr_cb3 = plt.colorbar(arr_pc3, orientation='vertical', ax=axf[2], fraction=cb_frac, pad=cb_pad)
arr_cb3.set_ticks(rraccum_bounds)
arr_cb3.set_ticklabels(rraccum_bounds)

axf[2].set_title('Last %d hours'%(accum_threshs[2]/60.0))


arr_pc4 = axf[3].pcolormesh(big_lon_array, big_lat_array, accums[3], cmap=rraccum_cmap, norm=rraccum_norm)

arr_cb4 = plt.colorbar(arr_pc4, orientation='vertical', ax=axf[3], fraction=cb_frac, pad=cb_pad)
arr_cb4.set_ticks(rraccum_bounds)
arr_cb4.set_ticklabels(rraccum_bounds)

axf[3].set_title('Last %d hours'%(accum_threshs[3]/60.0))


#radar_lats = radar_lats+np.linspace(0, 0.4, 4)

if len(radar_lats) > 1:

    dlat = np.diff(radar_lats)[-1]
    dlon = np.diff(radar_lons)[-1]
    dy = dlat*111.0
    dx = dlon*(111.0*np.cos(np.radians(radar_lats[-1])))
    speed_km = np.sqrt(dx**2 + dy**2)*60.0/dtime
    speed = km2nm*speed_km

# figure out a heading and throw it on there

    heading_physics = np.degrees(np.arctan2(dy, dx))
    heading = 90. - heading_physics
    if heading < 0.:
    	heading += 360.0

else:
    speed = 0.
    heading = 0.


#print 'speed: {}, heading: {}'.format(speed, heading)


for ia, a in enumerate(axf):

    a.set_xlim(*xl)
    a.set_ylim(*yl)


    a.grid(True)

    a.scatter(radar_lons, radar_lats, c=tdiffs, cmap=plt.cm.bone, s=15, vmin=-100, vmax=3600.0*6, edgecolors='none', alpha=0.7)

#             if False: # turning off the arrow plotting now, not really sure it's useful
# #            if speed >= 2.0:
#             	print 'plotting an arrow'
#             	a.arrow(radar_lat, radar_lon, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')

    axt = a.get_xticks().tolist()
    for j in range(len(axt)):
	axt[j] = convert_decimal_to_degree(axt[j])
    a.set_xticklabels(axt)

    axty = a.get_yticks().tolist()
#print 'y tick labels: {}'.format(axty)
    for j in range(len(axty)):
	axty[j] = convert_decimal_to_degree(axty[j])
    a.set_yticklabels(axty)



    for label in a.get_xticklabels():
	label.set_rotation(25) 


    for label in a.get_yticklabels():
	label.set_rotation(25) 

    a.tick_params(labelsize=8)

    plot_range_rings(a, radar_lats[-1], radar_lons[-1], major_circle_rads, minor_rads=minor_circle_rads)
    plot_azimuths(a, radar_lats[-1], radar_lons[-1], azimuths, minor_azs=minor_azimuths)
    a.set_aspect('equal')


#plt.xticks(rotation=45)
file_dt_string = radar_times[-1].strftime('%Y-%m-%d %H:%M:%S')
save_dt_string = radar_times[-1].strftime('%Y%m%d_%H%M%S')

fig.suptitle('SEAPOL: %s rainfall accumulation'%(file_dt_string), fontsize=16)
plt.tight_layout()
fig.subplots_adjust(top=0.92)



plt.savefig('%s/figures/longterm_rain/%s_longtermaccum.%s'%(base_path, save_dt_string, image_type), dpi=DPI)
plt.savefig('%s/figures/animated_gifs/longtermaccum_latest.%s'%(base_path, image_type), dpi=DPI)



plt.close(fig)



print 'Made it to the end of long term processing'







