# Brody Fuchs, Sept 2017
# brfuchs@atmos.colostate.edu

# March 2018: This is all of the rainfall plotting (and some analysis) code 
# that has been separated from the gridded processing

from __future__ import division

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
import cPickle as pickle
import Config
import gentools as gtools
import plot_tools as ptools
import radar_tools as rtools
import os
import scipy


try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception:
    base_path = os.path.dirname(os.path.realpath('__file__'))


cb_pad = 0.05 # colorbar options
cb_frac = 0.046
dt_fmt = '%Y%m%d_%H%M%S_rain.npz'

image_type = 'png'



parser = argparse.ArgumentParser(description='Put in a file to be processed')

parser.add_argument('--file', action="store", dest="file", default=None)
parser.add_argument('--config', action="store", dest="config", default=None)
parser.add_argument('--plot', action="store", dest="plot", default=1, type=int)
parser.add_argument('--save', action="store", dest="save", default=0, type=int)
parser.add_argument('--accumtimes', action="store", dest="accumtimes", type=int, nargs='+', default=[30, 60])



pargs = parser.parse_args()

accum_times = np.array(pargs.accumtimes)

if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v



if 'DPI' in ycfg.keys():
    DPI = ycfg['DPI']
else:
    DPI = 120

if 'largeDPI' in ycfg.keys():
    largeDPI = ycfg['largeDPI']
else:
    largeDPI = 200

if 'radar_title' in ycfg.keys():
    radar_title = ycfg['radar_title']
else:
    radar_title = 'SEAPOL'



major_circle_rads = np.arange(20, 100, 20)
minor_circle_rads = np.arange(10, 110, 20)

small_circle_rads = np.arange(10, 40, 10)

azimuths = np.arange(0, 360, 30)


# this will take in a filename (npz) to read in and then look backwards to do accumulation

fname = pargs.file
fbase = os.path.basename(fname)


# Try to read in the file
if fname is not None:
    pass
    raindata = np.load(fname)



# This says cs_colors here, but it's actually the rainfall
rr_colors = ['#FFFFFF', '#ccd1d1', '#99a3a4', '#707b7c', '#f4ec3f', '#f3bd4a', '#f39c12', '#ef776a', '#C85046',
            '#9B2422', '#600000']
rr_cmap = mpl.colors.ListedColormap(rr_colors)
rr_bounds = np.array([0, 2, 3, 5, 10, 20, 30, 50, 75, 150, 225, 300])
rr_norm = mpl.colors.BoundaryNorm(rr_bounds, rr_cmap.N)


rrage_colors = ['#505050', '#AE33DF', '#7E00BF', '#1D00C6', '#00BCD4', '#00DB81',
            '#E7F000', '#F78000', '#FF0000', '#777777', '#AAAAAA', '#CECECE', '#EDEDED']

rrage_cmap = mpl.colors.ListedColormap(rrage_colors)
rrage_bounds = np.array([0, 5, 10, 20, 30, 40, 60, 80, 100, 120, 180, 240, 360, 480])
rr_labs = ['0', '5', '10', '20', '30', '40', '60', '80', '100', '2h', '3h', '4h', '6h', '8h']
rrage_norm = mpl.colors.BoundaryNorm(rrage_bounds, rrage_cmap.N)


#rraccum_colors = ['#ffffff', '#dedede', '#aedeae', '#7effa1', '#27DA55', '#4Dc1FF', '#007cff', '#3300cf',
            #'#AE33DF', '#7E00BF', '#4C0083']

rraccum_colors = ['#ffffff', '#cecece', '#7effa1', '#27DA55', '#228e32', '#4Dc1FF', '#007cff', '#3333ef',
            '#E05BF3', '#A629DD', '#5E009F']

rraccum_cmap = mpl.colors.ListedColormap(rraccum_colors)
rraccum_bounds = np.array([0, 1, 2, 3, 5, 10, 20, 30, 50, 75, 150, 300])
rraccum_norm = mpl.colors.BoundaryNorm(rraccum_bounds, rraccum_cmap.N)


# accum_minute_thresh1 = 30. # accumulated rainfall in last 30 minutes
# accum_minute_thresh2 = 60. # accumulated rainfall in last 60 minutes

# # These thresholds are for a longer term plot that is separate from the realtime plot
# long_amt = 12*60.

accums = [] # this will hold the accumulation arrays, which need to be flexible

rainfall_age_thresh = 10.0 # 10 mm/hr
default_dtime = ycfg['default_dtime'] # time between scans to calculate accumulated rainfall
accum_lag_thresh = ycfg['accum_lag_thresh']

fig_path = ycfg['fig_path']


prefix = ycfg['qc_prefix']

mid_pt_lat = int(np.floor(raindata['lat'].shape[0]/2.0))
mid_pt_lon = int(np.floor(raindata['lon'].shape[0]/2.0))

radar_lat = np.average(raindata['lat'][mid_pt_lat,0])
radar_lon = np.average(raindata['lon'][0,mid_pt_lon])

# this is how wide the plot will be
lat_width = ycfg['lat_width'] # degrees
lon_width = ycfg['lon_width'] # degrees

xl = [radar_lon-lon_width/2.0, radar_lon+lon_width/2.0]
yl = [radar_lat-lat_width/2.0, radar_lat+lat_width/2.0]


first_lat = raindata['lat'][0,0]
first_lon = raindata['lon'][0,0]

dlat = np.average(np.diff(raindata['lat'][:, 0]))
dlon = np.average(np.diff(raindata['lon'][0]))


# Now go in and look for files that might be within 60 minutes of the specified file

dt_string = fbase[0:15]

file_time = datetime.datetime.strptime(fbase, dt_fmt)
file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')
date_string = file_time.strftime('%Y%m%d')

other_files = sorted(glob.glob('%s/*.npz'%(os.path.dirname(fname))))

file_time_strings = np.array( [gtools.parse_time_string(os.path.basename(_)) for _ in other_files if _ != fname] )
file_times = np.array( [datetime.datetime.strptime(_, dt_fmt) for _ in file_time_strings] )
ft_diffs = np.array([(_ - file_time).total_seconds() for _ in file_times])
#
#    # now need to figure out which files are between the beginning and the end
#

max_time_diff = np.max( [3*60.0, np.max(accum_times)] )


#valid_file_ind = np.where( (ft_diffs >= -60.0*long_amt) & (ft_diffs < 0) )[0]
valid_file_ind = np.where( (ft_diffs >= -1.0*60.0*max_time_diff) & (ft_diffs < 0) )[0]



if len(valid_file_ind): # means that there a files in the near past to add and do age and things like that
    #print valid_file_ind

    # basically need to loop thru valid files, read them in and then concatentate them
    rainfall_arr = []
    radar_lats = []
    radar_lons = []
    tdiffs = []


    for vi in valid_file_ind:
        this_array = np.load(other_files[vi])
        this_rainfall = this_array['rainfall']
        #if rainfall_arr.ndim == 2:
        #    rainfall_arr = np.concatenate((this_rainfall[np.newaxis, ...], rainfall_arr[np.newaxis, ...]), axis=0)
        #else:
        #    rainfall_arr = np.concatenate((this_rainfall[np.newaxis, ...], rainfall_arr), axis=0)
        rainfall_arr.append(this_rainfall)

        radar_lats.append(np.average(this_array['lat'][mid_pt_lat,0]))
        radar_lons.append(np.average(this_array['lon'][0,mid_pt_lon]))
        tdiffs.append(ft_diffs[vi])

        # the -1 index is the last one, which is the latest time (and the reference time)
    # after that, tack on the current one
    rainfall_arr.append(raindata['rainfall'])
    radar_lats.append(radar_lat)
    radar_lons.append(radar_lon)
    tdiffs.append(0.)

    rainfall_arr = np.array(rainfall_arr)
    radar_lats = np.array(radar_lats)
    radar_lons = np.array(radar_lons)
    tdiffs = np.array(tdiffs)

    lat_shifts = radar_lats - radar_lats[-1]
    lon_shifts = radar_lons - radar_lons[-1]

    lat_shifts_int = np.round(lat_shifts/dlat).astype(int)
    lon_shifts_int = np.round(lon_shifts/dlon).astype(int)

    # now do the shifting to be relative to the last radar frame
    for i_time in range(len(lat_shifts)):
        rainfall_arr[i_time] = np.roll(rainfall_arr[i_time], (lat_shifts_int[i_time], lon_shifts_int[i_time]), axis=(0,1))
        # only problem with np.roll is the values wrap around, which is not what we want, so we'll need to do
        # make those values 0.
    


    # I think what I'm gonna do here is make a list of accumulations based on the accum_times
    # Then just plot the first two?
    for _at in accum_times:
        accums.append(rtools.rainfall_accumulation(rainfall_arr, tdiffs, _at))


    # # Maybe let's just get the aggregation going, then we can do the shifts afterword
    # accum_60 = rtools.rainfall_accumulation(rainfall_arr, tdiffs, 60)
    # accum_30 = rtools.rainfall_accumulation(rainfall_arr, tdiffs, 30)
    rrage = rtools.rainfall_age(rainfall_arr, tdiffs, rainfall_thresh=rainfall_age_thresh, default_dtime=20.0)


    if pargs.plot:

    # ************* Now let's make the figure ****************************************

        fig, ax = plt.subplots(2,2, figsize=(10.5, 8.5))
        axf = ax.flatten()

        rr_pc = axf[0].pcolormesh(raindata['lon'], raindata['lat'], raindata['rainfall'], cmap=rr_cmap, norm=rr_norm)
        axf[0].set_title('Instantaneous rain rate (mm/hr)')

        rr_cb = plt.colorbar(rr_pc, orientation='vertical', ax=axf[0], fraction=cb_frac, pad=cb_pad)
        rr_cb.set_ticks(rr_bounds)
        rr_cb.set_ticklabels(rr_bounds)

   
    # # if True:
    # #     try:
    # #            # get the heading, subtract 180 to flip it, then get that + 60 and - 60
    # #         opp_hd = heading - 180.0
    # #         blanks = 90.0 - np.array([opp_hd-60.0, opp_hd+60.0])
    # #         blanks[blanks<0] += 360.0
    # #        print 'blanks: {}'.format(blanks)

    # #        plot_azimuths(axf[0], radar_lat, radar_lon, blanks, color='red', linewidth=1.0, alpha=0.4, linestyle='solid')

    # #    except Exception, bse:
    # #        print 'Error with plotting the blanked sector: {}'.format(bse)


        rrage_pc = axf[1].pcolormesh(raindata['lon'], raindata['lat'], rrage, 
                cmap=rrage_cmap, norm=rrage_norm)
        rrage_cb3 = plt.colorbar(rrage_pc, orientation='vertical', ax=axf[1], fraction=cb_frac, pad=cb_pad)
        rrage_cb3.set_ticks(rrage_bounds)
        rrage_cb3.set_ticklabels(rr_labs)
        axf[1].set_title('Rainfall age (min) using %d mm/hr threshold'%(rainfall_age_thresh))



        # Now down here, wanna make this more flexible. Check the length of the accums array. If it's less than 2
        # then just plot one in axf[2]. If 2 or more, just plot the first 2.

        if len(accums):

            arr_pc1 = axf[2].pcolormesh(raindata['lon'], raindata['lat'], accums[0], cmap=rraccum_cmap, 
                       norm=rraccum_norm)

            arr_cb1 = plt.colorbar(arr_pc1, orientation='vertical', ax=axf[2], fraction=cb_frac, pad=cb_pad)
            arr_cb1.set_ticks(rraccum_bounds)
            arr_cb1.set_ticklabels(rraccum_bounds)
            axf[2].set_title('%d min accum rainfall (mm)'%(accum_times[0]))

   
            if len(accums) >= 2:

                arr_pc2 = axf[3].pcolormesh(raindata['lon'], raindata['lat'], accums[1], 
                        cmap=rraccum_cmap, norm=rraccum_norm)
                arr_cb2 = plt.colorbar(arr_pc2, orientation='vertical', ax=axf[3], fraction=cb_frac, pad=cb_pad)

                arr_cb2.set_ticks(rraccum_bounds)
                arr_cb2.set_ticklabels(rraccum_bounds)
                axf[3].set_title('%d min accum rainfall (mm)'%(accum_times[1]))



        for ia, a in enumerate(axf):
            a.grid(True)
            a.set_xlim(*xl)
            a.set_ylim(*yl)

            a.scatter(radar_lons, radar_lats, c=tdiffs, cmap=plt.cm.bone_r, 
                       s=15, vmin=-1*3600.0*6, vmax=0, edgecolors='none', alpha=0.5)

#         if False: # turning off the arrow plotting now, not really sure it's useful
# #            if speed >= 2.0:
#           print 'plotting an arrow'
#           a.arrow(radar_lat, radar_lon, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')



            axt = a.get_xticks().tolist()
            for j in range(len(axt)):
                axt[j] = gtools.convert_decimal_to_degree(axt[j])
            a.set_xticklabels(axt)

            axty = a.get_yticks().tolist()
            for j in range(len(axty)):
                axty[j] = gtools.convert_decimal_to_degree(axty[j])
            a.set_yticklabels(axty)



            for label in a.get_xticklabels():
                label.set_rotation(25) 


            for label in a.get_yticklabels():
                label.set_rotation(25) 

            a.tick_params(labelsize=8)

            ptools.plot_range_rings(a, radar_lat, radar_lon, major_circle_rads, minor_rads=minor_circle_rads, text_flag=False)
            ptools.plot_azimuths(a, radar_lat, radar_lon, azimuths)
            a.set_aspect('equal')


        fig.suptitle('%s: %s rainfall'%(radar_title, file_dt_string), fontsize=16)
        plt.tight_layout()
        fig.subplots_adjust(top=0.92)

    #plt.savefig('%s/figures/rain/%s_rraccum.%s'%(base_path, dt_string, image_type), dpi=DPI)
        ptools.save_figure(fig, '%s/rain/%s'%(fig_path, date_string), '%s%s_rraccum.%s'%(prefix, dt_string, image_type), dpi=DPI)

    #plt.savefig('%s/figures/animated_gifs/rraccum_latest.%s'%(base_path, image_type), dpi=DPI)

        plt.close(fig)


# Okay, now need to check on a save flag, if so then save all of that!


    if pargs.save:
        
        out = {'time': [file_dt_string], 'accums': accums, 'accum_times': accum_times, 'lat': raindata['lat'], 
                                'lon': raindata['lon'], 'age': rrage, 'rainrate': raindata['rainfall']}
        accums_str = '_'.join(accum_times.astype(str))

        scipy.io.savemat('%s/mat_rainfiles/%s_rainfall_products_%s.mat'%(base_path, dt_string, accums_str), out)




        pass


else: # This means there is only one file around
    print 'Only one file'
    if True:

        print 'plotting only one file'
        fig, ax = plt.subplots(1,1, figsize=(8,6))
        rr_pc = ax.pcolormesh(raindata['lon'], raindata['lat'], raindata['rainfall'], cmap=rr_cmap, norm=rr_norm)
        ax.scatter(radar_lon, radar_lat, c='black', s=30)
#rr_pc = radar.cappi('RR', z=z_val, ax=axf[0], cmap=plt.cm.GnBu, coords='ll', xlim=xl, ylim=yl, norm=LogNorm(vmin=1e-1, vmax=100))
#    ax.set_title('Instantaneous rain rate (mm/hr)')
        ax.set_xlabel('') # override the default here
        ax.set_ylabel('')
        rr_cb = plt.colorbar(rr_pc, orientation='vertical', ax=ax, fraction=cb_frac, pad=cb_pad)
        rr_cb.set_ticks(rr_bounds)
        rr_cb.set_ticklabels(rr_bounds)

        ax.grid(True)
        ptools.ax_labels_minutes(ax)

        ptools.rotate_ticklabs(ax, x=True, y=True)

#    ptools.plot_range_rings(ax, radar_lat, radar_lon, major_circle_rads, minor_rads=minor_circle_rads)
#    ptools.plot_azimuths(ax, radar_lat, radar_lon, azimuths)

#    ax.set_aspect('equal')
        fig.suptitle('%s: %s rain rates (mm/hr)\n(no prior files available)'%(radar_title, file_dt_string), fontsize=14)


        plt.tight_layout()
        fig.subplots_adjust(top=0.90)

        ptools.save_figure(fig, '%s/rain/%s'%(fig_path, date_string), '%s%s_rr_solo.%s'%(prefix, dt_string, image_type), dpi=DPI)




#big_rain_array = np.zeros((nlat, nlon), np.float)
#big_lat_array = np.arange(map_limits['lat'][0], map_limits['lat'][1], dlat)
#big_lon_array = np.arange(map_limits['lon'][0], map_limits['lon'][1], dlon)



#if pargs.realtime:
#    #print 'Processing radar files in realtime'
#    # Check for the accumulated radar times pickled list
#    load_check = glob.glob('%s/%s/radar_times_accum.npz'%(base_path, ycfg['saved_data_path']))
#    if len(load_check): # if it's there
#        load_data = np.load('%s/%s/radar_times_accum.npz'%(base_path, ycfg['saved_data_path']))
#        rainfall = load_data['rainfall']
#        radar_times = load_data['times']
#        radar_lats = load_data['radar_lats']
#        radar_lons = load_data['radar_lons']
#        #print 'load data: {}'.format(load_data)
#    #   print 'Found an existing radar times pickle'
#    else:
#        radar_times = np.zeros(0, object)
#        rainfall = []
#        radar_lats = np.zeros(0, float)
#        radar_lons = np.zeros(0, float)
#
#        print 'Had to create a radar times pickle'
#
#
#else:
#    radar_times = np.zeros(0)
#    rainfall = []



#    xcut_n = i_rad_lon - int(nx/2)
#    xcut_x = i_rad_lon + int(nx/2) + 1

#    ycut_n = i_rad_lat - int(ny/2)
#    ycut_x = i_rad_lat + int(ny/2) + 1

#    print 'xcut: %d - %d, ycut: %d - %d'%(xcut_n, xcut_x, ycut_n, ycut_x)
#
#    big_rain_array[ycut_n: ycut_x, xcut_n: xcut_x] += trop_rr.filled(0.0)*dtime/60.0
#
#    print 'max RR: {}'.format(trop_rr.max())
#
#
#    #if not file_in_times and file_after:
#    #print 'add data flag: {}'.format(add_data_flag)
#    if add_data_flag:
#
#        # adding radar lat and lon
#        radar_lats = np.append(radar_lats, radar_lat)
#        radar_lons = np.append(radar_lons, radar_lon)
#
#        print 'This time is not already included in the saved data, so we are adding it'
#
#        if len(rainfall) == 0:
#        # if nothing in rainfall
#            print 'making a new rainfall array'
#            rainfall = deepcopy(big_rain_array[np.newaxis, ...])
#
#        else:
#        # if rainfall already has some values in it 
#            print 'concatenating the rainfall array'
#            rainfall = np.concatenate([rainfall, big_rain_array[np.newaxis, ...]])
#
#    else:
#        print 'This is not a new time, not adding it to the saved data'
#
#
#
#    #print 'rainfall array shape: {}'.format(rainfall.shape)
#
#    after_tdiffs = np.array([(file_time - _).total_seconds() for _ in radar_times])
#    #print 'after tdiffs max: {}'.format(after_tdiffs.max())
#    #print 'last after tdiffs: {}'.format(after_tdiffs[-1])
#
#    #print 'last time diff: {}'.format(after_tdiffs[5:])
#
#    if pargs.realtime:
#
#    # start by getting the time differences between the current file and previous ones
#
#        print 'len tdiffs: {}, len rainfall: {}'.format(len(after_tdiffs), len(rainfall))
#
#        if True:
#            if len(after_tdiffs) != len(rainfall):
#                len_diff = len(after_tdiffs) - len(rainfall)
#                after_tdiffs = after_tdiffs[len_diff:]
#
#    #print 'after tdiffs: {}'.format(after_tdiffs)
#    # Figure out which files are within the accum_minute_thresh1 (30 minutes)
#        good_files1 = np.where( (after_tdiffs <= accum_minute_thresh1*60.0) & (after_tdiffs >= 0) )
#    # Same with the 60 minutes. We could change these values if needed
#        good_files2 = np.where( (after_tdiffs <= accum_minute_thresh2*60.0) & (after_tdiffs >= 0) )
#
#        good_files_long = np.where( (after_tdiffs <= long_amt2*60.0) & (after_tdiffs >= 0) )
#
#        #long_amt2
#
#    # Here we sum up all the files that fall within the certain timeframes
#        accum_rainfall_good_times1 = np.ma.masked_invalid(np.sum(rainfall[good_files1], axis=0))
#        accum_rainfall_good_times2 = np.ma.masked_invalid(np.sum(rainfall[good_files2], axis=0))
#        accum_rainfall_good_times_long = np.ma.masked_invalid(np.sum(rainfall[good_files_long], axis=0))
#
#        #radar.add_field(accum_rainfall_good_times1, 'ACCRR30')
#        #radar.add_field(accum_rainfall_good_times2, 'ACCRR60')
#
#
#        before_ind = np.where(after_tdiffs >= 0.)[0]
#    #print 'before ind: {}'.format(before_ind)
#
#    # ***** Alright, now how do I do the rainfall age thing??
#
#
#    # Start by making an array with all masked values
#        rainfall_age = np.zeros_like(big_rain_array)
#        rainfall_age = np.ma.masked_where(rainfall_age == 0.0, rainfall_age)
#
#
#    # set to 0's here
#    #accrain_array = np.array(rainfall)
#        enough_rain = np.where(rainfall[before_ind]/(default_dtime/60.0) >= rainfall_age_thresh)
#
#    # now need to loop thru each entry in enough_rain
#        xy_pairs = zip(enough_rain[1], enough_rain[2])
#
#    # Only go thru unique x,y pairs
#        uniq_xy_pairs = list(set(xy_pairs))
#
#
#        for pt in uniq_xy_pairs:
#        # pt is an x,y pair
#        # find all the times when it rained more than the threshold value at the given point
#            all_pts = np.where( (enough_rain[1] == pt[0]) & (enough_rain[2] == pt[1]) )
#            t_inds_thresh = enough_rain[0][all_pts[0]]
#        # grab the last time
#            last_t = t_inds_thresh[-1] # figure out the time index where it was last raining
#        #print t_inds_thresh
#        # figure out what file we're looking at here. This logic is done to allow you to run a file
#        # that is not the most current and it will only look backwards relative to that file
#            this_file_ind = np.where(radar_times == file_time)[0][0]
#        # Here we just multiply the difference in index by the dtime, which will be set to 4 minutes
#            this_age = (radar_times[this_file_ind] - radar_times[last_t]).total_seconds()/60.0
#            #this_age = (this_file_ind - last_t)*default_dtime
#            if this_age <= ycfg['rainfall_age_thresh']: # only assign a value if it's under the threshold
#                    # if we didn't have any thresholds, the map would be covered and unintelligible
#                rainfall_age[pt] = this_age
#
#        #radar.add_field(rainfall_age, 'RRAGE')
#
#
#        print 'Making rainfall figure'
#        #print 'in realtime figure making: {}'.format(pargs.realtime)
#    # Now plot the rainfall products here. Since we've made them part of the radar object, we can use the plotting framework




##############################################################################################
##############################################################################################
        # making the same plot, but zooming in!
#        fig, ax = plt.subplots(2,2, figsize=(10.5, 8.5))
#        axf = ax.flatten()
#
#        #rr_pc = radar.cappi('RR', z=z_val, ax=axf[0], cmap=cs_cmap, norm=cs_norm, coords='ll', xlim=xl, ylim=yl)
#        rr_pc = axf[0].pcolormesh(radar_lon_array, radar_lat_array, trop_rr, cmap=cs_cmap, norm=cs_norm)
#    #rr_pc = radar.cappi('RR', z=z_val, ax=axf[0], cmap=plt.cm.GnBu, coords='ll', xlim=xl, ylim=yl, norm=LogNorm(vmin=1e-1, vmax=100))
#        axf[0].set_title('Instantaneous rain rate (mm/hr)')
#        axf[0].set_xlabel('') # override the default here
#        axf[0].set_ylabel('')
#        rr_cb = plt.colorbar(rr_pc, orientation='vertical', ax=axf[0], fraction=cb_frac, pad=cb_pad)
#        rr_cb.set_ticks(rr_bounds)
#        rr_cb.set_ticklabels(rr_bounds)
#
#        
#        rrage_pc = axf[1].pcolormesh(big_lon_array, big_lat_array, rainfall_age, 
#                     cmap=rrage_cmap, norm=rrage_norm)
#        rrage_cb3 = plt.colorbar(rrage_pc, orientation='vertical', ax=axf[1], fraction=cb_frac, pad=cb_pad)
#        rrage_cb3.set_ticks(rrage_bounds)
#        rrage_cb3.set_ticklabels(rr_labs)
#
#
#        arr_pc1 = axf[2].pcolormesh(big_lon_array, big_lat_array, accum_rainfall_good_times1, cmap=rraccum_cmap, 
#                                        norm=rraccum_norm)
#
#        arr_cb1 = plt.colorbar(arr_pc1, orientation='vertical', ax=axf[2], fraction=cb_frac, pad=cb_pad)
#        arr_cb1.set_ticks(rraccum_bounds)
#        arr_cb1.set_ticklabels(rraccum_bounds)
#
#        axf[1].set_title('Rainfall age (min) using %d mm/hr threshold'%(rainfall_age_thresh))
#
#        axf[2].set_title('30 min accum rainfall (mm)')
#        
#
#        arr_pc2 = axf[3].pcolormesh(big_lon_array, big_lat_array, accum_rainfall_good_times2, 
#                     cmap=rraccum_cmap, norm=rraccum_norm)
#        arr_cb2 = plt.colorbar(arr_pc2, orientation='vertical', ax=axf[3], fraction=cb_frac, pad=cb_pad)
#
#        arr_cb2.set_ticks(rraccum_bounds)
#        arr_cb2.set_ticklabels(rraccum_bounds)
#    
#        
#        axf[3].set_title('60 min accum rainfall (mm)')
#
#        for ia, a in enumerate(axf):
#            a.grid(True)
#
#
#            a.set_xlim(*xl_zoom)
#            a.set_ylim(*yl_zoom)
#
#            a.scatter(radar_lons[before_ind], radar_lats[before_ind], c=after_tdiffs[before_ind], cmap=plt.cm.bone, 
#                                        s=15, vmin=-100, vmax=3600.0*6, edgecolors='none', alpha=0.5)
#

#            axt = a.get_xticks().tolist()
#            for j in range(len(axt)):
#                axt[j] = convert_decimal_to_degree(axt[j])
#            a.set_xticklabels(axt)

#            axty = a.get_yticks().tolist()
#            #print 'y tick labels: {}'.format(axty)
#            for j in range(len(axty)):
#                axty[j] = convert_decimal_to_degree(axty[j])
#            a.set_yticklabels(axty)

#            for label in a.get_xticklabels():
#                label.set_rotation(25) 


#            for label in a.get_yticklabels():
#                label.set_rotation(25) 

#            a.tick_params(labelsize=8)

#            plot_range_rings(a, radar_lat, radar_lon, small_circle_rads)
#            plot_azimuths(a, radar_lat, radar_lon, azimuths)
#            a.set_aspect('equal')

#        fig.suptitle('SEAPOL: %s zoomed rainfall'%(file_dt_string), fontsize=16)
#        plt.tight_layout()
#        fig.subplots_adjust(top=0.92)

##        plt.savefig('%s/figures/rain/%s_rraccum.%s'%(base_path, dt_string, image_type), dpi=DPI)
#        plt.savefig('%s/figures/animated_gifs/rraccum_latest_zoom.%s'%(base_path, image_type), dpi=DPI)
#
#        plt.close(fig)
