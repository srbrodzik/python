# Brody Fuchs, Sept 2017
# brfuchs@atmos.colostate.edu

# Just some code to get started reading in radar data and doing some plotting and rainrate stuff
# This version will use the csuram stuff to make it easier

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
import gc
import scipy.ndimage as ndi
import gentools as gtools 
import plot_tools as ptools 
import radar_tools as radtools 


start = time.time()

base_path = os.path.dirname(os.path.realpath('__file__'))
radar_path = '%s/realtime_gridded'%(base_path)
file_type = 'nc'

parser = argparse.ArgumentParser(description='Put in a file to be processed')

#parser.add_argument('--noarg', action="store_true", default=False)
parser.add_argument('--file', action="store", dest="file")
parser.add_argument('--lastfile', action="store", dest="lastfile", type=bool, default=False)
parser.add_argument('--allfiles', action="store", dest="allfiles", type=bool, default=False)
parser.add_argument('--realtime', action="store", dest="realtime", type=bool, default=False)
parser.add_argument('--config', action="store", dest="config", default=None)


pargs = parser.parse_args()

#print pargs


if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v




# This is now separate and pertains to the file(s) being read in

if pargs.file is not None:
    radar_files = [pargs.file]

elif pargs.allfiles:
    radar_files = sorted(glob.glob('%s/*.%s'%(radar_path, file_type)))

elif pargs.lastfile:
    radar_files = [glob.glob('%s/*.%s'%(radar_path, file_type))[-1]]

else:
    print 'No files specified or found, exiting.'
    sys.exit()




# ******* Thresholds and settings


DPI = 150 # I chose this because it's the optimal size for PNGs on the monitor on the ship

image_type = 'png' # of tk looper only handles pngs
km2nm = 0.54 # 1 km is 0.54 nautical miles
nm2km = 1.0/km2nm # Just to have both conversions

# ****** Colorbar options
cb_pad = 0.05 # colorbar options
cb_frac = 0.046


major_circle_rads = np.arange(20, 100, 20)
minor_circle_rads = np.arange(10, 110, 20)

small_circle_rads = np.arange(10, 40, 10)

azimuths = np.arange(0, 360, 30)
#minor_azimuths = np.arange(22.5, 360, 45)
minor_azimuths = None



# ******** Some map options

#map_limits = {'lat': [38.0, 42.5], 'lon': [-107.0, -102.0]}
#map_limits = {'lat': [23.0, 30.5], 'lon': [-123.0, -112.0]}
map_limits = ycfg['map_limits']
default_dtime = ycfg['default_dtime']


# this is how wide the plot will be
lat_width = ycfg['lat_width'] # degrees
lon_width = ycfg['lon_width'] # degrees


lat_width_zoom = ycfg['lat_width_zoom']
lon_width_zoom = ycfg['lon_width_zoom']


# Here we're just going to use a sounding I ripped off from Hawaii on UWYO's website
# We need the sounding for the HID algorithm as well as the rainfall algorithm (which uses HID, among other things)
# We just need to know about where the freezing level is, which won't change more than a couple hundred meters
sounding_file = '%s/sample_tropical_sounding.txt'%(base_path)
# ******** Need to read this file and get it into the HID calculation **********

# using the 
snd = SkewT.Sounding(sounding_file)



# ***** defining the C/S algorithm parameters, these values were lifted from some of Brenda's code
# We could look into changing some of these if we feel it's necessary

minZdiff = 20 
deepcoszero = 40
shallowconvmin = 28
truncZconvthres = 43;
dBZformaxconvradius = 46;
weakechothres = 7
backgrndradius = 5       #(in km)
maxConvRadius = 10       #(in km)
minsize = 8              #(in km^2)
startslope = 50          #(in km^2)
maxsize = 2000           #(in km^2)


# ********* Read in file(s) ***************

add_data_flag = False

for irf, rf in enumerate(radar_files):
    # read in the file
#    print 'starting: {}'.format(time.time() - start)


    # The radar is moving but the domain over which we're accumulating rainfall info/stats will not
    # This means we need to construct the same lat/lon array no matter where the radar is.
    # We'll take the radar location info and use it to populate a portion of this larger array
    # later in the code



    # Here we're using the csuram module that I developed. You give it the radar file
    # as well as the names of the standard variables in the netcdf. 

    try:
    	radar = RadarData.RadarData(radar_file=rf, x='x0', y='y0', z='z0', band='C', squeeze=True, 
                    lat='lat0', lon='lon0', dz=ycfg['dz_name'], kdp=ycfg['kdp_name'], zdr=ycfg['zdr_name'], rho=ycfg['rho_name'])

    	    # Here we add the sounding object to the radar object
    	radar.add_sounding_object(snd)
    	# Interpolate the sounding heights to the radar heights
    	radar.interp_sounding()
    	# Calculate the HID and add that info to the radar object
    	#radar.set_hid(band=ycfg['radar_band'], use_temp=True) # this calls the appropriate functions and uses the sounding that we've
        this_hid = radar.get_hid(band=ycfg['radar_band'], use_temp=True)
        radar.add_field(this_hid, 'HID')
        # attributed to the radar object and uses the proper radar band

    except KeyError:
    	print 'Probably started out as a binary file with different names'
    	radar = RadarData.RadarData(radar_file=rf, x='x0', y='y0', z='z0', band='C', squeeze=True, 
                    lat='lat0', lon='lon0', dz='DZQC', zdr='ZDR2', kdp='KDP', rho='RHOHV2', hid='HID')    	

    	radar.add_sounding_object(snd)
    	# Interpolate the sounding heights to the radar heights
    	radar.interp_sounding()
    	# Calculate the HID and add that info to the radar object
    	radar.set_hid(band='C', use_temp=True) # this calls the appropriate functions and uses the sounding that we've
        # attributed to the radar object and uses the proper radar band

    map_res = radar.dx # km

    dlat = map_res/111.0
    dlon = map_res/(111.0*np.cos(np.radians(map_limits['lat'][0])))

    nlat = int(np.round((map_limits['lat'][1]-map_limits['lat'][0])/dlat))
    nlon = int(np.round((map_limits['lon'][1]-map_limits['lon'][0])/dlon))


    big_rain_array = np.zeros((nlat, nlon), np.float)
    big_lat_array = np.arange(map_limits['lat'][0], map_limits['lat'][1], dlat)
    big_lon_array = np.arange(map_limits['lon'][0], map_limits['lon'][1], dlon)


    # Just grabbing the time of the radar file from the name and appending it to the accum_times list
    last_dot = rf.rfind('.')
    d_index = rf.rfind('d')

    dt_string = rf[d_index+2:d_index+17]
    print dt_string


    file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
    radar.date = file_time

    file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')
    radar_times = np.zeros(0, object)

    before_tdiffs = np.array([(file_time - _).total_seconds() for _ in radar_times])

    #print 'before tdiffs: {}'.format(before_tdiffs)
    if (len(before_tdiffs) > 0) and (before_tdiffs[-1] <= accum_lag_thresh*60.0) and (before_tdiffs[-1] > 0.):
    	print 'before tdiff: {}'.format(before_tdiffs[-1])
    	dtime = before_tdiffs[-1]/60.0
    else:
    	dtime = default_dtime

    print 'dtime: {} minutes'.format(dtime)

    file_in_times = file_time in radar_times
    file_after = (before_tdiffs>0).all() or len(before_tdiffs) == 0

    #print 'file in times: {}, file after: {}'.format(file_in_times, file_after)

    if not file_in_times and pargs.realtime:
    	print 'time not in saved data, going to add it'
    	if file_after:
    	    print 'time is at the end'
    	    add_data_flag = True
    	    radar_times = np.append(radar_times, file_time)
    	else:
    	    print 'file is in between radar times, might try to insert it?'

    else:
    	add_data_flag = False


    lat = radar.data['lat0'][:]
    lon = radar.data['lon0'][:]
    z = radar.data['z0'][:]

    z_val = ycfg['z_val']
    z_ind = np.argmin(np.abs(z-z_val))

    z_top = ycfg['z_top']
    z_bot = ycfg['z_bot']

    z_top_ind = np.argmin(np.abs(z-z_top))
    z_bot_ind = np.argmin(np.abs(z-z_bot))

    radar_lat_array = radar.data[radar.lat_name]
    radar_lon_array = radar.data[radar.lon_name]

    #print 'z index: {}'.format(z_ind)

    # Back out location of the radar using the x/y/lat/lon info
    x_0 = np.where(radar.data[radar.x_name] == 0)[0][0]
    y_0 = np.where(radar.data[radar.y_name] == 0)[0][0]

    radar_lat = radar_lat_array[y_0, x_0]
    radar_lon = radar_lon_array[y_0, x_0]

    radar_dx = radar.dx
    radar_dy = radar.dy

   
    nx = len(radar.data[radar.x_name])
    ny = len(radar.data[radar.y_name])
    nz = len(radar.data[radar.z_name])
    #print 'nx: {}, ny: {}, nz: {}, dx: {}, dy: {}, dz: {}'.format(nx, ny, nz, radar.dx, radar.dy, radar.dz)
    xrg = [radar.data[radar.x_name].min(), radar.data[radar.x_name].max()]
    yrg = [radar.data[radar.y_name].min(), radar.data[radar.y_name].max()]
    zrg = [radar.data[radar.z_name].min(), radar.data[radar.z_name].max()]
    #print 'xrange: {}, yrange: {}, zrange: {}'.format(xrg, yrg, zrg)



    # These are the x (lon) and y (lat) limits for plotting later on
    xl = [radar_lon-lon_width/2.0, radar_lon+lon_width/2.0]
    yl = [radar_lat-lat_width/2.0, radar_lat+lat_width/2.0]

    xl_zoom = [radar_lon-lon_width_zoom/2.0, radar_lon+lon_width_zoom/2.0]
    yl_zoom = [radar_lat-lat_width_zoom/2.0, radar_lat+lat_width_zoom/2.0]


    # this finds the index in the larger array for putting the current radar data based on the current location
    i_rad_lat = np.argmin(np.abs(radar_lat-big_lat_array))
    i_rad_lon = np.argmin(np.abs(radar_lon-big_lon_array))



    if ycfg['rain_method'] == 'slant': # Here is where the calc for rainfall is
    		# This is the slant/stepped approach here

    	dbz_slant = np.ma.masked_invalid(radar_slant_values(radar.data['x0'], radar.data['y0'], z, radar.data[radar.dz_name]))
    	zdr_slant = np.ma.masked_invalid(radar_slant_values(radar.data['x0'], radar.data['y0'], z, radar.data[radar.zdr_name]))
    	rho_slant = np.ma.masked_invalid(radar_slant_values(radar.data['x0'], radar.data['y0'], z, radar.data[radar.rho_name]))
    	kdp_slant = np.ma.masked_invalid(radar_slant_values(radar.data['x0'], radar.data['y0'], z, radar.data[radar.kdp_name]))
    	hid_slant = np.ma.masked_invalid(radar_slant_values(radar.data['x0'], radar.data['y0'], z, radar.data[radar.hid_name]))
    #cs_slant = np.ma.masked_invalid(radar_slant_values(radar.data['x0'], radar.data['y0'], z, cs_arr))
    	try:

    	    cs_raintype, cs_rtypes = raintype.raintype(dbz_slant, refl_missing_val=dbz_slant.min(), 
                                   refl_dx=1.0, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                   shallowconvmin=shallowconvmin, truncZconvthres=truncZconvthres,
                                   dBZformaxconvradius=dBZformaxconvradius,
                                   weakechothres=weakechothres, backgrndradius=backgrndradius,
                                   maxConvRadius=maxConvRadius, minsize=minsize,
                                   startslope=startslope, maxsize=maxsize)

    # The output of this is a little weird, it outputs numbers from 0-9 and then a key for what they mean
    # Our rainfall algorithm only takes convective, stratiform and mixed so we need to convert
    # the output to output we can use

    # ***** Probably need to convert his categories into ours, maybe here before we plot??
    # ISO_CONV_CORE: 4, 'CONVECTIVE': 2, ISO_CS_CORE: 9 ARE CONVECTIVE WHICH WILL BE ASSINGED A 1
    # WEAK_ECHO: 6, STRATIFORM: 1, ISO_CONV_FRINGE: 5 ARE STRATIFORM AND GET A -1
    # UNCERTAIN: 3 IS mixed and will get a 0

	    # first doing the convection
	    conv = (cs_raintype == 4) | (cs_raintype == 2) | (cs_raintype == 9)
	    strat = (cs_raintype == 6) | (cs_raintype == 1) | (cs_raintype == 5)
	    mixed = cs_raintype == 3


	    cs_arr = deepcopy(cs_raintype)
	    cs_arr[conv] = 2
	    cs_arr[strat] = 1
	    cs_arr[mixed] = 3

	    cs_arr[dbz_slant<0] = 0


    	except Exception as cse:
    	    print 'Error with C/S classification: {}'.format(cse)
    	    logging.error(cse, exc_info=True)

    	    cs_arr = np.zeros(dbz_slant.shape, np.float)


   
    	fig, ax = plt.subplots(3, 2, figsize=(8,10))
    	axf = ax.flatten()

    	axf[0].pcolormesh(lon, lat, dbz_slant, cmap=radar.cmaps[radar.dz_name], vmin=0, vmax=80)
    	axf[0].set_title('Slanted reflectivity')
    	axf[1].pcolormesh(lon, lat, radar.data[radar.dz_name][z_ind], cmap=radar.cmaps[radar.dz_name], vmin=0, vmax=80)
    	axf[1].set_title('CAPPI reflectivity: %.1f km'%(z_val))

    	axf[2].pcolormesh(lon, lat, zdr_slant, cmap=radar.cmaps[radar.zdr_name], vmin=-1, vmax=3)
    	axf[2].set_title('Slanted ZDR')
    	axf[3].pcolormesh(lon, lat, radar.data[radar.zdr_name][z_ind], cmap=radar.cmaps[radar.zdr_name], vmin=-1, vmax=3)
    	axf[3].set_title('CAPPI ZDR: %.1f km'%(z_val))

    	axf[4].pcolormesh(lon, lat, rho_slant, cmap=radar.cmaps[radar.rho_name], vmin=0.7, vmax=1.0)
    	axf[4].set_title('Slanted rho')
    	axf[5].pcolormesh(lon, lat, radar.data[radar.rho_name][z_ind], cmap=radar.cmaps[radar.rho_name], vmin=0.7, vmax=1.0)
    	axf[5].set_title('CAPPI rho: %.1f km'%(z_val))

    	plt.tight_layout()

    	plt.savefig('dbzslant.png')


    	trop_rr, trop_rr_meth = csu_blended_rain_tropical.calc_blended_rain_tropical(dz=dbz_slant, zdr=zdr_slant, kdp=kdp_slant, 
                        cs=cs_arr, fhc=hid_slant)

    	trop_rr = np.ma.masked_invalid(trop_rr)

    	ranges = [0.0, 30.0, 60.0, 80.0, 120.0]
    	heights = [0.5, 1.0, 1.5, 2.0, 2.5]


    	shi = slant_height_indices(radar.data['x0'], radar.data['y0'], z, ranges, heights)
    	slant_height = np.zeros(dbz_slant.shape, np.float)
    
    	for i_s in range(shi.shape[0]):
    	    for j_s in range(shi.shape[1]):
    	    	slant_height[i_s, j_s] = z[shi[i_s, j_s]]

    	fig, ax = plt.subplots(1, 1, figsize=(8,6))
    	shi_pc = ax.pcolormesh(lon, lat, slant_height)
    	shi_cb = plt.colorbar(shi_pc, ax=ax)

    	ax.scatter(radar_lon, radar_lat, c='black', s=50)

    	plot_range_rings(ax, radar_lat, radar_lon, major_circle_rads, minor_rads=minor_circle_rads)

    	ax.set_aspect('equal')

    	ax.set_xlabel('Lon')
    	ax.set_ylabel('Lat')
    
    	ax.set_xlim(*xl)
    	ax.set_ylim(*yl)

    	ax.set_title('Slanted beam height rain calculation')

    	plt.tight_layout()

    	plt.savefig('zindarray.png')


    elif ycfg['rain_method'] == 'composite': # Here is where the composite approach will be used

    	z_top_ind = np.argmin(np.abs(z-z_top))
    	z_bot_ind = np.argmin(np.abs(z-z_bot))


    	dbz = radar.data[radar.dz_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]
    	zdr = radar.data[radar.zdr_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]
    	kdp = radar.data[radar.kdp_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]
    	rho = radar.data[radar.rho_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]

    	hid = radar.data[radar.hid_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]

	dbz_comp = radtools.composite(dbz, bad_val=-9.9, func=np.max)
	print 'composite dBZ shape: {}'.format(dbz_comp.shape)


    	try:

    	    cs_raintype, cs_rtypes = raintype.raintype(dbz_comp, refl_missing_val=dbz.min(), 
                                   refl_dx=1.0, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                   shallowconvmin=shallowconvmin, truncZconvthres=truncZconvthres,
                                   dBZformaxconvradius=dBZformaxconvradius,
                                   weakechothres=weakechothres, backgrndradius=backgrndradius,
                                   maxConvRadius=maxConvRadius, minsize=minsize,
                                   startslope=startslope, maxsize=maxsize)



    # The output of this is a little weird, it outputs numbers from 0-9 and then a key for what they mean
    # Our rainfall algorithm only takes convective, stratiform and mixed so we need to convert
    # the output to output we can use

    # ***** Probably need to convert his categories into ours, maybe here before we plot??
    # ISO_CONV_CORE: 4, 'CONVECTIVE': 2, ISO_CS_CORE: 9 ARE CONVECTIVE WHICH WILL BE ASSINGED A 1
    # WEAK_ECHO: 6, STRATIFORM: 1, ISO_CONV_FRINGE: 5 ARE STRATIFORM AND GET A -1
    # UNCERTAIN: 3 IS mixed and will get a 0

	    # first doing the convection
	    conv = (cs_raintype == 4) | (cs_raintype == 2) | (cs_raintype == 9)
	    strat = (cs_raintype == 6) | (cs_raintype == 1) | (cs_raintype == 5)
	    mixed = cs_raintype == 3


	    cs_arr = deepcopy(cs_raintype)
	    cs_arr[conv] = 2
	    cs_arr[strat] = 1
	    cs_arr[mixed] = 3

	    cs_arr[dbz_comp<0] = 0


    	except Exception as cse:
    	    print 'Error with C/S classification: {}'.format(cse)
    	    logging.error(cse, exc_info=True)

    	    cs_arr = np.zeros(dbz.shape, np.float)


        #print 'C/S max: {}'.format(cs_arr.max())
        #print 'C/S shape: {}'.format(cs_arr.shape)

    	trop_rr, trop_rr_meth = csu_blended_rain_tropical.calc_blended_rain_tropical(dz=dbz, zdr=zdr, kdp=kdp, 
                        cs=cs_arr, fhc=hid)


    	#print 'orig RR shape: {}'.format(trop_rr.shape)

    	trop_rr = np.ma.masked_invalid(trop_rr)
    	trop_rr = np.ma.max(trop_rr, axis=0)

    	#print 'final RR shape: {}'.format(trop_rr.shape)

    	pass

    else:
    	print 'Enter either slant or composite for the rain_method in the yaml config file'

    # Adding the convective stratiform field to the radar object
    radar.add_field(cs_arr, 'CS')


    # ****** NOW PLOTTING STUFF **********

    # this will do the 6 panel plot of the variables in the varlist

    fig, ax = plt.subplots(3,2, figsize=(9, 10.5))

    fig, ax = radar.cappi_multiplot(z=z_val, varlist=[radar.dz_name, radar.zdr_name, radar.kdp_name, radar.rho_name, 
    				'CS', 'HID'], coords='ll', xlim=xl, ylim=yl, fig=fig, ax=ax)

    axf = ax.flatten()
#    major_circle_rads = np.arange(20, 100, 20)
#    minor_circle_rads = np.arange(10, 110, 20)


    #print 'Radar Lat: %.2f, Lon: %.2f'%(radar_lat, radar_lon)

    for a in axf:

	axt = a.get_xticks().tolist()
	for j in range(len(axt)):
	    axt[j] = gtools.convert_decimal_to_degree(axt[j])
	a.set_xticklabels(axt)

	axty = a.get_yticks().tolist()
	#print 'y tick labels: {}'.format(axty)
	for j in range(len(axty)):
	    axty[j] = gtools.convert_decimal_to_degree(axty[j])
	a.set_yticklabels(axty)


        a.scatter(radar_lon, radar_lat, s=25, facecolor='black')
        a.grid(True)

	for label in a.get_xticklabels():
            label.set_rotation(25) 

	for label in a.get_yticklabels():
            label.set_rotation(25) 


	a.tick_params(labelsize=8)
	a.set_xlabel('')
	a.set_ylabel('')

        ptools.plot_range_rings(a, radar_lat, radar_lon, major_circle_rads, minor_rads=minor_circle_rads)
        ptools.plot_azimuths(a, radar_lat, radar_lon, azimuths)

	a.set_aspect('equal')


    plt.tight_layout()
    plt.subplots_adjust(top=0.93)

    plt.savefig('%s/figures/cappi/%s_variables.%s'%(base_path, dt_string, image_type), dpi=DPI)
#    print 'figure saved: {}'.format(time.time() - start)

    plt.close(fig)

    ############### Self consistency plot ##############################
    try:

        fig, ax = plt.subplots(1,2, figsize=(10,5))

        ax[0].scatter(radar.data['DZQC'].ravel(), radar.data['KDP'].ravel(), c='black', alpha=0.3, s=8, edgecolors='none')

        ax[0].set_xlabel('Z')
        ax[0].set_ylabel('Kdp')

        ax[0].set_xlim(0, 65)
        ax[0].set_ylim(-1, 7)

        # Second plot will be ZDR/Kdp
        ax[1].scatter(radar.data[radar.zdr_name].ravel(), radar.data['KDP'].ravel(), c='black', alpha=0.3, s=8, edgecolors='none')
        ax[1].set_xlabel('ZDR')
        ax[1].set_ylabel('Kdp')

        ax[1].set_xlim(0, 6)
        ax[1].set_ylim(0, 7)



        for a in ax:
            a.grid(True)

        fig.suptitle('SEAPOL: %s self consistency'%(file_dt_string), fontsize=16)
        plt.tight_layout()
        fig.subplots_adjust(top=0.92)
    
        plt.savefig('%s/figures/self_consistency/%s_z_kdp.%s'%(base_path, dt_string, image_type), dpi=DPI)


    except Exception, sce:
        print 'Error with self consistency plot: {}'.format(sce)



    quit()


    #######################################################33



    # ***************** Going to delete the radar object here to try to save memory
    #del radar
    gc.collect() # Doing some garbage collecting here to try to save some memory as well.



#     if len(radar_lats) > 1:

#             dlat = np.diff(radar_lats)[-1]
#             dlon = np.diff(radar_lons)[-1]
#             dy = dlat*111.0
#             dx = dlon*(111.0*np.cos(np.radians(radar_lats[-1])))
#             speed_km = np.sqrt(dx**2 + dy**2)*60.0/dtime
#             speed = km2nm*speed_km

#     # figure out a heading and throw it on there

#             heading_physics = np.degrees(np.arctan2(dy, dx))
#             heading = 90. - heading_physics
#             if heading < 0.:
#                 heading += 360.0

#         else:
#             speed = 0.
#             heading = 0.
# #
    
#         print 'speed: {}, heading: {}'.format(speed, heading)




##############################################################################################
##############################################################################################
    if True:

        try: 
        
            fig, ax = plt.subplots(1,1, figsize=(8, 6))
        #axf = ax.flatten()
            #print 'big lon array: {}, big lat array: {}, acc rain: {}'.format(big_lon_array.shape, 
                                                       # big_lat_array.shape, accum_rainfall_good_times_long.shape)
            #print 'working acc rain: {}'.format(accum_rainfall_good_times1.shape)

            #print 'radar lon array: {}, radar lat array: {}, cs array: {}'.format(radar_lon_array.shape, radar_lat_array.shape,
                                                        #cs_arr.shape)

            argtl = accum_rainfall_good_times_long[:len(big_lat_array), :len(big_lon_array)].filled(0.0)

            argtl_sm = ndi.filters.gaussian_filter(argtl, 2)

            fit_lon_array = big_lon_array[:argtl.shape[1]]

            #print 'argtl: {}'.format(argtl.shape)


#            acc_ct = ax.contourf(fit_lon_array, big_lat_array, argtl_sm, levels=np.arange(25, 250, 50), 
#                                colors='black', linewidths=0.8)

            #acc_ct = ax.contourf(fit_lon_array, big_lat_array, argtl)

            #acc_cb = plt.colorbar(acc_ct, ax=ax[0])


            accs_ct = ax.contour(fit_lon_array, big_lat_array, argtl_sm, levels=np.arange(20, 200, 20), 
                                        colors='black', linewidths=0.8)
            #accs_cb = plt.colorbar(accs_ct, ax=ax[1])

            cs_ma = np.ma.masked_where(cs_arr != 2, cs_arr)
            print 'cs ma max: {}'.format(cs_ma.max())
            cs_ct = ax.pcolormesh(radar_lon_array, radar_lat_array, cs_ma, cmap=radar.cmaps['CS'], vmin=0, vmax=3)


            #for a in ax:

            ax.set_xlim(*xl)
            ax.set_ylim(*yl)
            ax.grid(True)

            plot_range_rings(ax, radar_lat, radar_lon, small_circle_rads)
            plot_azimuths(ax, radar_lat, radar_lon, azimuths)
            ax.set_aspect('equal')

            fig.suptitle('SEAPOL: %s Puddles and convective cells'%(file_dt_string), fontsize=16)

            plt.tight_layout()
            fig.subplots_adjust(top=0.92)

            plt.savefig('%s/figures/puddle_analysis/%s_accum_cs.%s'%(base_path, dt_string, image_type), dpi=DPI)
        
        except Exception, e2:
            print 'Error with plotting puddle stuff: {}'.format(e2)
            logging.error(e2, exc_info=True)




##############################################################################################
##############################################################################################


    	total_time_since_beginning_of_list = (file_time-radar_times[0]).total_seconds()/60.0
    	#print 'after tdiffs max: {}'.format(after_tdiffs.max())
    	#print 'time since beginning: {}'.format(total_time_since_beginning_of_list)

    	try:

    	    if total_time_since_beginning_of_list > ycfg['save_data_thresh']*60.0: # save data thresh is in hours
    	    	print 'Need to chop off some times from the saved data'
    	    # now need to figure out where to cut it off, which index
    	    	after_tdiff_minutes = after_tdiffs/60.0
    	    	print 'length of radar times: {}'.format(len(radar_times))
    	    	print 'beginning: {}'.format(after_tdiff_minutes[:10])
    	    	past_limit_ind = np.where(after_tdiff_minutes >= ycfg['save_data_thresh']*60.0)[0]
    	    	print 'past limit indices: {}'.format(past_limit_ind)
    	    	if len(past_limit_ind) > 0: # this means we need to remove some data on the front
    	    	    remove_ind = past_limit_ind[-1]
    	    	    print 'Need to remove data past index: {}'.format(remove_ind)
    	    	    radar_times = radar_times[remove_ind:]
    	    	    rainfall = rainfall[remove_ind:]
    	    	    radar_lats = radar_lats[remove_ind:]
    	    	    radar_lons = radar_lons[remove_ind:]

    	except Exception, de:
    	    print 'Error with trying to remove old data: {}'.format(de)

    	np.savez_compressed('%s/%s/radar_times_accum.npz'%(base_path, ycfg['saved_data_path']), times=radar_times, rainfall=rainfall, 
    					radar_lats=radar_lats, radar_lons=radar_lons)


    	try:
    	    print 'Gonna try to back up the saved data'
    	    backup_path = '/home/rmet/SPURS2/data_backup'
    	    bp_check = sorted(glob.glob('%s/*.npz'%(backup_path)))
    	    print 'len of dataup backup check: {}'.format(len(bp_check))
    	    bsave_flag = False # set this up as the default

    	    if len(bp_check) == 0:
    	    	print 'No backups'
    	    	bsave_flag = True
    	    else:
    	    	# if there are files, need to check them for times
    	    	save_file_time_strings = np.array([os.path.basename(_) for _ in bp_check])
    	    	print 'save files: {}'.format(save_file_time_strings[-1])
    	    	sskip = len('radar_times_accum_')
    	    	sf_times = np.array([ datetime.datetime.strptime(_[sskip:sskip+15], '%Y%m%d-%H%M%S') for _ in save_file_time_strings])
    	    	sf_time_diffs = np.array([ (file_time-_).total_seconds() for _ in sf_times ])

    	    	if sf_time_diffs.min() > 2.0*3600.0: # if it's been more than 12 hours, save a new file
    	    	    bsave_flag = True

    	    if bsave_flag:

    	    	os.system('cp -v %s/%s/radar_times_accum.npz %s/radar_times_accum_%s.npz'%(
    	    				base_path, ycfg['saved_data_path'], backup_path, file_time.strftime('%Y%m%d-%H%M%S')))


    	except Exception, bue:
    	    print 'Error with trying to back up: {}'.format(bue)
            logging.error(bue, exc_info=True)


    print 'figures and files saved'

   # ********* Now pick a radius? or a range of x and y values to check max rain rate and accum rainfall and age

    # gonna be centered around i_rad_lat and i_rad_lon
    # put in a giant try except just so it doesn't choke out the rest of this!
    try:
    	rain_search_range = ycfg['rain_search_range'] # km
    	rain_sr_lon = int(np.round(rain_search_range/radar_dx))
    	rain_sr_lat = int(np.round(rain_search_range/radar_dy))


    	near_radar_rain = big_rain_array[i_rad_lat-rain_sr_lat:i_rad_lat+rain_sr_lat, \
    			i_rad_lon-rain_sr_lon:i_rad_lon+rain_sr_lon]

    	accum_search_range = ycfg['accum_search_range'] # km
    	accum_sr_lon = int(np.round(accum_search_range/radar_dx))
    	accum_sr_lat = int(np.round(accum_search_range/radar_dy))

    	near_radar_accum = accum_rainfall_good_times2[i_rad_lat-accum_sr_lat:i_rad_lat+accum_sr_lat, \
    			i_rad_lon-accum_sr_lon:i_rad_lon+accum_sr_lon]

    	rain_value = near_radar_rain.max()/(default_dtime/60.0) 
    	accum_value = near_radar_accum.max()

    	print 'near ship rain max: %.1f, accum max: %1f'%(rain_value, accum_value)    			

    	if ( (rain_value >= ycfg['rain_search_thresh']) or (accum_value >= ycfg['accum_search_thresh']) ):
    	# basically, if it's raining or rain has accumulated near the ship, log it
    	    print 'Going thru some rain or a freshwater puddle!'
    	    logfile = open('/home/rmet/SPURS2/logs/ship/%s.log'%(file_time.strftime('%Y%m%d')), 'a+')
    	    logfile.write('%s\t%.2f\t%.2f\n'%(file_time.strftime('%Y%m%d-%H%M%S'), rain_value, accum_value))
    	    logfile.close()
    	    pass

    except Exception, loge:
    	print 'ERROR with the ship rainfall logging: {}'.format(loge)




print 'Made it to the end of the processing'

        # src3d = src3d[np.newaxis, ...]

        # for key in sorted(src_hist.keys()[1:]):
        #     src3d = np.concatenate((src3d,src_hist[key][np.newaxis,...]), axis=0)






