# Brody Fuchs, Sept 2017
# brfuchs@atmos.colostate.edu

# Just some code to get started reading in radar data and doing some plotting and rainrate stuff
# This version will use the csuram stuff to make it easier

# March 2018: This code will be used with the DROPS code and will do plotting of gridded data
# and do instantaneous rainfall calculations and save rainfall data, but will not do the 
# aggregation of rainfall data for things like 30/60 minute accumulation


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
from skewPy import SkewT
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
import cPickle as pickle
import steiner_houze_yuter_cs as shy
from csuram import analysis_tools as atools



start = time.time()

try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception:
    base_path = os.path.dirname(os.path.realpath('__file__'))



radar_path = '%s/realtime_gridded'%(base_path)
file_type = 'nc'

parser = argparse.ArgumentParser(description='Put in a file to be processed')

#parser.add_argument('--noarg', action="store_true", default=False)
parser.add_argument('--file', action="store", dest="file", default=None)
parser.add_argument('--lastfile', action="store", dest="lastfile", type=int, default=0)
parser.add_argument('--allfiles', action="store", dest="allfiles", type=int, default=0)
parser.add_argument('--wildcard', action="store", dest="wildcard", type=str, default='')
parser.add_argument('--realtime', action="store", dest="realtime", type=int, default=0)
parser.add_argument('--analysis', action="store", dest="analysis", type=int, default=0)
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
    radar_files = sorted(glob.glob('%s/*%s*.%s'%(radar_path, pargs.wildcard, file_type)))

elif pargs.lastfile:
    radar_files = [glob.glob('%s/*.%s'%(radar_path, file_type))[-1]]

else:
    print 'No files specified or found, exiting.'
    sys.exit()

print 'len radar files: {}'.format(len(radar_files))


# ******* Thresholds and settings ****************

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



image_type = 'png' # of tk looper only handles pngs
km2nm = 0.54 # 1 km is 0.54 nautical miles
nm2km = 1.0/km2nm # Just to have both conversions

# ****** Colorbar options ***************************
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
prefix = ycfg['qc_prefix']


fig_path = ycfg['fig_path']


# Here we're just going to use a sounding I ripped off from Hawaii on UWYO's website
# We need the sounding for the HID algorithm as well as the rainfall algorithm (which uses HID, among other things)
# We just need to know about where the freezing level is, which won't change more than a couple hundred meters
# ******** Need to read this file and get it into the HID calculation **********

snd_data = np.genfromtxt('%s/sample_tropical_sounding.txt'%(base_path), dtype=None, skip_header=6)
snd_ht = snd_data['f1'].astype(float)
snd_temp = snd_data['f2']
snd_pres = snd_data['f0']
snd_dwpt = snd_data['f3']

snd = SkewT.Sounding(data={'temp': snd_temp, 'hght': snd_ht, 'dwpt': snd_dwpt, 'pres': snd_pres})


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
    	radar = RadarData.RadarData(radar_file=rf, x='x0', y='y0', z='z0', band=ycfg['radar_band'], squeeze=True, 
                    lat='lat0', lon='lon0', dz=ycfg['dz_name'], kdp=ycfg['kdp_name'], zdr=ycfg['zdr_name'], rho=ycfg['rho_name'])

    	    # Here we add the sounding object to the radar object
    	radar.add_sounding_object(snd)
    	# Interpolate the sounding heights to the radar heights
    	radar.interp_sounding()
    	# Calculate the HID and add that info to the radar object
    	#radar.set_hid(band=ycfg['radar_band'], use_temp=True) # this calls the appropriate functions and uses the sounding that we've
        this_hid = radar.get_hid(band=ycfg['radar_band'], use_temp=True)
        # attributed to the radar object and uses the proper radar band

    except KeyError, ke:
    	print 'Probably started out as a binary file with different names: {}'.format(ke)
    	radar = RadarData.RadarData(radar_file=rf, x='x0', y='y0', z='z0', band='C', squeeze=True, 
                    lat='lat0', lon='lon0', dz='DZQC', zdr='ZDR2', kdp='KDP', rho='RHOHV2', hid='HID')    	

    	radar.add_sounding_object(snd)
    	# Interpolate the sounding heights to the radar heights
    	radar.interp_sounding()
    	# Calculate the HID and add that info to the radar object
    	this_hid = radar.get_hid(band=ycfg['radar_band'], use_temp=True) # this calls the appropriate functions and uses the sounding that we've
        # attributed to the radar object and uses the proper radar band

    this_hid[radar.data[radar.dz_name].mask] = -1

    radar.add_field(this_hid, 'HID')

    map_res = radar.dx # km

    # Just grabbing the time of the radar file from the name and appending it to the accum_times list
    last_dot = rf.rfind('.')
    d_index = rf.rfind('d')

    dt_string = rf[d_index+2:d_index+17]
    print dt_string

    file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
    radar.date = file_time

    file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')
    date_string = file_time.strftime('%Y%m%d')
#    radar_times = np.zeros(0, object)


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


    if ycfg['rain_beam_method'] == 'slant': # Here is where the calc for rainfall is
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
	    traceback.print_exc()
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


    	rr_array, rr_meth = csu_blended_rain_tropical.calc_blended_rain_tropical(dz=dbz_slant, zdr=zdr_slant, kdp=kdp_slant, 
                        cs=cs_arr, fhc=hid_slant)

    	rr_array = np.ma.masked_invalid(rr_array)

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


    elif ycfg['rain_beam_method'] == 'composite': # Here is where the composite approach will be used

    	z_top_ind = np.argmin(np.abs(z-z_top))
    	z_bot_ind = np.argmin(np.abs(z-z_bot))


    	dbz = radar.data[radar.dz_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]
    	zdr = radar.data[radar.zdr_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]
    	kdp = radar.data[radar.kdp_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]
    	rho = radar.data[radar.rho_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]

    	hid = radar.data[radar.hid_name][:].squeeze()[z_bot_ind:z_top_ind, :, :]

	dbz_comp = radtools.composite(dbz, bad_val=-9.9, func=np.max)
	print 'composite dBZ shape: {}'.format(dbz_comp.shape)

	print 'C/S z height: {}'.format(ycfg['cs_z_val'])

        z_cs_ind = np.argmin(np.abs(z-ycfg['cs_z_val']))
        dbz_cs_lev = radar.data[radar.dz_name][z_cs_ind]

    	try:

	    yh_cs, yh_cc, yh_bkgnd = shy.conv_strat_latlon(dbz_cs_lev, lat, lon, 40.0, method='SYH', a=8, b=64)

	    yh_cc_ma = np.ma.masked_where(yh_cc < 0, yh_cc)


	    radar.add_field(yh_cs, 'SHYCS')
	    radar.add_field(yh_cc, 'SHYCC')
	    radar.add_field(yh_bkgnd, 'SHYBG')

	    radar.add_field(dbz_cs_lev-yh_bkgnd, 'DBZdiff')


	    # let's make 3-5 convective, 1-2 mixed, and 0 stratiform
        # This is where we could make the change to either include mixed as a separate category
        # or just move it into convective

	    cs_arr = np.full(yh_cs.shape, np.nan)
	    yh_conv = (yh_cs == 3) | (yh_cs == 4) | (yh_cs == 5)
	    yh_mixed = (yh_cs == 1) | (yh_cs == 2)
	    yh_strat = yh_cs == 0

	    cs_arr[yh_conv] = 2
	    cs_arr[yh_mixed] = 3
	    cs_arr[yh_strat] = 1


    	except Exception as cse:
    	    print 'Error with C/S classification: {}'.format(cse)
    	    logging.error(cse, exc_info=True)
	    traceback.print_exc()
    	    cs_arr = np.zeros(dbz.shape, np.float)


        # Need to put a check in here for rrmethod

        if 'rainmethod' in ycfg.keys():
            print 'Rain method: {}'.format(ycfg['rainmethod'])
        else:
            print 'No rain method selected, defaulting to csublended'

        if 'rainmethod' in ycfg.keys() and ycfg['rainmethod'] == 'csutropical':
    	    rr_array, rr_meth = csu_blended_rain_tropical.calc_blended_rain_tropical(dz=dbz, zdr=zdr, kdp=kdp, 
                        cs=cs_arr, fhc=hid)

        elif 'rainmethod' in ycfg.keys() and ycfg['rainmethod'] == 'csuhidro':
            rr_array, rr_meth = csu_blended_rain.csu_hidro_rain(dz=dbz, zdr=zdr, kdp=kdp, fhc=hid, band=ycfg['radar_band'],
                        thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3, thresh_nexrad=53.0)


        elif 'rainmethod' in ycfg.keys() and ycfg['rainmethod'] == 'csublended':
            rr_array, rr_meth = csu_blended_rain.calc_blended_rain(dz=dbz, zdr=zdr, kdp=kdp, ice_flag=False, band=ycfg['radar_band'],
                        thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3)

        else: # if nothing specified, just go with the default of csublended
            rr_array, rr_meth = csu_blended_rain.calc_blended_rain(dz=dbz, zdr=zdr, kdp=kdp, ice_flag=False, band=ycfg['radar_band'],
                        thresh_dz=38.0, thresh_zdr=0.5, thresh_kdp=0.3)



    	rr_array = np.ma.masked_invalid(rr_array)
    	rr_array = np.ma.max(rr_array, axis=0)


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

        ptools.plot_range_rings(a, radar_lat, radar_lon, major_circle_rads, minor_rads=minor_circle_rads, text_flag=False)
        ptools.plot_azimuths(a, radar_lat, radar_lon, azimuths)

	a.set_aspect('equal')


    plt.tight_layout()
    plt.subplots_adjust(top=0.93)



    #plt.savefig('%s/figures/cappi/%s_variables.%s'%(base_path, dt_string, image_type), dpi=DPI)
#    print 'figure saved: {}'.format(time.time() - start)

    ptools.save_figure(fig, '%s/cappi/%s'%(fig_path, date_string), '%s%s_variables.%s'%(prefix, dt_string, image_type), dpi=DPI)
    plt.close(fig)

    ############### Self consistency plot ##############################
    z_rain_top_ind = np.argmin(np.abs(z-4.0))
    z_rain_bot_ind = np.argmin(np.abs(z-1.0))

    try:

        fig, ax = plt.subplots(2,2, figsize=(10,10))
	axf = ax.flatten()

	logmin = 5e-4
	logmax = 100

        dz_bins = np.arange(0, 70, 1.25)
        kdp_bins = np.arange(-1, 7, 0.25)
        zdr_bins = np.arange(-1, 6, 0.25)

        dz_kdp_hist = np.histogramdd((radar.data[radar.dz_name].ravel(), radar.data[radar.kdp_name].ravel()), 
                                        bins=(dz_bins, kdp_bins), normed=True)[0]
        dz_kdp_hist *= 100.0/dz_kdp_hist.sum()
        dk_pc = axf[0].pcolormesh(dz_bins[:-1], kdp_bins[:-1], dz_kdp_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        dk_cb = plt.colorbar(dk_pc, ax=axf[0])

        #ax[0].scatter(radar.data[radar.dz_name].ravel(), radar.data[radar.kdp_name].ravel(), c='black', alpha=0.3, s=8, edgecolors='none')

        axf[0].set_xlabel('Z')
        axf[0].set_ylabel('Kdp')
	axf[0].set_title('All points Z vs Kdp')

        zdr_kdp_hist = np.histogramdd((radar.data[radar.zdr_name].ravel(), radar.data[radar.kdp_name].ravel()), 
                                        bins=(zdr_bins, kdp_bins), normed=True)[0]
        zdr_kdp_hist *= 100.0/zdr_kdp_hist.sum()
        zk_pc = axf[1].pcolormesh(zdr_bins[:-1], kdp_bins[:-1], zdr_kdp_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        zk_cb = plt.colorbar(zk_pc, ax=axf[1])

        # Second plot will be ZDR/Kdp
        axf[1].set_xlabel('ZDR')
        axf[1].set_ylabel('Kdp')
	axf[1].set_title('All points ZDR vs Kdp')

        #ax[1].set_xlim(0, 6)
        #ax[1].set_ylim(0, 7)

	# ********** Now let's do the raining points only *************

    	z_rain_top_ind = np.argmin(np.abs(z-4.0))
    	z_rain_bot_ind = np.argmin(np.abs(z-1.0))


        dz_kdp_rain_hist = np.histogramdd((radar.data[radar.dz_name][z_rain_bot_ind:z_rain_top_ind].ravel(), 
						radar.data[radar.kdp_name][z_rain_bot_ind:z_rain_top_ind].ravel()), 
                                        bins=(dz_bins, kdp_bins), normed=True)[0]
        dz_kdp_rain_hist *= 100.0/dz_kdp_rain_hist.sum()
        dkr_pc = axf[2].pcolormesh(dz_bins[:-1], kdp_bins[:-1], dz_kdp_rain_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        dkr_cb = plt.colorbar(dkr_pc, ax=axf[2])

        axf[2].set_xlabel('Z')
        axf[2].set_ylabel('Kdp')
	axf[2].set_title('Raining points')


        zdr_kdp_rain_hist = np.histogramdd((radar.data[radar.zdr_name][z_rain_bot_ind:z_rain_top_ind].ravel(), 
				radar.data[radar.kdp_name][z_rain_bot_ind:z_rain_top_ind].ravel()), 
                                        bins=(zdr_bins, kdp_bins), normed=True)[0]
        zdr_kdp_rain_hist *= 100.0/zdr_kdp_rain_hist.sum()
        zkr_pc = axf[3].pcolormesh(zdr_bins[:-1], kdp_bins[:-1], zdr_kdp_rain_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        zkr_cb = plt.colorbar(zk_pc, ax=axf[3])

        # Second plot will be ZDR/Kdp
        axf[3].set_xlabel('ZDR')
        axf[3].set_ylabel('Kdp')
	axf[3].set_title('Raining points')



        for a in axf:
            a.grid(True)

        fig.suptitle('%s: %s self consistency'%(radar_title, file_dt_string), fontsize=16)
        plt.tight_layout()
        fig.subplots_adjust(top=0.92)
    
        #plt.savefig('%s/figures/self_consistency/%s_z_kdp.%s'%(base_path, dt_string, image_type), dpi=DPI)
    	ptools.save_figure(fig, '%s/self_consistency/%s'%(fig_path, date_string), '%s%s_z_kdp.%s'%(prefix, dt_string, image_type), dpi=DPI)


    except Exception, sce:
        print 'Error with self consistency plot: {}'.format(sce)
	traceback.print_exc()


    ############### Self consistency scatter plot ##############################
    try:

        fig, ax = plt.subplots(2,2, figsize=(10,10))
	axf = ax.flatten()
	scat_color = '0.5'

        axf[0].scatter(radar.data[radar.dz_name].ravel(), radar.data[radar.kdp_name].ravel(), c=scat_color, alpha=0.3, s=8, edgecolors='none')

        axf[0].set_xlabel('Z')
        axf[0].set_ylabel('Kdp')
	axf[0].set_title('All points Z vs Kdp')

	axf[0].set_xlim(0, 60)
	axf[0].set_ylim(-1, 6)


        axf[1].scatter(radar.data[radar.zdr_name].ravel(), radar.data[radar.kdp_name].ravel(), c=scat_color, alpha=0.5, s=8, edgecolors='none')

        # Second plot will be ZDR/Kdp
        axf[1].set_xlabel('ZDR')
        axf[1].set_ylabel('Kdp')
	axf[1].set_title('All points ZDR vs Kdp')

        axf[1].set_xlim(-1, 6)
        axf[1].set_ylim(-1, 6)

	# ********** Now let's do the raining points only *************




        axf[2].scatter(radar.data[radar.dz_name][z_rain_bot_ind:z_rain_top_ind].ravel(), 
						radar.data[radar.kdp_name][z_rain_bot_ind:z_rain_top_ind].ravel(), 
						c=scat_color, alpha=0.5, s=8, edgecolors='none')

        axf[2].set_xlabel('Z')
        axf[2].set_ylabel('Kdp')
	axf[2].set_title('Raining points')

	axf[2].set_xlim(0, 60)
	axf[2].set_ylim(-1, 6)

        axf[3].scatter(radar.data[radar.zdr_name][z_rain_bot_ind:z_rain_top_ind].ravel(), 
				radar.data[radar.kdp_name][z_rain_bot_ind:z_rain_top_ind].ravel(), 
				c=scat_color, alpha=0.5, s=8, edgecolors='none')


        axf[3].set_xlim(-1, 6)
        axf[3].set_ylim(-1, 6)

        # Second plot will be ZDR/Kdp
        axf[3].set_xlabel('ZDR')
        axf[3].set_ylabel('Kdp')
	axf[3].set_title('Raining points')



        for a in axf:
            a.grid(True)

        fig.suptitle('%s: %s self consistency'%(radar_title, file_dt_string), fontsize=16)
        plt.tight_layout()
        fig.subplots_adjust(top=0.92)
    
        #plt.savefig('%s/figures/self_consistency/%s_z_kdp.%s'%(base_path, dt_string, image_type), dpi=DPI)
    	ptools.save_figure(fig, '%s/self_consistency_scatter/%s'%(fig_path, date_string), 
			'%s%s_z_kdp.%s'%(prefix, dt_string, image_type), dpi=DPI)


    except Exception, sce:
        print 'Error with self consistency plot: {}'.format(sce)
	traceback.print_exc()

        # Now gonna make a similar plot as was done in some earlier experiments
        # where we plot Kdp (on a log scale) vs ZDR for conv/strat in different colors
    try:
#            cs_arr[conv] = 2 cs_arr[strat] = 1  cs_arr[mixed] = 3

        fig, ax = plt.subplots(1, 1, figsize=(9,9))
        conv_pt = np.where(cs_arr == 2)
        strat_pt = np.where(cs_arr == 1)
        mixed_pt = np.where(cs_arr == 3)

        kdp_rain = radar.data[radar.kdp_name][z_rain_bot_ind:z_rain_top_ind]
        zdr_rain = radar.data[radar.zdr_name][z_rain_bot_ind:z_rain_top_ind]

        # kdp_rain = radar.data[radar.kdp_name][0:4]
        # zdr_rain = radar.data[radar.zdr_name][0:4]

        kdp_rain_mask = kdp_rain.mask
        zdr_rain_mask = zdr_rain.mask 
        combo_mask = np.logical_or(kdp_rain_mask, zdr_rain_mask)


        # They need to have the same mask, otherwise the .ravel() won't have the same sizes
        kdp_rain.mask = combo_mask
        zdr_rain.mask = combo_mask

        conv_kdp = kdp_rain[:, conv_pt[0], conv_pt[1]].compressed().ravel()
        strat_kdp = kdp_rain[:, strat_pt[0], strat_pt[1]].compressed().ravel()
        mixed_kdp = kdp_rain[:, mixed_pt[0], mixed_pt[1]].compressed().ravel()

        conv_zdr = zdr_rain[:, conv_pt[0], conv_pt[1]].compressed().ravel()
        strat_zdr = zdr_rain[:, strat_pt[0], strat_pt[1]].compressed().ravel()
        mixed_zdr = zdr_rain[:, mixed_pt[0], mixed_pt[1]].compressed().ravel()


        strat_sc = ax.scatter(strat_kdp, strat_zdr, color='dodgerblue', edgecolors='none', alpha=0.25, s=20, label='stratiform')
        mixed_sc = ax.scatter(mixed_kdp, mixed_zdr, color='gold', edgecolors='none', alpha=0.25, s=20, label='mixed')
        conv_sc = ax.scatter(conv_kdp, conv_zdr, color='red', edgecolors='none', alpha=0.5, s=20, label='convective')


        ax.legend(loc='lower right', scatterpoints=1)

        ax.set_xlim(0.01, 10)
        ax.set_xscale('log')

        ax.set_xlabel('K$_{dp}$')
        ax.set_ylabel('Z$_{DR}$')


        ax.set_ylim(0, 2)


        ax.grid(True)

        fig.suptitle('%s: %s self consistency'%(radar_title, file_dt_string), fontsize=16)
        plt.tight_layout()
        fig.subplots_adjust(top=0.92)
    
        #plt.savefig('%s/figures/self_consistency/%s_z_kdp.%s'%(base_path, dt_string, image_type), dpi=DPI)
        ptools.save_figure(fig, '%s/self_consistency_scatter_cs/%s'%(fig_path, date_string), 
                        '%s%s_z_kdp_cs.%s'%(prefix, dt_string, image_type), dpi=DPI)



    except Exception, cse:
        print 'ERROR with convective/strat self consistency: {}'.format(cse)
	traceback.print_exc()


#####################################################################################
#####################################################################################
	# This is where we'll create the rainfall dictionary to be pickled for use later on

    rain_dict = {}
    rain_dict['rainfall'] = rr_array
    rain_dict['lat'] = lat
    rain_dict['lon'] = lon
    rain_dict['time'] = file_time

    rain_path = ycfg['rain_data_path']
	# make the path if it isn't there
    if len(glob.glob(rain_path)) == 0:
	os.system('mkdir -v %s'%(rain_path))

	# for some weird reason, lat and lon can be a masked array sometimes
    if isinstance(lat, np.ma.masked_array):
	lat = lat.filled(0.0)
	lon = lon.filled(0.0)

    np.savez_compressed('%s/%s_rain.npz'%(rain_path, dt_string), lat=lat, lon=lon, rainfall=rr_array.filled(0.0), 
    					time=file_time)

#####################################################################################
#####################################################################################

# Now just do some quick logging of basic variables so we can get an idea of where to start looking for cases


    if pargs.analysis:
 
	if rr_array.max() <= 1.0:
	    rrmax = 0.0
	else:
	    rrmax = np.percentile(rr_array[rr_array>0], 99)

	logfile = open('%s/radar_stats/%s_rr.log'%(base_path, date_string), 'a+')
	logfile.write('%s  %.2f\n'%(file_time.strftime('%Y%m%d-%H%M%S'), rrmax))
	logfile.close()

#####################################################################################


    # Now log the radar's lat/lon position using "radar_lat" and "radar_lon"
	print 'Logging ship lat/lon'


	ll_logfile = open('%s/radar_stats/%s_latlon.log'%(base_path, date_string), 'a+')
	ll_logfile.write('%s  %.3f %.3f\n'%(file_time.strftime('%Y%m%d-%H%M%S'), radar_lat, radar_lon))
	ll_logfile.close()


    else:
	pass


        # okay, now make an option to output files in a netcdf
    if 'output_nc' in ycfg.keys() and ycfg['output_nc']:

        if 'nc_path' in ycfg.keys():
            nc_path = ycfg['nc_path']
        else:
            nc_path = '.'
        
        nc_outfile = '%s/%s_%s_processed_gridded.nc'%(nc_path, ycfg['nc_prefix'], dt_string)

        print 'outputting processed data to a netcdf: {}'.format(nc_outfile)

        # create the output dictionary, then loop thru the fields and attach them
        out = {
                 'lat': {'data': lat, 'dims': ('y', 'x'), 'type': 'f'},
                 'lon': {'data': lon, 'dims': ('y', 'x'), 'type': 'f'},
                 'z-heights': {'data': z, 'dims': ('z'), 'type': 'f'},
		'rr': {'data': rr_array, 'dims': ('y', 'x'), 'type': 'f'},
        }

        for ncf in ycfg['nc_fields']:
            # first, let's check if ncf is actually in the radar data
            if ncf in radar.data.keys():
                out[ncf] = {'data': radar.data[ncf], 'type': 'f'}
            # next need to figured out the dimensions. If it's 2D, it's gonna be x and y
            # If it's 3D, x, y, z
                if radar.data[ncf].ndim == 2:
                    out[ncf]['dims'] = ('y', 'x')
                elif radar.data[ncf].ndim == 3:
                    out[ncf]['dims'] = ('z', 'y', 'x')
                else:
                    print 'check on the dimensions of the {} variable, does not seem to be 2D or 3D'.format(ncf)



        # this will remain constant no matter what
        dim_dict = {'x': rr_array.shape[1], 'y': rr_array.shape[0], 'z': len(z)}

        attr_dict = {'validtime': file_dt_string,
                        'created': 'Brody Fuchs, CSU, brfuchs@atmos.colostate.edu'}

        atools.write_dict_to_netcdf(nc_outfile, out, dim_dict, attr_dict=attr_dict)


    else:
        pass



    #######################################################33


    # ***************** Going to delete the radar object here to try to save memory
    #del radar
    gc.collect() # Doing some garbage collecting here to try to save some memory as well.


print 'Made it to the end of the processing'



