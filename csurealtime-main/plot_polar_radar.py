# Brody Fuchs, Sept 2017, brfuchs@atmos.colostate.edu

# code to read in the raw uf files and use PyART
# to do some QC'ing in real time

import pyart

from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, 
                            csu_dsd, csu_kdp, csu_misc, fundamentals)

import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import argparse
import os
import glob
from copy import deepcopy
from csuram import RadarConfig
import datetime
import Config
import ctables
import sys
from csuram import analysis_tools as atools
import gentools as gtools 
import plot_tools as ptools 
import radar_tools as radtools 
from collections import OrderedDict
from matplotlib.colors import LogNorm
from skewPy import SkewT
import traceback


km2nm = 0.54 # 1 km is 0.54 nautical miles
nm2km = 1.0/km2nm # Just to have both conversions

print 'start'

# files coming in from inpath and the updated files will be output to outpath
try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception:
    base_path = os.path.dirname(os.path.realpath('__file__'))





parser = argparse.ArgumentParser(description='Put in a file to be processed')

parser.add_argument('--file', action="store", dest="file")
parser.add_argument('--realtime', action="store", dest="realtime", type=int, default=0)
parser.add_argument('--config', action="store", dest="config", default=None)
parser.add_argument('--analysis', action="store", dest="analysis", type=int, default=0)

pargs = parser.parse_args()

if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v

if pargs.realtime:
    print 'QC Radar in realtime'
    pass
elif pargs.file is not None:
    rb = pargs.file
else:
    rb = 'SEA20170901_000018.uf'

prefix = ycfg['qc_prefix']

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




fig_path = ycfg['fig_path']


# this is how wide the plot will be
lat_width = ycfg['lat_width'] # degrees
lon_width = ycfg['lon_width'] # degrees


radarfile = deepcopy(rb)

rbase = os.path.basename(radarfile)


a_index = rbase.rfind(prefix[-1])

dt_string = rbase[a_index+1:a_index+16]


nsweeps_thresh = ycfg['nsweeps_thresh']



#print dt_string
file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')

date_string = file_time.strftime('%Y%m%d')

#major_circle_rads = np.arange(20, 140, 20)
major_circle_rads = np.arange(20, 100, 20)
minor_circle_rads = np.arange(10, 110, 20)
azimuths = np.arange(0, 360, 30)
rhi_azimuths = np.arange(0, 60, 10)
#minor_azimuths = np.arange(22.5, 360, 45)
minor_azimuths = None


rhi_y_max = 18.0
rhi_x_max = 120


radar = pyart.io.cfradial.read_cfradial(radarfile, file_field_names=True)

cfg = RadarConfig.RadarConfig(dz=ycfg['dz_name'], zdr=ycfg['zdr_name'], kdp=ycfg['kdp_name'], 
                        rho=ycfg['rho_name'], hid='HID', ph=ycfg['phase_name'], vel=ycfg['vel_name'])

radar_lat = radar.latitude['data'][0]
radar_lon = radar.longitude['data'][0]

# This checks for large jump in azimuthal data and blocks them out so pyart doesn't try to interpolate them

# Try to get the nyquist from the radar file, but if not, check to see if vmax is specified in the file
# otherwise just use a default value

try:
	# round it to the nearest 0.5 m/s
    vmax = 0.5*np.round(radar.get_nyquist_vel(0)/0.5)
except:
    if 'vmax' in ycfg.keys():
	vmax = ycfg['vmax']
    else:
	vmax = 20.0

print 'vmax used for plotting: {}'.format(vmax)

PRF = 1000.0
mur = 3.0e8/(2.0*PRF)
gate_spacing = mur/radar.ngates


if radar.scan_type == 'ppi':

    az_diff = np.diff(radar.azimuth['data'])
    jumps = np.where(az_diff >= 30.0)[0]

    if len(jumps):
        for f in radar.fields.keys():
	    try:
		for j in jumps:

		#print 'masking {}'.format(j)
		#radar.fields[f]['data'][j-1].mask = True
		    radar.fields[f]['data'][j].mask = True
		    radar.fields[f]['data'][j+1].mask = True
	    except Exception:
		pass


print 'radar fields: {}'.format(radar.fields.keys())

phase_name = ycfg['phase_name']
vel_name = ycfg['vel_name']
sq_name = ycfg['sq_name']

rscantype = radar.scan_type
nsweeps = radar.nsweeps

# This is a check of whether or not to go on...
valid_check = ( (rscantype == 'ppi') and (nsweeps >= nsweeps_thresh) ) or ( (rscantype == 'rhi') ) 


print 'Radar scan type: {}, number of sweeps: {}, valid: {}'.format(rscantype, nsweeps, valid_check)

print 'Max reflectivity: {}'.format(radar.fields[cfg.dz_name]['data'].max())



radar_base = os.path.basename(radarfile)
dot_loc = radar_base.rfind('.')
new_radar_name = '%s%s_%s.nc'%(prefix, dt_string, radar.scan_type)
#prefix, dt_string, radar.scan_type


# Wanna think about how to plot over the top of a map...


if radar.scan_type == 'ppi':

    if 'map_flag' in ycfg.keys() and ycfg['map_flag']:
        # this means we want to plot on a map
        display = pyart.graph.RadarMapDisplay(radar)
        display_func = display.plot_ppi_map
        limx = None
        limy = None
        extra_args = {'min_lon': radar_lon-lon_width/2.0, 'max_lon': radar_lon+lon_width/2.0, 
                        'min_lat': radar_lat-lat_width/2.0, 'max_lat': radar_lat+lat_width/2.0,
                     'projection': 'lcc', 'resolution': 'i', 'lat_0': radar_lat, 'lon_0': radar_lon}
        print 'PPI on a map'

    else:
        limx = [-120, 120]
        limy = [-120, 120]
        display = pyart.graph.RadarDisplay(radar)
        display_func = display.plot_ppi
        extra_args = {}
        print 'PPI not on a map'

    # regardless of map or not, going to only plot the first PPI sweep, at least in real time
    if 'ppi_elev_plot' in ycfg.keys():
        sweeps = [ycfg['ppi_elev_plot']]
    else:
        sweeps = [0]


elif radar.scan_type == 'rhi':
    limx = [0, rhi_x_max]
    limy = [0, rhi_y_max]
    sweeps = np.arange(radar.nsweeps)
    display = pyart.graph.RadarDisplay(radar)
    display_func = display.plot
    print 'RHI plot'
    extra_args = {}


for i, ns in enumerate(sweeps):
    try:
    #print 'i: {}, ns: {}'.format(i, ns)
        fig = plt.figure(figsize=(8, 10))

        ax0 = fig.add_subplot(321)
        display_func(cfg.dz_name, sweep=ns, vmin=cfg.plot_params[cfg.dz_name]['lims'][0], vmax=cfg.plot_params[cfg.dz_name]['lims'][1], 
		cmap=cfg.plot_params[cfg.dz_name]['cmap'], colorbar_label='Reflectivity (dBZ)', mask_outside=True, **extra_args)
        display.set_limits(xlim=limx, ylim=limy)


        ax1 = fig.add_subplot(322)
        display_func(vel_name, sweep=ns, vmin=-1*vmax, vmax=vmax, cmap=ctables.Carbone11, 
             colorbar_label='Radial velocity (m/s)', mask_outside=False, **extra_args)
        display.set_limits(xlim=limx, ylim=limy)


        ax2 = fig.add_subplot(323)
        display_func(cfg.zdr_name, sweep=ns,  vmin=cfg.plot_params[cfg.zdr_name]['lims'][0], vmax=cfg.plot_params[cfg.zdr_name]['lims'][1], 
	       cmap=cfg.plot_params[cfg.zdr_name]['cmap'], colorbar_label='%s %s'%(cfg.plot_params[cfg.zdr_name]['name'], 
	    cfg.plot_params[cfg.zdr_name]['units']), mask_outside=False, **extra_args)
        display.set_limits(xlim=limx, ylim=limy)


        ax3 = fig.add_subplot(324)
        display_func(cfg.rho_name, sweep=ns, vmin=0.92, vmax=cfg.plot_params[cfg.rho_name]['lims'][1], 
	    cmap=cfg.plot_params[cfg.rho_name]['cmap'], colorbar_label='%s %s'%(cfg.plot_params[cfg.rho_name]['name'], 
		    cfg.plot_params[cfg.rho_name]['units']), mask_outside=False, **extra_args)
        display.set_limits(xlim=limx, ylim=limy)


        ax4 = fig.add_subplot(325)

        low_phase = ycfg['low_phase']
        high_phase = ycfg['high_phase']

        display_func(cfg.ph_name, sweep=ns, vmin=low_phase, vmax=high_phase, cmap=cfg.plot_params[cfg.ph_name]['cmap'],
	     colorbar_label='Differential phase ($^{\circ}$)', mask_outside=False, **extra_args)
        display.set_limits(xlim=limx, ylim=limy)



        ax5 = fig.add_subplot(326)

        display_func(cfg.kdp_name, sweep=ns, vmin=-1.0, vmax=6.0, cmap=cfg.plot_params[cfg.kdp_name]['cmap'],
		 colorbar_label='%s %s'%(cfg.plot_params[cfg.kdp_name]['name'], cfg.plot_params[cfg.kdp_name]['units']), 
                mask_outside=False, **extra_args)
        display.set_limits(xlim=limx, ylim=limy)




        all_axes = [ax0, ax1, ax2, ax3, ax4, ax5]


        if radar.scan_type == 'ppi':

	    for a in all_axes:
	       	a.set_xlabel('')
	       	a.set_ylabel('')
	       	a.set_title('')
	    #a.axis('off')
	       	a.tick_params(
		  axis='x',          # changes apply to the x-axis
		  which='both',      # both major and minor ticks are affected
		  bottom='off',      # ticks along the bottom edge are off
		  top='off',         # ticks along the top edge are off
		  labelbottom='off') # labels along the bottom edge are off

	       	a.tick_params(
		  axis='y',          # changes apply to the x-axis
		  which='both',      # both major and minor ticks are affected
		  left='off',      # ticks along the bottom edge are off
		  right='off',         # ticks along the top edge are off
		  labelleft='off') # labels along the bottom edge are off


	       	ptools.plot_range_rings(a, 0.0, 0.0, major_circle_rads, conv=nm2km, text_flag=False)
	       	ptools.plot_azimuths(a, 0.0, 0.0, azimuths, conv=nm2km, max_range=120.0, minor_azs=minor_azimuths)

	       	a.set_aspect('equal')
	       	a.scatter([0], [0], c='black', s=25, alpha=0.5, edgecolors='none')

		if 'map_flag' in ycfg.keys() and ycfg['map_flag']:
		    display.basemap.drawstates()  

	    angle_title = np.average(radar.get_elevation(ns))



        elif radar.scan_type == 'rhi':
	    for a in all_axes:
	       a.set_title('')


	    angle_title = np.average(radar.get_azimuth(ns))

        plt.tight_layout()

        fig.suptitle('%s %s %s %.1f$^{\circ}$'%(radar_title, file_dt_string, radar.scan_type.upper(), angle_title), fontsize=16)
        fig.subplots_adjust(top=0.93)

        print 'saving %s figure %.1f'%(radar.scan_type.upper(), angle_title)

        ptools.save_figure(fig, '%s/raw_%s/%s'%(fig_path, radar.scan_type, date_string), 
                '%s%s_%s_sweep%d.png'%(prefix, dt_string, radar.scan_type, i), dpi=DPI)

        plt.close(fig)

    except Exception, fe:
        print 'Error with the figure making: {}'.format(fe)
	traceback.print_exc()


# Plotting a test figure with the good data flag and the signal to noise ratio
# testing DROPS output
if False:

    fig = plt.figure(figsize=(8, 10))

    ax0 = fig.add_subplot(211)
    display_func('FL', sweep=0, vmin=0, vmax=1, 
        cmap=cfg.plot_params[cfg.ph_name]['cmap'], colorbar_label='Data flag', mask_outside=False)
    display.set_limits(xlim=limx, ylim=limy)

    ax2 = fig.add_subplot(212)
    display_func('SR', sweep=0,  vmin=0, vmax=130, cmap=cfg.plot_params[cfg.ph_name]['cmap'], colorbar_label='SNR', mask_outside=False)
    display.set_limits(xlim=limx, ylim=limy)

    plt.tight_layout()


    plt.savefig('testdropsfig.png')


# *****************************************************************
# *************** now make a larger figure with only 2 plots on it
# *****************************************************************


# try:
if True:

    try:

        lsweep = 0

    #display2 = pyart.graph.RadarDisplay(radar)

        fig = plt.figure(figsize=(11, 5))
        axs = []


        axs.append(fig.add_subplot(1, 2, 1))

        display_func(cfg.dz_name, sweep=lsweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0], vmax=cfg.plot_params[cfg.dz_name]['lims'][1], 
                cmap=cfg.plot_params[cfg.dz_name]['cmap'], colorbar_label='Reflectivity (dBZ)', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)


        axs.append(fig.add_subplot(1, 2, 2))
        display_func(vel_name, sweep=lsweep, vmin=-1*vmax, vmax=vmax, cmap=ctables.Carbone11, 
                 colorbar_label='Radial velocity (m/s)', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)



        if radar.scan_type == 'ppi':

            for a in axs:
        # a.set_xlabel('W-E distance (km)')
        # a.set_ylabel('N-S distance (km)')
                a.set_xlabel('')
                a.set_ylabel('')
                a.set_title('')
        #a.axis('off')
                a.tick_params(
                        axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom='off',      # ticks along the bottom edge are off
                        top='off',         # ticks along the top edge are off
                        labelbottom='off') # labels along the bottom edge are off

                a.tick_params(
                        axis='y',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        left='off',      # ticks along the bottom edge are off
                        right='off',         # ticks along the top edge are off
                        labelleft='off') # labels along the bottom edge are off




                ptools.plot_range_rings(a, 0.0, 0.0, major_circle_rads, conv=nm2km, text_flag=False)
                ptools.plot_azimuths(a, 0.0, 0.0, azimuths, conv=nm2km, max_range=120.0, minor_azs=minor_azimuths)

                a.set_aspect('equal')
                a.scatter([0], [0], c='black', s=25, alpha=0.5, edgecolors='none')

            angle_title = np.average(radar.get_elevation(0))



        elif radar.scan_type == 'rhi':
            for a in axs:
                a.set_title('')


            angle_title = np.average(radar.get_azimuth(0))

        plt.tight_layout()

        fig.suptitle('%s %s %s %.1f$^{\circ}$'%(radar_title, file_dt_string, radar.scan_type.upper(), angle_title), fontsize=13)
        fig.subplots_adjust(top=0.90)

        print 'saving large %s figure'%(radar.scan_type.upper())

    #plt.savefig('%s/figures/large_dbz_vel_%s/%s%s_%s.png'%(base_path, radar.scan_type, prefix, dt_string, radar.scan_type), dpi=largeDPI)
        ptools.save_figure(fig, '%s/large_dbz_vel_%s/%s'%(fig_path, radar.scan_type, date_string), 
                    '%s%s_%s.png'%(prefix, dt_string, radar.scan_type), dpi=largeDPI)

        plt.close(fig)


    except Exception, lfe:
        print 'Error with making large figure: {}'.format(lfe)
	traceback.print_exc()

# # now writing the updated file to a cfradial file to the outpath. Trying to keep the same datetime stamps


if radar.scan_type == 'rhi':
    pass
else:
    sys.exit()

# HID should already be in the radar object/file
# but if not, then go ahead and calculate it here
if 'HID' not in radar.fields.keys():

    snd_data = np.genfromtxt('%s/sample_tropical_sounding.txt'%(base_path), dtype=None, skip_header=6)
    snd_ht = snd_data['f1'].astype(float)
    snd_temp = snd_data['f2']

    radar_T, radar_z = radtools.interpolate_sounding_to_radar(snd_temp, snd_ht, radar)



    scores = csu_fhc.csu_fhc_summer(dz=radar.fields[cfg.dz_name]['data'], zdr=radar.fields[cfg.zdr_name]['data'], 
		rho=radar.fields[cfg.rho_name]['data'], kdp=radar.fields[cfg.kdp_name]['data'], use_temp=True, band='C', T=radar_T)
    fh = np.argmax(scores, axis=0) + 1


    radar = radtools.add_field_to_radar_object(fh, radar, field_name='HID', units='', 
				   long_name='Dominant hydrometeor classification',
				   standard_name='Hydrometeor type', 
				   dz_field=cfg.dz_name)




    # On the left, plot the HID confidence
    max_score = np.max(scores, axis=0)
    tot_score = np.sum(scores, axis=0)
    hid_conf = 100.0*max_score/tot_score
    hid_conf = np.ma.masked_invalid(hid_conf)

    # Now need to add it to the radar object

    radar = radtools.add_field_to_radar_object(hid_conf, radar, field_name='HIDCONF', units='%', 
				   long_name='HID confidence',
				   standard_name='HID confidence', 
                               dz_field=cfg.dz_name)

else:
    print 'HID already in stored in the radar object/file'
    pass


if False:
    print 'Plotting RHI HID'
    fig = plt.figure(figsize=(11,6))


    ax1 = fig.add_subplot(121)

    display_func('HIDCONF', sweep=this_sweep, cmap=plt.cm.hot_r, vmin=0, vmax=100, 
             colorbar_label='HID confidence (%)', mask_outside=False)
    display.set_limits(xlim=limx, ylim=limy)

    #cb = cfg.HID_barplot_colorbar(fig)
    ax1.set_title('')


    ax2 = fig.add_subplot(122)

    display_func('HID', sweep=this_sweep, cmap=cfg.plot_params['HID']['cmap'], norm=cfg.plot_params['HID']['norm'], 
             colorbar_label='HID', mask_outside=False, colorbar_flag=False)
    display.set_limits(xlim=limx, ylim=limy)

    cb = cfg.HID_barplot_colorbar(fig)
    ax2.set_title('')

    #ptools.plot_azimuths(ax1, 0.0, 0.0, rhi_azimuths, conv=nm2km, max_range=12000.0)

    plt.tight_layout()
    fig.subplots_adjust(top=0.92, right=0.89)

    fig.suptitle('%s: %s RHI HID'%(radar_title, file_dt_string), fontsize=16)

    ptools.save_figure(fig, '%s/rhi_HID/%s'%(fig_path, date_string), 
        '%s%s_%s_HID.png'%(prefix, dt_string, radar.scan_type), dpi=DPI)





# Now if it's an RHI, do an elongated plot
for i, ns in enumerate(sweeps):
    try:
    

        fig = plt.figure(figsize=(7, 10))


        ax0 = fig.add_subplot(611)
        display_func(cfg.dz_name, sweep=ns, vmin=cfg.plot_params[cfg.dz_name]['lims'][0], vmax=cfg.plot_params[cfg.dz_name]['lims'][1], 
                cmap=cfg.plot_params[cfg.dz_name]['cmap'], 
                 colorbar_label='Reflectivity (dBZ)', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)

        ax1 = fig.add_subplot(612)

        display_func(vel_name, sweep=ns, vmin=-1*vmax, vmax=vmax, cmap=ctables.Carbone11, 
                 colorbar_label='Radial velocity (m/s)', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)


        ax2 = fig.add_subplot(613)
        display_func(cfg.zdr_name, sweep=ns,  vmin=cfg.plot_params[cfg.zdr_name]['lims'][0], vmax=cfg.plot_params[cfg.zdr_name]['lims'][1], 
            cmap=cfg.plot_params[cfg.zdr_name]['cmap'],
                 colorbar_label='%s %s'%(cfg.plot_params[cfg.zdr_name]['name'], cfg.plot_params[cfg.zdr_name]['units']), mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)


        ax3 = fig.add_subplot(614)

        display_func(cfg.rho_name, sweep=ns, vmin=0.92, vmax=cfg.plot_params[cfg.rho_name]['lims'][1], cmap=cfg.plot_params[cfg.rho_name]['cmap'],
                 colorbar_label='%s %s'%(cfg.plot_params[cfg.rho_name]['name'], cfg.plot_params[cfg.rho_name]['units']), mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)


        ax4 = fig.add_subplot(615)


        # blah = display_func(cfg.ph_name, sweep=ns, vmin=low_phase, vmax=high_phase, cmap=cfg.plot_params[cfg.ph_name]['cmap'],
        #          colorbar_label='Differential phase ($^{\circ}$)', mask_outside=False, colorbar_flag=True)
        # display.set_limits(xlim=limx, ylim=limy)

        display_func(cfg.kdp_name, sweep=ns, vmin=-1.25, vmax=6.0, cmap=cfg.plot_params[cfg.kdp_name]['cmap'],
             colorbar_label='%s %s'%(cfg.plot_params[cfg.kdp_name]['name'], cfg.plot_params[cfg.kdp_name]['units']), mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)


        ax5 = fig.add_subplot(616)
# display.plot_ppi(cfg.kdp_name, sweep=sweep, vmin=0, vmax=1, cmap=cfg.cmaps[cfg.kdp_name],
#                  colorbar_label='%s %s'%(cfg.names[cfg.kdp_name], cfg.units[cfg.kdp_name]), mask_outside=True)

        display_func('HID', sweep=ns, cmap=cfg.plot_params['HID']['cmap'], norm=cfg.plot_params['HID']['norm'], 
                colorbar_label='HID', mask_outside=False, colorbar_flag=True)
        display.set_limits(xlim=limx, ylim=limy)

        hidcb_old = fig.axes[-1]
        #hidcb_old.set_xticks(np.arange(0, 11))
        hidcb_old.set_visible(False)


        hidcb = cfg.HID_barplot_colorbar(fig, location=[0.845, 0.027, 0.010, 0.125])
        hidcb.ax.tick_params(labelsize=9)





        all_axes = [ax0, ax1, ax2, ax3, ax4, ax5]

        for a in all_axes:
            a.set_title('')
            a.set_xlabel('')
            a.set_ylabel('')
            a.set_yticks(np.arange(0, rhi_y_max, 5))
            a.set_yticklabels(np.arange(0, rhi_y_max, 5))
            a.grid(True)


        angle_title = np.average(radar.get_azimuth(ns))

        #plt.tight_layout()

        fig.suptitle('%s %s %s %.1f$^{\circ}$'%(radar_title, file_dt_string, radar.scan_type.upper(), angle_title), fontsize=16)
        fig.subplots_adjust(top=0.95, right=0.98, bottom=0.02, left=0.08)

        print 'saving %s true aspect figure %.1f'%(radar.scan_type.upper(), angle_title)

#        plt.savefig('%s/figures/rhi_true_aspect/%s%s_%s_sweep%d.png'%(base_path, prefix, dt_string, 
#                                radar.scan_type, i), dpi=DPI)
	ptools.save_figure(fig, '%s/rhi_true_aspect/%s'%(fig_path, date_string), 
                    '%s%s_%s_sweep%d.png'%(prefix, dt_string, radar.scan_type, i), dpi=DPI)

    except Exception, te:
        print 'error with true aspect RHI: {}'.format(te)
        traceback.print_exc()


    # now try to do self consistency with RHIs
    



### Okay, now calculate dVr/dr and plot (at least to start with)
# Then also plot HID as well

vel_filled = radar.fields[vel_name]['data'].filled(0.0)
diverg = np.zeros_like(vel_filled)
diverg_mask = np.zeros(diverg.shape, bool)
smooth_width = 49
for ie in range(vel_filled.shape[0]):
    #mask_ne = [i for i, j in zip(radar.fields[vel_name]['data'].mask[1:], radar.fields[vel_name]['data'].mask[:-1]) if i != j]
    mask_ne = [i for i in range(radar.fields[vel_name]['data'].mask.shape[1]-1) if radar.fields[vel_name]['data'].mask[ie, i] != \
                                        radar.fields[vel_name]['data'].mask[ie, i+1]]
    #print mask_ne
    #quit()
    sm_vel = atools.convolve_smooth(range(len(vel_filled[ie])), vel_filled[ie], mode='extend', window='121', width=smooth_width)[0]
    diverg[ie] = np.gradient(sm_vel)/gate_spacing
    # Now need to loop thru each mask_ne and set the divergence values to 0 there (and within the smooth width!)

    for mn in mask_ne:
        minmask = int(np.max([0, mn-smooth_width/3.0]))
        maxmask = int(np.min([radar.fields[vel_name]['data'].mask.shape[1], mn+smooth_width/3.0]))
        diverg[ie][minmask:maxmask] = 0.0



# Okay, now let's add it to the radar (PyART) object so it's easier to plot
radar = radtools.add_field_to_radar_object(diverg, radar, field_name='DIV', units='', 
                               long_name='Range divergence (1/s)',
                               standard_name='Range divergence (1/s)', 
                               dz_field=cfg.dz_name)

print 'Plotting velocity and ray-based divergence'

# ******* Okay, now plotting the velocity and the ray-based divergence ********

this_sweep = 1
fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(121)

display_func(vel_name, sweep=this_sweep, vmin=-1*vmax, vmax=vmax, cmap=ctables.Carbone11, 
         colorbar_label='Radial velocity (m/s)', mask_outside=False)
display.set_limits(xlim=limx, ylim=limy)

ax2 = fig.add_subplot(122)

blah = display_func('DIV', sweep=this_sweep, vmin=-3e-3, vmax=3e-3, cmap=plt.cm.PiYG,
         colorbar_label='Range divergence (1/s)', mask_outside=False, colorbar_flag=True)
display.set_limits(xlim=limx, ylim=limy)

for a in [ax1, ax2]:
    a.set_title('')


fig.suptitle('%s %s %s %.1f$^{\circ}$'%(radar_title, file_dt_string, radar.scan_type.upper(), angle_title), fontsize=16)

plt.tight_layout()
fig.subplots_adjust(top=0.91)


ptools.save_figure(fig, '%s/rhi_velocity/%s'%(fig_path, date_string), 
    '%s%s_%s_vel_div_sweep%d.png'%(prefix, dt_string, radar.scan_type, this_sweep), dpi=DPI)

############################################################################################

# Now let's try self consistency

if 'radar_dsd_analysis' in ycfg.keys() and ycfg['radar_dsd_analysis']:



        #xyz = radar.get_gate_x_y_z(s)
        #beam_range = np.sqrt(xyz[0]**2 + xyz[1]**2)
        #zvals = xyz[2]
        #z_wh = zvals[wh_ht]

        print 'Plotting RHI self consistency'

        logmin = 5e-4
        logmax = 100

        dz_bins = np.arange(0, 70, 1.25)
        kdp_bins = np.arange(-1, 7, 0.25)
        zdr_bins = np.arange(-1, 6, 0.25)

        zdr_vals = []
        kdp_vals = []
        dz_vals = []

        zdr_vals_rain = []
        kdp_vals_rain = []
        dz_vals_rain = []
        rad_xyz = []


        for nsw in range(radar.nsweeps):
            # loop thru each sweep
            sw_xyz = radar.get_gate_x_y_z(nsw)
            rad_xyz.append(sw_xyz)
            sw_z = sw_xyz[2]
            dz_sw_data = radar.get_field(nsw, cfg.dz_name)
            zdr_sw_data = radar.get_field(nsw, cfg.zdr_name)
            kdp_sw_data = radar.get_field(nsw, cfg.kdp_name)

        #    kdp_rain = radar.data[radar.kdp_name][z_rain_bot_ind:z_rain_top_ind]
        #    zdr_rain = radar.data[radar.zdr_name][z_rain_bot_ind:z_rain_top_ind]
        #
        #    # kdp_rain = radar.data[radar.kdp_name][0:4]
        #    # zdr_rain = radar.data[radar.zdr_name][0:4]
        #
            dbz_mask = dz_sw_data.mask
            kdp_mask = kdp_sw_data.mask
            zdr_mask = zdr_sw_data.mask
            combo_mask = np.logical_or(np.logical_or(kdp_mask, zdr_mask), dbz_mask)

            dz_sw_data.mask = combo_mask
            kdp_sw_data.mask = combo_mask
            zdr_sw_data.mask = combo_mask


            dz_vals.extend(dz_sw_data.compressed().ravel())
            zdr_vals.extend(zdr_sw_data.compressed().ravel())
            kdp_vals.extend(kdp_sw_data.compressed().ravel())    

            in_rain = np.where(sw_z <= 4000.0)
            dz_vals_rain.extend(dz_sw_data[in_rain].compressed())   
            zdr_vals_rain.extend(zdr_sw_data[in_rain].compressed())   
            kdp_vals_rain.extend(kdp_sw_data[in_rain].compressed())   
         

        dz_vals = np.array(dz_vals)
        zdr_vals = np.array(zdr_vals)
        kdp_vals = np.array(kdp_vals)

        dz_vals_rain = np.array(dz_vals_rain)
        zdr_vals_rain = np.array(zdr_vals_rain)
        kdp_vals_rain = np.array(kdp_vals_rain)


        # ******* Now we're plotting the self consistency *************************

        fig, ax = plt.subplots(2, 2, figsize=(9, 9))
        axf = ax.flatten()

        dz_kdp_hist = np.histogramdd((dz_vals, kdp_vals),
        				bins=(dz_bins, kdp_bins), normed=True)[0]
        dz_kdp_hist *= 100.0/dz_kdp_hist.sum()
        dk_pc = axf[0].pcolormesh(dz_bins[:-1], kdp_bins[:-1], dz_kdp_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        dk_cb = plt.colorbar(dk_pc, ax=axf[0])

        #ax[0].scatter(radar.data[radar.dz_name].ravel(), radar.data[radar.kdp_name].ravel(), c='black', alpha=0.3, s=8, edgecolors='none')

        axf[0].set_xlabel('Z')
        axf[0].set_ylabel('Kdp')
        axf[0].set_title('All points Z vs Kdp')

        zdr_kdp_hist = np.histogramdd((zdr_vals, kdp_vals),
        				bins=(zdr_bins, kdp_bins), normed=True)[0]
        zdr_kdp_hist *= 100.0/zdr_kdp_hist.sum()
        zk_pc = axf[1].pcolormesh(zdr_bins[:-1], kdp_bins[:-1], zdr_kdp_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        zk_cb = plt.colorbar(zk_pc, ax=axf[1])

        # Second plot will be ZDR/Kdp
        axf[1].set_xlabel('ZDR')
        axf[1].set_ylabel('Kdp')
        axf[1].set_title('All points ZDR vs Kdp')


        dz_kdp_rain_hist = np.histogramdd((dz_vals_rain, kdp_vals_rain),
        				bins=(dz_bins, kdp_bins), normed=True)[0]
        dz_kdp_rain_hist *= 100.0/dz_kdp_rain_hist.sum()
        dkr_pc = axf[2].pcolormesh(dz_bins[:-1], kdp_bins[:-1], dz_kdp_rain_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        dkr_cb = plt.colorbar(dkr_pc, ax=axf[2])

        axf[2].set_xlabel('Z')
        axf[2].set_ylabel('Kdp')
        axf[2].set_title('Rain points Z vs Kdp')


        zdr_kdp_rain_hist = np.histogramdd((zdr_vals_rain, kdp_vals_rain),
        				bins=(zdr_bins, kdp_bins), normed=True)[0]
        zdr_kdp_rain_hist *= 100.0/zdr_kdp_rain_hist.sum()
        zkr_pc = axf[3].pcolormesh(zdr_bins[:-1], kdp_bins[:-1], zdr_kdp_rain_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        zkr_cb = plt.colorbar(zkr_pc, ax=axf[3])

        # Second plot will be ZDR/Kdp
        axf[3].set_xlabel('ZDR')
        axf[3].set_ylabel('Kdp')
        axf[3].set_title('Rain points ZDR vs Kdp')

        for a in axf:
            a.grid(True)


        fig.suptitle('%s: %s RHI self consistency'%(radar_title, file_dt_string), fontsize=16)

        plt.tight_layout()
        fig.subplots_adjust(top=0.92)

        ptools.save_figure(fig, '%s/rhi_self_consistency/%s'%(fig_path, date_string), 
            '%s%s_%s.png'%(prefix, dt_string, radar.scan_type), dpi=DPI)

        # *************** Now can we do a HID plot **********************
        # This means I need to load in the sounding....



        # ******************** Then do Nw/D0 retreivals **************************

        print 'Nw/D0 retreivals'

        d0, Nw, mu = csu_dsd.calc_dsd(dz=radar.fields[cfg.dz_name]['data'], zdr=radar.fields[cfg.zdr_name]['data'], \
        				kdp=radar.fields[cfg.kdp_name]['data'], band='C')

        Nw = np.log10(Nw)


        radar = radtools.add_field_to_radar_object(d0, radar, field_name='D0', units='', 
                                       long_name='Median drop diameter',
                                       standard_name='Median drop diameter', 
                                       dz_field=cfg.dz_name)

        radar = radtools.add_field_to_radar_object(Nw, radar, field_name='NW', units='', 
                                       long_name='Normalized intercept parameter',
                                       standard_name='Normalized intercept parameter', 
                                       dz_field=cfg.dz_name)


        #ocean_ncolors = 8
        #ocean_oldcmap = plt.cm.get_cmap("ocean_r", ocean_ncolors) #generate a jet map with 10 values
        #ocean_old_vals = ocean_oldcmap(np.arange(ocean_ncolors)) #extract those values as an array
        #ocean_old_vals[0] = [0.92, 0.92, 0.92, 1]
        #ocean_old_vals[1] = [0.75, 0.75, 0.75, 1]
        #ocean_cmap = mpl.colors.LinearSegmentedColormap.from_list("newspec", ocean_old_vals)

        new_stern_cmap = ptools.modify_cmap("gist_stern", {0: [0.3, 0.0, 0.0, 1.0], 1: [0.9, 0.0, 0.0, 1.0]}, 
        			ncolors=8, new_name='newcmap')


        fig = plt.figure(figsize=(14,5))


        ax1 = fig.add_subplot(131)

        display_func('D0', sweep=this_sweep, cmap=new_stern_cmap, vmin=0.7, vmax=2.5, 
                 colorbar_label='Median volume drop diameter (mm)', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)

        ax1.set_title('')

        ax2 = fig.add_subplot(132)


        display_func('NW', sweep=this_sweep, cmap=plt.cm.Spectral_r, vmin=2, vmax=6, 
                 colorbar_label='Normalized intercept parameter', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)

        ax2.set_title('')

        ax3 = fig.add_subplot(133)

        display_func(cfg.dz_name, sweep=this_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0], vmax=cfg.plot_params[cfg.dz_name]['lims'][1], 
        		cmap=cfg.plot_params[cfg.dz_name]['cmap'], colorbar_label='Reflectivity (dBZ)', mask_outside=True)
        display.set_limits(xlim=limx, ylim=limy)
        ax3.set_title('')

        plt.tight_layout()
        fig.subplots_adjust(top=0.92)

        fig.suptitle('%s: %s RHI DSD retreivals'%(radar_title, file_dt_string), fontsize=16)


        ptools.save_figure(fig, '%s/experimental/%s'%(fig_path, date_string), 
            '%s%s_%s_DSDparams.png'%(prefix, dt_string, radar.scan_type), dpi=DPI)

        rad_x_arr = rad_xyz[0][0]
        rad_y_arr = rad_xyz[0][1]
        rad_z_arr = rad_xyz[0][2]


        for rx in range(1, len(rad_xyz)):
            rad_x_arr = np.vstack((rad_x_arr, rad_xyz[rx][0]))
            rad_y_arr = np.vstack((rad_y_arr, rad_xyz[rx][1]))
            rad_z_arr = np.vstack((rad_z_arr, rad_xyz[rx][2]))

        # now grab the points that are below 4000 m in Z to start

        rain_z_ind = np.where(rad_z_arr <= 4000.0)

        d0_nw_combo_mask = np.logical_or(d0.mask, Nw.mask)
        d0.mask = d0_nw_combo_mask
        Nw.mask = d0_nw_combo_mask

        d0_rain = d0[rain_z_ind].compressed()
        nw_rain = Nw[rain_z_ind].compressed()

        d0_bins = np.linspace(0, 5, 25)
        nw_bins = np.linspace(2, 6, 25)


        nw_d0_rain_hist = np.histogramdd((d0_rain, nw_rain),
        				bins=(d0_bins, nw_bins), normed=True)[0]
        nw_d0_rain_hist *= 100.0/nw_d0_rain_hist.sum()


        fig, ax = plt.subplots(1, 1, figsize=(8,6))


        #ndr_pc = ax.pcolormesh(d0_bins[:-1], nw_bins[:-1], nw_d0_rain_hist.T, cmap=plt.cm.cubehelix_r, norm=LogNorm(vmin=logmin, vmax=logmax))
        #ndr_cb = plt.colorbar(ndr_pc, ax=ax)
        #ndr_cb.set_label('Frequency (%)')

        ax.scatter(d0_rain, nw_rain, color='0.5', alpha=0.2, edgecolors='none', s=6)


        ax.set_xlabel('D$_0$ (mm)')
        ax.set_ylabel('log N$_w$')
        ax.grid(True)
        ax.set_xlim(0, 4)
        ax.set_ylim(1, 6)



        plt.tight_layout()
        fig.subplots_adjust(top=0.92)

        fig.suptitle('%s: %s RHI DSD scatter'%(radar_title, file_dt_string), fontsize=16)

        ptools.save_figure(fig, '%s/experimental/%s'%(fig_path, date_string), 
            '%s%s_%s_nwd0scatter.png'%(prefix, dt_string, radar.scan_type), dpi=DPI)



####################################################################################
# ************* Down here at the end, get echo top heights and log them
#####################################################################################

if pargs.analysis:


    dbz_levs = [0, 10, 20, 30, 40, 50]

    #hts = OrderedDict()
    hts = OrderedDict()

    for dl in dbz_levs:
	hts[dl] = 0.0




    # we'll need to go thru each sweep and look for where values are above 30 dBZ




    for s in range(radar.nsweeps):
	dbz_data = radar.get_field(s, cfg.dz_name)

	for h in dbz_levs:
			    #print '{} dBZ'.format(h)
	    wh_ht = np.where(dbz_data >= float(h))
	    if len(wh_ht[0]):
		xyz = radar.get_gate_x_y_z(s)
		beam_range = np.sqrt(xyz[0]**2 + xyz[1]**2)
		zvals = xyz[2]
		z_wh = zvals[wh_ht]

		z_rep_val = np.percentile(z_wh/1000.0, 95)
			    #       print '{} dBZ, z value: {}, current max z val: {}'.format(h, z_rep_val, hts[h][-1])

		if z_rep_val > hts[h]:
		    hts[h] = z_rep_val
			    #               print 'reassigning the max z val: {}'.format(hts[h][-1])
					    #time.sleep(0.05)

		else:
			    #               print 'not reassigning'
		    pass



    logfile = open('%s/radar_stats/%s_dbzhts.log'%(base_path, date_string), 'a+')
    logfile.write('%s  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f\n'%(file_time.strftime('%Y%m%d-%H%M%S'), hts[0], hts[10], hts[20], 
					hts[30], hts[40], hts[50]))
    logfile.close()




