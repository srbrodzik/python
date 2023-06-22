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


km2nm = 0.54 # 1 km is 0.54 nautical miles
nm2km = 1.0/km2nm # Just to have both conversions

print 'start'

# files coming in from inpath and the updated files will be output to outpath
base_path = os.path.dirname(os.path.realpath('__file__'))


rc_flag = True
# This parses command line arguments, which is how a lot of this will be done

parser = argparse.ArgumentParser(description='Put in a file to be processed')

parser.add_argument('--file', action="store", dest="file")
parser.add_argument('--realtime', action="store", dest="realtime", type=bool, default=False)
parser.add_argument('--config', action="store", dest="config", default=None)

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

#print 'config data: {}'.format(ycfg)

outpath = '%s/%s'%(base_path, ycfg['path2'])
#inpath = '%s/radar_data/seapol'%(base_path)
prefix = ycfg['prefix2']
DPI = 120
largeDPI = 200


# this is how wide the plot will be
lat_width = ycfg['lat_width'] # degrees
lon_width = ycfg['lon_width'] # degrees


radarfile = deepcopy(rb)

rbase = os.path.basename(radarfile)


a_index = rbase.rfind('A')

dt_string = rbase[a_index+1:a_index+16]



# if '.' in rbase:
# 	last_dot = rbase.rfind('.')
# 	dt_string = rbase[a_index+1:last_dot]
# else:
# 	last_dot = None
# 	dt_string = rbase[a_index+1:]

# last_dot = rf.rfind('.')
# d_index = rf.rfind('d')

# dt_string = rf[d_index+2:d_index+17]
# print dt_string

nsweeps_thresh = ycfg['nsweeps_thresh']



#print dt_string
file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')

date_string = file_time.strftime('%Y%m%d')

#major_circle_rads = np.arange(20, 140, 20)
major_circle_rads = np.arange(20, 100, 20)
minor_circle_rads = np.arange(10, 110, 20)
azimuths = np.arange(0, 360, 30)
#minor_azimuths = np.arange(22.5, 360, 45)
minor_azimuths = None


vmax = 20.0
rhi_y_max = 15.0
rhi_x_max = 120


# Read in the UF file with pyart here, we keep the file field names for simplicity: these are the 2
# character names mandated by UF convention


radar = pyart.io.cfradial.read_cfradial(radarfile, file_field_names=True)
#cfg = RadarConfig.RadarConfig(dz='DBZ2', zdr='ZDR2', kdp='KDP2', rho='RHOHV2', hid='HID')


cfg = RadarConfig.RadarConfig(dz=ycfg['dz_name'], zdr=ycfg['zdr_name'], kdp=ycfg['kdp_name'], 
                        rho=ycfg['rho_name'], hid='HID', ph=ycfg['phase_name'], vel=ycfg['vel_name'])

# This checks for large jump in azimuthal data and blocks them out so pyart doesn't try to interpolate them



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


display = pyart.graph.RadarDisplay(radar)

if radar.scan_type == 'ppi':

    limx = [-120, 120]
    limy = [-120, 120]
    sweeps = [0]
    display_func = display.plot_ppi
    #cb_orient = 'vertical'

elif radar.scan_type == 'rhi':
    limx = [0, rhi_x_max]
    limy = [0, rhi_y_max]
    sweeps = np.arange(radar.nsweeps)
    display_func = display.plot


for i, ns in enumerate(sweeps):
    try:
    #print 'i: {}, ns: {}'.format(i, ns)
        fig = plt.figure(figsize=(8, 10))

        ax0 = fig.add_subplot(321)
        display_func(cfg.dz_name, sweep=ns, vmin=cfg.plot_params[cfg.dz_name]['lims'][0], vmax=cfg.plot_params[cfg.dz_name]['lims'][1], 
		cmap=cfg.plot_params[cfg.dz_name]['cmap'], colorbar_label='Reflectivity (dBZ)', mask_outside=True)
        display.set_limits(xlim=limx, ylim=limy)

        ax2 = fig.add_subplot(323)
        display_func(cfg.zdr_name, sweep=ns,  vmin=cfg.plot_params[cfg.zdr_name]['lims'][0], vmax=cfg.plot_params[cfg.zdr_name]['lims'][1], 
	       cmap=cfg.plot_params[cfg.zdr_name]['cmap'], colorbar_label='%s %s'%(cfg.plot_params[cfg.zdr_name]['name'], 
	    cfg.plot_params[cfg.zdr_name]['units']), mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)

        ax1 = fig.add_subplot(322)

        display_func(vel_name, sweep=ns, vmin=-1*vmax, vmax=vmax, cmap=ctables.Carbone11, 
	     colorbar_label='Radial velocity (m/s)', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)


        ax3 = fig.add_subplot(324)

        display_func(cfg.rho_name, sweep=ns, vmin=0.92, vmax=cfg.plot_params[cfg.rho_name]['lims'][1], 
	    cmap=cfg.plot_params[cfg.rho_name]['cmap'], colorbar_label='%s %s'%(cfg.plot_params[cfg.rho_name]['name'], 
		    cfg.plot_params[cfg.rho_name]['units']), 
	    mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)


        ax4 = fig.add_subplot(325)


        good_phase = radar.fields[cfg.ph_name]['data'][radar.fields[cfg.ph_name]['data']>= -360.0]

    #low_phase = atools.round_to_value(np.percentile(good_phase, 10), 10)
    #high_phase = atools.round_to_value(np.percentile(good_phase, 90), 10)

    #print low_phase, high_phase
        med_phase = atools.round_to_value(np.percentile(good_phase, 50), 10)
        print med_phase

        low_phase = 30.0
        high_phase = 120.0


        display_func(cfg.ph_name, sweep=ns, vmin=low_phase, vmax=high_phase, cmap=cfg.plot_params[cfg.ph_name]['cmap'],
	     colorbar_label='Differential phase ($^{\circ}$)', mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)



        ax5 = fig.add_subplot(326)
    # display.plot_ppi(cfg.kdp_name, sweep=sweep, vmin=0, vmax=1, cmap=cfg.cmaps[cfg.kdp_name],
    #                  colorbar_label='%s %s'%(cfg.names[cfg.kdp_name], cfg.units[cfg.kdp_name]), mask_outside=True)

        display_func(cfg.kdp_name, sweep=ns, vmin=-1.0, vmax=6.0, cmap=cfg.plot_params[cfg.kdp_name]['cmap'],
		 colorbar_label='%s %s'%(cfg.plot_params[cfg.kdp_name]['name'], cfg.plot_params[cfg.kdp_name]['units']), mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)

        all_axes = [ax0, ax1, ax2, ax3, ax4, ax5]


        if radar.scan_type == 'ppi':

	    for a in all_axes:
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

	    angle_title = np.average(radar.get_elevation(ns))



        elif radar.scan_type == 'rhi':
	    for a in all_axes:
	       a.set_title('')


	    angle_title = np.average(radar.get_azimuth(ns))

        plt.tight_layout()

        fig.suptitle('SEAPOL %s %s %.1f$^{\circ}$'%(file_dt_string, radar.scan_type.upper(), angle_title), fontsize=16)
        fig.subplots_adjust(top=0.93)

        print 'saving %s figure %.1f'%(radar.scan_type.upper(), angle_title)

        ptools.save_figure(fig, '%s/figures/raw_%s/%s'%(base_path, radar.scan_type, date_string), 
                '%s%s_%s_sweep%d.png'%(prefix, dt_string, radar.scan_type, i), dpi=DPI)

        plt.close(fig)

    except Exception, fe:
        print 'Error with the figure making: {}'.format(fe)


# Plotting a test figure with the good data flag and the signal to noise ratio
# testing DROPS output

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

    fig.suptitle('SEAPOL %s %s %.1f$^{\circ}$'%(file_dt_string, radar.scan_type.upper(), angle_title), fontsize=13)
    fig.subplots_adjust(top=0.90)

    print 'saving large %s figure'%(radar.scan_type.upper())

    #plt.savefig('%s/figures/large_dbz_vel_%s/%s%s_%s.png'%(base_path, radar.scan_type, prefix, dt_string, radar.scan_type), dpi=largeDPI)
    ptools.save_figure(fig, '%s/figures/large_dbz_vel_%s/%s'%(base_path, radar.scan_type, date_string), 
                    '%s%s_%s.png'%(prefix, dt_string, radar.scan_type), dpi=largeDPI)

    plt.close(fig)


# except Exception, lfe:
#     print 'Error with making large figure: {}'.format(lfe)

# # now writing the updated file to a cfradial file to the outpath. Trying to keep the same datetime stamps


if radar.scan_type == 'rhi':
    pass
else:
    sys.exit()

# Now if it's an RHI, do an elongated plot
for i, ns in enumerate(sweeps):
    #try:
    if True:

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


        blah = display_func(cfg.ph_name, sweep=ns, vmin=low_phase, vmax=high_phase, cmap=cfg.plot_params[cfg.ph_name]['cmap'],
                 colorbar_label='Differential phase ($^{\circ}$)', mask_outside=False, colorbar_flag=True)
        display.set_limits(xlim=limx, ylim=limy)


        ax5 = fig.add_subplot(616)
# display.plot_ppi(cfg.kdp_name, sweep=sweep, vmin=0, vmax=1, cmap=cfg.cmaps[cfg.kdp_name],
#                  colorbar_label='%s %s'%(cfg.names[cfg.kdp_name], cfg.units[cfg.kdp_name]), mask_outside=True)

        display_func(cfg.kdp_name, sweep=ns, vmin=-1.25, vmax=6.0, cmap=cfg.plot_params[cfg.kdp_name]['cmap'],
             colorbar_label='%s %s'%(cfg.plot_params[cfg.kdp_name]['name'], cfg.plot_params[cfg.kdp_name]['units']), mask_outside=False)
        display.set_limits(xlim=limx, ylim=limy)

        all_axes = [ax0, ax1, ax2, ax3, ax4, ax5]

        for a in all_axes:
            a.set_title('')
            a.set_xlabel('')
            a.set_ylabel('')
            a.set_yticks(np.arange(0, rhi_y_max, 5))
            a.set_yticklabels(np.arange(0, rhi_y_max, 5))
            a.grid(True)


        angle_title = np.average(radar.get_azimuth(ns))

        plt.tight_layout()

        fig.suptitle('SEAPOL %s %s %.1f$^{\circ}$'%(file_dt_string, radar.scan_type.upper(), angle_title), fontsize=16)
        fig.subplots_adjust(top=0.93, right=0.95)

        print 'saving %s true aspect figure %.1f'%(radar.scan_type.upper(), angle_title)

#        plt.savefig('%s/figures/rhi_true_aspect/%s%s_%s_sweep%d.png'%(base_path, prefix, dt_string, 
#                                radar.scan_type, i), dpi=DPI)
	ptools.save_figure(fig, '%s/figures/rhi_true_aspect/%s'%(base_path, date_string), 
                    '%s%s_%s_sweep%d.png'%(prefix, dt_string, radar.scan_type, i), dpi=DPI)

    # except Exception, te:
    #     print 'error with true aspect RHI: {}'.format(te)
    #     pass


    # now try to do self consistency with RHIs
    


####################################################################################
# ************* Down here at the end, get echo top heights and log them
#####################################################################################

PRF = 1200.0
mur = 3.0e8/(2.0*PRF)

dbz_levs = [0, 10, 20, 30, 40, 50]

#hts = OrderedDict()
hts = OrderedDict()

for dl in dbz_levs:
    hts[dl] = 0.0



gate_spacing = mur/radar.ngates

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
logfile.write('%s  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f\n'%(file_time.strftime('%Y%m%d-%H%M%S'), hts[0], hts[10], hts[20], hts[30], hts[40], hts[50]))
logfile.close()






