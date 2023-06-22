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
import traceback
import scipy.ndimage as ndi
import attenuation_corr as acorr


km2nm = 0.54 # 1 km is 0.54 nautical miles
nm2km = 1.0/km2nm # Just to have both conversions

# files coming in from inpath and the updated files will be output to outpath

try:
    base_path = os.path.dirname(os.path.realpath(__file__))
except Exception:
    base_path = os.path.dirname(os.path.realpath('__file__'))



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


outpath = ycfg['qc_data_path']
#inpath = '%s/radar_data/seapol'%(base_path)
prefix = ycfg['qc_prefix']


# this is how wide the plot will be
lat_width = ycfg['lat_width'] # degrees
lon_width = ycfg['lon_width'] # degrees


rc_flag = ycfg['rhi_range_correction']
print 'Range correction flag: {}'.format(rc_flag)

radarfile = deepcopy(rb)

rbase = os.path.basename(radarfile)


a_index = rbase.rfind(prefix[-1])

if '.' in rbase:
    last_dot = rbase.rfind('.')
    dt_string = rbase[a_index+1:last_dot]
else:
    last_dot = None
    dt_string = rbase[a_index+1:]

nsweeps_thresh = ycfg['nsweeps_thresh']

#print dt_string
file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')
date_string = file_time.strftime('%Y%m%d')



major_circle_rads = np.arange(20, 100, 20)
minor_circle_rads = np.arange(10, 110, 20)
azimuths = np.arange(0, 360, 30)
minor_azimuths = None


rhi_y_max = 12.0
rhi_x_max = 60


qc_flag = ycfg['qc_flag']
print 'QC flag: {}'.format(qc_flag)

unfold_flag = ycfg['unfold_flag']
print 'unfold flag: {}'.format(unfold_flag)


if 'DPI' in ycfg.keys():
    DPI = ycfg['DPI']
else:
    DPI = 120

if 'largeDPI' in ycfg.keys():
    largeDPI = ycfg['largeDPI']
else:
    largeDPI = 200



# Read in the UF file with pyart here, we keep the file field names for simplicity: these are the 2
# character names mandated by UF convention

if 'uf' in rbase:
    radar = pyart.io.read_uf(radarfile, file_field_names=True)
    cfg = RadarConfig.RadarConfig(dz='DZ', zdr='DR', kdp='KD', rho='RH', hid='HID')

else:
    radar = pyart.io.read(radarfile, file_field_names=True)
        #cfg = RadarConfig.RadarConfig(dz='DBZ2', zdr='ZDR2', kdp='KDP2', rho='RHOHV2', hid='HID')

    cfg = RadarConfig.RadarConfig(dz=ycfg['qc_dz_name'], zdr=ycfg['qc_zdr_name'], kdp=ycfg['qc_kdp_name'], 
                                rho=ycfg['qc_rho_name'], hid='HID', ph=ycfg['qc_phase_name'])


print 'scan type: {}'.format(radar.scan_type)

# This checks for large jump in azimuthal data and blocks them out so pyart doesn't try to interpolate them

az_diff = np.diff(radar.azimuth['data'])
jumps = np.where(np.abs(az_diff) >= 30.0)[0]

if len(jumps):
    for f in radar.fields.keys():
        for j in jumps:

            #print 'masking {}'.format(j)
            #radar.fields[f]['data'][j-1].mask = True
            radar.fields[f]['data'][j].mask = True
            radar.fields[f]['data'][j+1].mask = True


try:
    if len(jumps) > 0:
	az_ends = radar.azimuth['data'][jumps]
	az_starts = radar.azimuth['data'][jumps+1]

	logfile = open('%s/radar_stats/%s_sectorlims.log'%(base_path, date_string), 'a+')
	logfile.write('%s  %.1f %.1f\n'%(file_time.strftime('%Y%m%d-%H%M%S'), np.median(az_ends), np.median(az_starts)))
	logfile.close()

except Exception, aze:
    print 'Error with logging sector lims: {}'.format(aze)


# **************** Apply range correction if True **********************

if rc_flag and radar.scan_type == 'rhi': # this only needs to happen for the RHIs, which is just insane
    range_correction = 20*np.log10(radar.range['data']/1000.0) - 35.0
    print 'Original reflectivity max: {}'.format(radar.fields[ycfg['orig_dz_name']]['data'].max())
    radar.fields[ycfg['orig_dz_name']]['data'] += range_correction
    print 'Corrected reflectivity max: {}'.format(radar.fields[ycfg['orig_dz_name']]['data'].max())



print 'original radar.fields.keys: {}'.format(radar.fields.keys())


phase_name = ycfg['orig_phase_name']
vel_name = ycfg['orig_vel_name']
sq_name = ycfg['orig_sq_name']

rscantype = radar.scan_type
nsweeps = radar.nsweeps

# This is a check of whether or not to go on...
valid_check = ( (rscantype == 'ppi') and (nsweeps >= nsweeps_thresh) ) or ( (rscantype == 'rhi') ) 


print 'Radar scan type: {}, number of sweeps: {}, valid: {}'.format(rscantype, nsweeps, valid_check)


if valid_check:


    radar_lat = radar.latitude['data'][0]
    radar_lon = radar.longitude['data'][0]


    dz = radar.fields[ycfg['orig_dz_name']]['data']
    dr = radar.fields[ycfg['orig_zdr_name']]['data']
    kd = radar.fields[ycfg['orig_kdp_name']]['data']
    rh = radar.fields[ycfg['orig_rho_name']]['data']


    cc_flag = True
    cc_var = 'CC'

    try:
        ccorr = deepcopy(radar.fields['UNKNOWN_82']['data'])
        ccorr[ccorr<0] += 65535
        ccorr = (ccorr-1)/65533.0
    except Exception, cce:
        print 'Error with CCORR: {}'.format(cce)
	traceback.print_exc()
        cc_flag = False
        cc_var = ycfg['orig_rho_name']


    if cc_flag:
    # Use this total_mask on other fields as well, like rho and zdr to clean those up as well
        radar = radtools.add_field_to_radar_object(ccorr, radar, field_name='CC', units='', 
                                       long_name='Correlation coefficient',
                                       standard_name='Correlation coefficient', 
                                       dz_field=ycfg['orig_dz_name'])

    else:
        cc_var = 'RHOHV2'



    if qc_flag:
        


#             spec_at, cor_z = pyart.correct.calculate_attenuation(radar, 0, refl_field=ycfg['orig_dz_name'],
#                                                             ncp_field=ycfg['orig_sq_name'], rhv_field=ycfg['orig_rho_name'],
#                                                             phidp_field=ycfg['orig_phase_name'], rhv_min=0.7, ncp_min=0.4, 
# #                                                            a_coef=0.05, beta=0.9, fzl=4500.0)
#                                                             a_coef=0.045, beta=0.001, fzl=4500.0, doc=5)

#             radar.add_field('ZATT', spec_at)
#             radar.add_field('ZCATT', cor_z)


#             diff_spec_at, cor_zdr = pyart.correct.calculate_attenuation(radar, 0, refl_field=ycfg['orig_zdr_name'],
#                                                             ncp_field=ycfg['orig_sq_name'], rhv_field=ycfg['orig_rho_name'],
#                                                             phidp_field=ycfg['orig_phase_name'], rhv_min=0.7, ncp_min=0.4, 
# #                                                            a_coef=0.0157, beta=0.8, fzl=4500.0)
#                                                             a_coef=0.004, beta=0.006, fzl=4500.0, doc=5)
#             radar.add_field('ZDRATT', diff_spec_at)
#             radar.add_field('ZDRCATT', cor_zdr)





#             fig = plt.figure(figsize=(12, 8))
#             axs = []

#             display = pyart.graph.RadarDisplay(radar)
#             limx = [0, 120]
#             limy = [0, 15]
#             largeDPI = 200

#             ax1 = fig.add_subplot(3, 2, 1)
#                 # first, plot original dBZ
#             display.plot(ycfg['orig_dz_name'], sweep=0, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
#                                         vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
#                                         colorbar_label='Original Z$_H$ (dBZ)', mask_outside=False)
#             ax2 = fig.add_subplot(3, 2, 3)
#             display.plot('ZATT', sweep=0, vmin=0, vmax=0.1, colorbar_label='attenuation (dB/km)', mask_outside=False, cmap=plt.cm.gist_stern_r)


#             ax3 = fig.add_subplot(3, 2, 5)

#             display.plot('ZCATT', sweep=0, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
#                                         vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
#                                         colorbar_label='Att corr Z$_H$ (dBZ)', mask_outside=False)


#             ax4 = fig.add_subplot(3, 2, 2)
#             display.plot(ycfg['orig_zdr_name'], sweep=0, vmin=-1,
#                                         vmax=4, cmap=cfg.plot_params[cfg.zdr_name]['cmap'],
#                                         colorbar_label='Original Z$_{DR}$ (dBZ)', mask_outside=False)

#             ax5 = fig.add_subplot(3, 2, 4)
#             display.plot('ZDRATT', sweep=0, vmin=0, vmax=0.01, colorbar_label="diff'l attenuation (dB/km)", mask_outside=False, cmap=plt.cm.gist_stern_r)


#             ax6 = fig.add_subplot(3, 2, 6)
#             display.plot('ZDRCATT', sweep=0, vmin=-1,
#                                         vmax=4, cmap=cfg.plot_params[cfg.zdr_name]['cmap'],
#                                         colorbar_label='Att corr Z$_{DR}$ (dBZ)', mask_outside=False)



#             for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
#                 a.set_xlim(*limx)
#                 a.set_ylim(*limy)
#                 a.set_title('')

#             plt.tight_layout()

#             plt.savefig('z_attenuation_test.png')


#             quit()




        # I think we wanna keep this on, but just put the original phase (masked) field in DROPS if using that

        dz_um = radtools.extract_unmasked_data(radar, ycfg['orig_dz_name'])
        dp_um = radtools.extract_unmasked_data(radar, phase_name)
    # Range needs to be supplied as a variable, and it needs to be the same shape as dzN, etc.
        rng2d, az2d = np.meshgrid(radar.range['data'], radar.azimuth['data'])


    # Here we're using some CSU radartools code to simultaneously calculate Kdp and standard deviation of the phase
    # which will be used for QCing the data
        kdN, fdN, sdN = csu_kdp.calc_kdp_bringi(dp=dp_um, dz=dz_um, rng=rng2d/1000.0, thsd=12, gs=250.0, window=7, nfilter=2)

        # This is a mask (filter) based on the standard dev of the phase
        sdp_mask = csu_misc.differential_phase_filter(sdN, thresh_sdp=ycfg['sdp_thresh'])
	# sometimes the sdp goes nuts in high reflectivity regions, so let's make sure we pair
	# this mask with relatively low dBZ regions
	sdp_mask = np.logical_and(sdp_mask, radar.fields[ycfg['orig_dz_name']]['data'] <= 10)

	zdr_um = radtools.extract_unmasked_data(radar, ycfg['orig_zdr_name'])

	if 'insect_mask' in ycfg.keys() and ycfg['insect_mask']:

        # Here we're using ZDR and DBZ to calculate an insect filter
            insect_mask = csu_misc.insect_filter(dz_um, zdr_um,  
				dz_range=[[-100, 10], [10, 15], [15, 20], [20, 25], [25, 30], [30, 35]],
				dr_thresh=[1, 1.3, 1.7, 2.1, 2.5, 2.8])

	else:
	    insect_mask = np.zeros((zdr_um.shape), bool)


        bad = -32768
        dz_insect = 1.0 * dz_um
        dz_insect[insect_mask] = bad
        dz_sdp = 1.0 * dz_um
        dz_sdp[sdp_mask] = bad


    # Here we're now combining the masks together
        new_mask = np.logical_or(insect_mask, sdp_mask)
        dz_qc = 1.0 * dz_um
        dz_qc[new_mask] = bad
        # Also doing some despeckling for further QC
        mask_ds = csu_misc.despeckle(dz_qc, ngates=ycfg['despeckle_gates'])
        dz_qc[mask_ds] = bad

    # this total_mask has the 3 masks combined together

        total_mask = np.logical_or(mask_ds, new_mask) # this combines the masks
	# this needs to be paired with the signal to noise ratio
	
#        sq_mask = np.logical_and(radar.fields[sq_name]['data'] < ycfg['sq_thresh'], \
#				radar.fields[ycfg['orig_sn_name']]['data'] < ycfg['sn_thresh'])

        sq_mask = radar.fields[sq_name]['data'] < ycfg['sq_thresh']

	# start with the signal quality mask, and if the signal to noise is there as well, then string that together
	if (ycfg['orig_sn_name'] in radar.fields.keys()) and ('sq_thresh' in ycfg.keys()):
	    sq_mask = np.logical_and(sq_mask, radar.fields[ycfg['orig_sn_name']]['data'] < ycfg['sn_thresh'])

	# now add that to the mask

        total_mask = np.logical_or(total_mask, sq_mask)


        if cc_flag: # if the correlation coefficient is around, use it as an additional mask
            rho_mask = ccorr < ycfg['rho_thresh']
            total_mask = np.logical_or(total_mask, rho_mask)


	if 'plot_masks' in ycfg.keys() and ycfg['plot_masks']:
    #	# this is for diagnostics and plotting the different masks being applied to the data

	    try:

		radar = radtools.add_field_to_radar_object(insect_mask, radar, field_name='INSECTMASK', units='',
					       long_name='Insect mask',
					       standard_name='Insect mask',
					       dz_field=ycfg['orig_dz_name'])

		radar = radtools.add_field_to_radar_object(sdp_mask, radar, field_name='SDPMASK', units='',
					       long_name='SDP mask',
					       standard_name='SDP mask',
					       dz_field=ycfg['orig_dz_name'])

    #	    radar = radtools.add_field_to_radar_object(fdN, radar, field_name='FDPTEMP', units='',
    #					   long_name='Filtered dphase',
    #					   standard_name='Filtered dphase',
    #					   dz_field=ycfg['orig_dz_name'])

		radar = radtools.add_field_to_radar_object(sdN, radar, field_name='SDPTEMP', units='',
					       long_name='sigma dphase',
					       standard_name='sigma dphase',
					       dz_field=ycfg['orig_dz_name'])



		fig = plt.figure(figsize=(10, 12))
		axs = []

		display = pyart.graph.RadarDisplay(radar)
		this_sweep = 1


		ax1 = fig.add_subplot(3, 2, 1)
		# first, plot original dBZ
		display.plot(ycfg['orig_dz_name'], sweep=this_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
					vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
					colorbar_label='Original Z$_H$ (dBZ)', mask_outside=False)

		ax2 = fig.add_subplot(3, 2, 2)
		display.plot(ycfg['orig_zdr_name'], sweep=this_sweep, vmin=cfg.plot_params[cfg.zdr_name]['lims'][0],
					vmax=cfg.plot_params[cfg.zdr_name]['lims'][1], cmap=cfg.plot_params[cfg.zdr_name]['cmap'],
					colorbar_label='Original Z$_{DR}$ (dB)', mask_outside=False)


		ax3 = fig.add_subplot(3, 2, 3)
		display.plot(ycfg['orig_phase_name'], sweep=this_sweep, vmin=0.0, vmax=360.0, cmap=plt.cm.jet,
					colorbar_label='Filtered phase', mask_outside=False)


		ax4 = fig.add_subplot(3, 2, 4)
		display.plot('SDPTEMP', sweep=this_sweep, vmin=0.0, vmax=100.0, cmap=plt.cm.jet,
					colorbar_label='std dev phase', mask_outside=False)


		ax5 = fig.add_subplot(3, 2, 5)
		display.plot('INSECTMASK', sweep=this_sweep, vmin=0, vmax=1, cmap=plt.cm.Purples,
					colorbar_label='Insect mask', mask_outside=False)


		ax6 = fig.add_subplot(3, 2, 6)

		display.plot('SDPMASK', sweep=this_sweep, vmin=0, vmax=1, cmap=plt.cm.Purples,
					colorbar_label='SDP mask', mask_outside=False)




		for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
		    a.grid(True)
		    if radar.scan_type == 'rhi':
			a.set_ylim(0, 18)


		plt.tight_layout()

		ptools.save_figure(fig, base_path, 'raw_masks.png', dpi=largeDPI)

	    except:
		pass




    #print 'rho mask: {}'.format(rho_mask)

        dz_qc[total_mask] = bad
        kdN[total_mask] = bad


        # Doing an initial DZQC and adding it to the radar object
        # because I'll need it in the second trip echo code

        radar = radtools.add_field_to_radar_object(dz_qc, radar, field_name='DZQCI', units='dBZ',
                                       long_name='Reflectivity (Combo Filtered)',
                                       standard_name='Reflectivity (Combo Filtered)',
                                       dz_field=ycfg['orig_dz_name'])



	# make a plot of the raw variables as a diagnostic if needed
	if 'plot_raw' in ycfg.keys() and ycfg['plot_raw']:

	    try:
		fig = plt.figure(figsize=(10, 12))
		axs = []

		this_sweep = 1

		display = pyart.graph.RadarDisplay(radar)

		ax1 = fig.add_subplot(3, 2, 1)
		# first, plot original dBZ
		display.plot(ycfg['orig_dz_name'], sweep=this_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
					vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
					colorbar_label='Original Z$_H$ (dBZ)', mask_outside=False)

		ax2 = fig.add_subplot(3, 2, 2)
		display.plot('DZQCI', sweep=this_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
					vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
					colorbar_label='Original Z$_H$ (dBZ)', mask_outside=False)



		ax3 = fig.add_subplot(3, 2, 3)
		display.plot(ycfg['orig_zdr_name'], sweep=this_sweep, vmin=cfg.plot_params[cfg.zdr_name]['lims'][0],
					vmax=cfg.plot_params[cfg.zdr_name]['lims'][1], cmap=cfg.plot_params[cfg.zdr_name]['cmap'],
					colorbar_label='ZDR (dB)', mask_outside=False)


		ax4 = fig.add_subplot(3, 2, 4)
		display.plot(ycfg['orig_rho_name'], sweep=this_sweep, vmin=ycfg['rho_thresh'],
					vmax=cfg.plot_params[cfg.rho_name]['lims'][1], cmap=cfg.plot_params[cfg.rho_name]['cmap'],
					colorbar_label='CC', mask_outside=False)


		ax5 = fig.add_subplot(3, 2, 5)
		display.plot(ycfg['orig_sq_name'], sweep=this_sweep, vmin=0, vmax=0.8, cmap=plt.cm.nipy_spectral_r,
					colorbar_label='SQI', mask_outside=False)

		try:
		    ax6 = fig.add_subplot(3, 2, 6)
		    display.plot('SNR16', sweep=this_sweep, vmin=-10, vmax=20, cmap=plt.cm.jet,
					colorbar_label='SNR (dB)', mask_outside=False)
		except:
		    pass


		for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
		    a.grid(True)
		    if radar.scan_type == 'rhi':
			a.set_ylim(0, 18)


		plt.tight_layout()

		ptools.save_figure(fig, base_path, 'raw_variables.png', dpi=largeDPI)

	    except:
		pass



        # Okay, now here goes the 2nd trip filtering
        # we'll first do a check to make sure it's in the config file and that it's true

        if 'second_trip_flag' in ycfg.keys() and ycfg['second_trip_flag']:

            print 'Doing second trip echo filtering'
        # let's do a smoothed dBZ field after it's been QCed
            dz_qc_filled = deepcopy(dz_qc)
            dz_qc_filled[dz_qc_filled < -999] = -20.0
            sm_rad = 2

            dz_qc_smooth = ndi.gaussian_filter(dz_qc_filled, (sm_rad, sm_rad), mode='nearest')


            sm_g_rad = 5
        # now let's do a plot of smoothed spectral width
            width_smooth = ndi.gaussian_filter(radar.fields[ycfg['orig_width_name']]['data'].filled(2.0), (sm_g_rad, sm_g_rad), mode='nearest')


        # now let's do a smoothed SQI
            sqi_smooth = ndi.gaussian_filter(radar.fields[ycfg['orig_sq_name']]['data'], (sm_g_rad, sm_g_rad), mode='nearest')


        # now let's try a combo filter between SQI and spectrum width to see if we can 
        # isolate 2nd trip echoes


            second_trip_flag = np.ma.getdata((sqi_smooth <= 0.12))
            stf_int = second_trip_flag.astype(int)


            total_mask = np.logical_or(total_mask, second_trip_flag)

            # can plot diagnostics here to make sure everything is working
            if 'plot_second_trip' in ycfg.keys() and ycfg['plot_second_trip']:


		try:

		    radar = radtools.add_field_to_radar_object(dz_qc_smooth, radar, field_name='DZQCSMOOTH', units='dBZ',
					   long_name='Reflectivity (Combo Filtered)',
					   standard_name='Reflectivity (Combo Filtered)',
					   dz_field='DZQCI')


		    radar = radtools.add_field_to_radar_object(width_smooth, radar, field_name='WIDTHSMOOTH', units='---',
					   long_name='-------', standard_name='------', dz_field='DZQCI')


		    radar = radtools.add_field_to_radar_object(sqi_smooth, radar, field_name='SQISMOOTH', units='---',
					   long_name='-------', standard_name='------', dz_field='DZQCI')

	    # now let's try a combo filter between SQI and spectrum width to see if we can 
	    # isolate 2nd trip echoes

		    radar = radtools.add_field_to_radar_object(stf_int, radar, field_name='STF', units='---',
					   long_name='-------', standard_name='------', dz_field='DZQCI')


	    # okay, now try to plot the qc'ed data and then plot the dBZ of just the second trip
		    fig = plt.figure(figsize=(10, 12))
		    axs = []

		    display = pyart.graph.RadarDisplay(radar)
		    #limx = [0, 120]
		    #limy = [0, 15]
		    this_sweep = 1 

		    ax1 = fig.add_subplot(3, 2, 1)
		    # first, plot original dBZ
		    display.plot(ycfg['orig_dz_name'], sweep=this_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
					    vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
					    colorbar_label='Original Z$_H$ (dBZ)', mask_outside=False)
		    ax2 = fig.add_subplot(3, 2, 2)
		    display.plot('DZQCSMOOTH', sweep=this_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
					    vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
					    colorbar_label='Smoothed Z$_H$ (dBZ)', mask_outside=False)



		    ax3 = fig.add_subplot(3, 2, 3)
		    # # now plot the QCed dBZ
		    # display.plot('WIDTH2', sweep=0, vmin=0, vmax=4, colorbar_label='Spec width', mask_outside=False)
		    display.plot('SQISMOOTH', sweep=this_sweep, vmin=0, vmax=1, cmap=plt.cm.nipy_spectral_r,
                                        colorbar_label='SQI', mask_outside=False)


		    ax4 = fig.add_subplot(3, 2, 4)
		    # now plot the QCed dBZ
		    display.plot('WIDTHSMOOTH', sweep=this_sweep, vmin=0, vmax=2, colorbar_label='Spec width smooth', mask_outside=False)


		    ax5 = fig.add_subplot(3, 2, 5)
		    display.plot('STF', sweep=this_sweep, vmin=0, vmax=1, colorbar_label='Second trip echoes', mask_outside=False)

		    if ycfg['orig_sn_name'] in radar.fields.keys():
			ax6 = fig.add_subplot(3, 2, 6)
			display.plot(ycfg['orig_sn_name'], sweep=0, vmin=-10, vmax=20, cmap=plt.cm.jet, 
					     colorbar_label='SNR (dB)', mask_outside=False)


		    plt.tight_layout()

		    ptools.save_figure(fig, base_path, 'secondtrip_test.png', dpi=largeDPI)

		except:
		    pass


        dz_qc[total_mask] = bad
        kdN[total_mask] = bad








    # ********** try attenuation correction on reflectivity??? 
    # *************

        radar = radtools.add_field_to_radar_object(dz_qc, radar, field_name='DZQC', units='dBZ', 
                                       long_name='Reflectivity (Combo Filtered)',
                                       standard_name='Reflectivity (Combo Filtered)', 
                                       dz_field=ycfg['orig_dz_name'])




    # now alter certain fields
        filter_fields = [cc_var, ycfg['orig_zdr_name'], vel_name, ycfg['orig_kdp_name'], cfg.ph_name]

        for ff in filter_fields:
            old_var_mask = deepcopy(radar.fields[ff]['data'].mask)
            new_var_mask = np.logical_or(old_var_mask, total_mask)
            radar.fields[ff]['data'].mask = new_var_mask

    # Adding some of the things we've calculated to the radar object


    # making sure the phase is between -180 and +180, will be plotting now
        fdN[fdN<-180.0] += 360.0
        fdN[fdN>180.0] -= 360.0
        fdN[total_mask] = bad


        radar = radtools.add_field_to_radar_object(kdN, radar, field_name='KDP', units='deg/km', 
                                       long_name='CSU Specific Differential Phase',
                                       standard_name='CSU Specific Differential Phase', 
                                       dz_field='DZQC')
        radar = radtools.add_field_to_radar_object(fdN, radar, field_name='FDP', units='deg', 
                                       long_name='CSU Filtered Differential Phase',
                                       standard_name='CSU Filtered Differential Phase', 
                                       dz_field='DZQC')
        radar = radtools.add_field_to_radar_object(sdN, radar, field_name='SDP', units='deg', 
                                       long_name='CSU Standard Deviation of Differential Phase',
                                       standard_name='CSU Standard Deviation of Differential Phase', 
                                       dz_field='DZQC')


    if 'att_corr' in ycfg.keys() and ycfg['att_corr']:
    # first, let's start off by doing some reflectivity (and differential reflectivity)
    # attenuation correction, cuz it makes the most sense to do that before doing anything else

        gate_heights = deepcopy(radar.get_gate_x_y_z(0)[-1]/1000.0)
        for ins in range(1, radar.nsweeps):
            gate_heights = np.concatenate((gate_heights, radar.get_gate_x_y_z(ins)[-1]/1000.0)) # want it in km

	# these gate heights are arranged the same way as the radar.fields are

	# **** NEED TO LOOK INTO THE A VALUES THAT I WANT HERE

        zcorr, satt, zaa = acorr.carey_atten_corr(radar.fields[ycfg['qc_dz_name']]['data'], radar.fields[ycfg['phase_name']]['data'], 
                dr=radar.fields[ycfg['qc_zdr_name']]['data'], kd=radar.fields[ycfg['kdp_name']]['data'], 
                rh=radar.fields[ycfg['rho_name']]['data'], height=gate_heights, a_default=-0.15, deva=0.10, 
            bad=np.nan, force_default=False, phase_diff_calc=True, start_ind=radar.sweep_start_ray_index['data'],
            end_ind=radar.sweep_end_ray_index['data'], min_pts=150)


        zcorr = np.ma.masked_invalid(zcorr)
        satt = np.ma.masked_invalid(satt)


        radar = radtools.add_field_to_radar_object(zcorr, radar, field_name='DZQCC', units='dBZ',
                                   long_name='Reflectivity (Filtered/Attn corr)',
                                   standard_name='Reflectivity (Filtered/Attn corr)',
                                   dz_field='DZQC')


        zdrcorr, zdr_satt, zdraa = acorr.carey_atten_corr(radar.fields[ycfg['qc_zdr_name']]['data'], radar.fields[ycfg['phase_name']]['data'], 
                dr=radar.fields[ycfg['qc_zdr_name']]['data'], kd=radar.fields[ycfg['kdp_name']]['data'], 
                rh=radar.fields[ycfg['rho_name']]['data'], height=gate_heights, a_default=-0.018, deva=0.005, 
            bad=np.nan, force_default=False, phase_diff_calc=True, start_ind=radar.sweep_start_ray_index['data'],
            end_ind=radar.sweep_end_ray_index['data'], min_pts=150)


        zdrcorr = np.ma.masked_invalid(zdrcorr)
        zdr_satt = np.ma.masked_invalid(zdr_satt)


        radar = radtools.add_field_to_radar_object(zdrcorr, radar, field_name='ZDRC', units='dBZ',
                                   long_name='Diffl Reflectivity (Attn corr/Combo Filtered)',
                                   standard_name='Diffl Reflectivity (Attn corr/Combo Filtered)',
                                   dz_field='DZQC')



        zatt_sweep = 0



        if radar.scan_type == 'rhi':

            display = pyart.graph.RadarDisplay(radar)
            limx = [0, 120]
            limy = [0, 15]
            nrows = 6
            ncols = 1
            figx = 12
            figy = 12
            

        else:
            display = pyart.graph.RadarDisplay(radar)
            limx = [-120, 120]
            limy = [-120, 120]
            nrows = 3
            ncols = 2
            figx = 10
            figy = 12

        fig = plt.figure(figsize=(figx, figy))
        axs = []


        ax1 = fig.add_subplot(nrows, ncols, 1)
            # first, plot original dBZ
        display.plot(ycfg['qc_dz_name'], sweep=zatt_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
                                    vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
                                    colorbar_label='Orig Z$_H$', mask_outside=False)
        ax2 = fig.add_subplot(nrows, ncols, 2)
        display.plot('DZQCC', sweep=zatt_sweep, vmin=cfg.plot_params[cfg.dz_name]['lims'][0],
                                    vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'],
                                    colorbar_label='Corr Z$_H$', mask_outside=False)

        ax3 = fig.add_subplot(nrows, ncols, 3)
        display.plot(ycfg['phase_name'], sweep=zatt_sweep, vmin=-180,
                                    vmax=0, cmap=cfg.plot_params[cfg.ph_name]['cmap'],
                                    colorbar_label='Diffl phase', mask_outside=False)


        ax4 = fig.add_subplot(nrows, ncols, 4)
        display.plot(ycfg['kdp_name'], sweep=zatt_sweep, vmin=-1,
                                    vmax=6, cmap=cfg.plot_params[cfg.kdp_name]['cmap'],
                                    colorbar_label='Kdp', mask_outside=False)


        ax5 = fig.add_subplot(nrows, ncols, 5)
        display.plot(ycfg['qc_zdr_name'], sweep=zatt_sweep, vmin=cfg.plot_params[cfg.zdr_name]['lims'][0],
                                    vmax=cfg.plot_params[cfg.zdr_name]['lims'][1], cmap=cfg.plot_params[cfg.zdr_name]['cmap'],
                                    colorbar_label='Orig Z$_DR$', mask_outside=False)

        ax6 = fig.add_subplot(nrows, ncols, 6)
        display.plot('ZDRC', sweep=zatt_sweep, vmin=cfg.plot_params[cfg.zdr_name]['lims'][0],
                                    vmax=cfg.plot_params[cfg.zdr_name]['lims'][1], cmap=cfg.plot_params[cfg.zdr_name]['cmap'],
                                    colorbar_label='Corr Z$_DR$', mask_outside=False)


        #print np.unique(radar.fields['FDP']['data'][:10])
        #print np.unique(radar.fields['FDP']['data'][-10:])


        for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
            a.set_xlim(*limx)
            a.set_ylim(*limy)
            a.set_title('')
            a.set_ylabel('')

        fig.suptitle('%s %s Attenuation correction'%(prefix, file_dt_string), fontsize=18)

        plt.tight_layout()
        fig.subplots_adjust(top=0.94)

        plt.savefig('z_carey_attenuation_test.png', dpi=DPI)




    if unfold_flag:
        # try to unfold the velocities
        corrected_velocity = pyart.correct.dealias_unwrap_phase(radar, vel_field=vel_name, corr_vel_field='CV')

        radar = radtools.add_field_to_radar_object(corrected_velocity['data'], radar, field_name='CV', units='m/s', 
                                       long_name='Unfolded radial velocity',
                                       standard_name='Corrected velocity', 
                                       dz_field='DZQC')

        if radar.scan_type == 'rhi':


            lsweep = 0
            display = pyart.graph.RadarDisplay(radar)
            limx = [0, 120]
            limy = [0, 15]
            vmax = 20.0


            fig = plt.figure(figsize=(11, 5))
            axs = []

            axs.append(fig.add_subplot(1, 2, 1))
            display.plot(vel_name, sweep=lsweep, vmin=-1*vmax, vmax=vmax, cmap=ctables.Carbone11, 
                 colorbar_label='Radial velocity (m/s)', mask_outside=False)
            display.set_limits(xlim=limx, ylim=limy)

            axs.append(fig.add_subplot(1, 2, 2))

            display.plot('CV', sweep=lsweep, vmin=-1*vmax, vmax=vmax, 
                cmap=ctables.Carbone11, colorbar_label='Corrected velocity (m/s)', mask_outside=False)
            display.set_limits(xlim=limx, ylim=limy)



            for a in axs:
                a.set_title('')


            angle_title = np.average(radar.get_azimuth(0))

            plt.tight_layout()

#fig.suptitle('SEAPOL %s %s %.1f$^{\circ}$'%(file_dt_string, radar.scan_type.upper(), angle_title), fontsize=13)
            fig.subplots_adjust(top=0.90)

            #print 'saving large %s figure'%(radar.scan_type.upper())

    #plt.savefig('%s/figures/large_dbz_vel_%s/%s%s_%s.png'%(base_path, radar.scan_type, prefix, dt_string, radar.scan_type), dpi=largeDPI)
            ptools.save_figure(fig, base_path, 'veltest.png', dpi=largeDPI)

    #plt.close(fig)


	# this is basically for diagnostic purposes, will plot original/filtered dBZ and phase
	# trying to figure out how to get rid of second-trip echoes
    if 'plot_dz_diag' in ycfg.keys() and ycfg['plot_dz_diag'] and radar.scan_type == 'ppi':
        # if it's there and it's true
        fig = plt.figure(figsize=(11, 10))
        axs = []

        display = pyart.graph.RadarDisplay(radar)
        #limx = [0, 120]
        #limy = [0, 15]

        ax1 = fig.add_subplot(2, 2, 1)
        # first, plot original dBZ
        display.plot(ycfg['orig_dz_name'], sweep=0, vmin=cfg.plot_params[cfg.dz_name]['lims'][0], 
                                        vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'], 
                                        colorbar_label='Original Z$_H$ (dBZ)', mask_outside=False)

        ax2 = fig.add_subplot(2, 2, 2)
        # now plot the QCed dBZ
        display.plot(ycfg['qc_dz_name'], sweep=0, vmin=cfg.plot_params[cfg.dz_name]['lims'][0], 
                                        vmax=cfg.plot_params[cfg.dz_name]['lims'][1], cmap=cfg.plot_params[cfg.dz_name]['cmap'], 
                                        colorbar_label='QC Z$_H$ (dBZ)', mask_outside=False)

        ax3 = fig.add_subplot(2, 2, 3)
        # now plot the QCed dBZ
        display.plot(ycfg['orig_phase_name'], sweep=0, vmin=cfg.plot_params[cfg.ph_name]['lims'][0], 
                                        vmax=cfg.plot_params[cfg.ph_name]['lims'][1], cmap=cfg.plot_params[cfg.ph_name]['cmap'], 
                                        colorbar_label='Original phase (deg)', mask_outside=False)

        ax4 = fig.add_subplot(2, 2, 4)
        # now plot the QCed dBZ
        display.plot('FDP', sweep=0, vmin=cfg.plot_params[cfg.ph_name]['lims'][0], 
                                        vmax=cfg.plot_params[cfg.ph_name]['lims'][1], cmap=cfg.plot_params[cfg.ph_name]['cmap'], 
                                        colorbar_label='Filtered phase (deg)', mask_outside=False)


        plt.tight_layout()

        ptools.save_figure(fig, base_path, 'dbz_comparison.png', dpi=largeDPI)



        # now remove some fields from the radar object before writing it out
    if 'qc_remove_fields' in ycfg.keys():
        # now go thru each one and remove it from radar if it's in there
        for qrf in ycfg['qc_remove_fields']:
            if qrf in radar.fields.keys():
                radar.fields.pop(qrf)        

    print 'Adding HID'
    try:
	snd_data = np.genfromtxt('%s/sample_tropical_sounding.txt'%(base_path), dtype=None, skip_header=6)
	snd_ht = snd_data['f1'].astype(float)
	snd_temp = snd_data['f2']

	radar_T, radar_z = radtools.interpolate_sounding_to_radar(snd_temp, snd_ht, radar)



	scores = csu_fhc.csu_fhc_summer(dz=radar.fields[cfg.dz_name]['data'], zdr=radar.fields[cfg.zdr_name]['data'],
		    rho=radar.fields[cfg.rho_name]['data'], kdp=radar.fields[cfg.kdp_name]['data'], 
			use_temp=True, band='C', T=radar_T)
	fh = np.argmax(scores, axis=0) + 1


	radar = radtools.add_field_to_radar_object(fh, radar, field_name='HID', units='',
				       long_name='CSU Dolan 2013 dominant hydrometeor classification',
				       standard_name='CSU Hydrometeor type',
				       dz_field=cfg.dz_name)

    except Exception, hide:
	print 'Error with HID: {}'.format(hide)




    #print 'radar.fields as about to save: {}'.format(radar.fields)
    print 'final radar fields: {}'.format(radar.fields.keys())

    radar_base = os.path.basename(radarfile)
    dot_loc = radar_base.rfind('.')
    new_radar_name = '%s%s_%s.nc'%(prefix, dt_string, radar.scan_type)
    #prefix, dt_string, radar.scan_type

    pyart.io.cfradial.write_cfradial('%s/%s/%s'%(outpath, radar.scan_type, new_radar_name), radar, 
                            format='NETCDF4', time_reference=False)
    print 'saved cfradial %s/%s/%s'%(outpath, radar.scan_type, new_radar_name)


    # ************* Now adding a UF file for DROPS purposes *******************
#        new_uf_name = new_radar_name.replace('nc', 'uf')
#        uf_path = outpath.replace('cfradial', 'uf')
#
#        print 'before uf radar fields: {}'.format(radar.fields.keys())
#
#        pyart.io.write_uf('%s/%s/%s'%(uf_path, radar.scan_type, new_uf_name), radar)


else:
        print 'File is not valid, not plotting or saving'











