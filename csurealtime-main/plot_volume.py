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
import gentools as gtools
import plot_tools as ptools
import radar_tools as radtools


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
parser.add_argument('--var', action="store", dest="var", default=None)
parser.add_argument('--config', action="store", dest="config", default=None)

pargs = parser.parse_args()

if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v



if pargs.var is not None:
    plot_var = ycfg[pargs.var]
else:
    plot_var = 'DBZ2'


print 'plot variable: {}'.format(plot_var)

outpath = '%s/figures/volume'%(base_path)
#inpath = '%s/radar_data/seapol'%(base_path)
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


radarfile = deepcopy(pargs.file)

rbase = os.path.basename(radarfile)


a_index = rbase.rfind(prefix[-1])

if '.' in rbase:
    last_dot = rbase.rfind('.')
    try:
	dt_string = rbase[a_index+1:last_dot]
        file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
    except ValueError:
        dt_string = rbase[a_index+1:last_dot-4]
        file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')

else:
	last_dot = None
	dt_string = rbase[a_index+1:]

#nsweeps_thresh = ycfg['nsweeps_thresh']

#print dt_string
file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')
date_string = file_time.strftime('%Y%m%d')


major_circle_rads = np.arange(20, 140, 20)
minor_circle_rads = np.arange(10, 110, 20)
azimuths = np.arange(0, 360, 45)
minor_azimuths = np.arange(22.5, 360, 45)




# Read in the UF file with pyart here, we keep the file field names for simplicity: these are the 2
# character names mandated by UF convention

cfg = RadarConfig.RadarConfig(dz=ycfg['dz_name'], zdr=ycfg['zdr_name'], kdp=ycfg['kdp_name'], rho=ycfg['rho_name'], 
				hid='HID', vel=ycfg['vel_name'])

if 'uf' in rbase:
	radar = pyart.io.read_uf(radarfile, file_field_names=True)
	#cfg = RadarConfig.RadarConfig(dz='DZQC', zdr='DR', kdp='KD', rho='RH', hid='HID')


else:
	radar = pyart.io.read(radarfile, file_field_names=True)
	#cfg = RadarConfig.RadarConfig(dz='DBZ2', zdr='ZDR2', kdp='KDP2', rho='RHOHV2', hid='HID')

print 'Plotting {} volume'.format(radar.scan_type)


phase_name = ycfg['phase_name']
vel_name = ycfg['vel_name']
sq_name = ycfg['sq_name']


radar_lat = radar.latitude['data'][0]
radar_lon = radar.longitude['data'][0]


radar_base = os.path.basename(radarfile)
dot_loc = radar_base.rfind('.')
new_radar_name = '%s%s_%s.nc'%(prefix, dt_string, radar.scan_type)
#prefix, dt_string, radar.scan_type


display = pyart.graph.RadarDisplay(radar)

if radar.scan_type == 'ppi':

    limx = [-120, 120]
    limy = [-120, 120]
    sweep = 0
    display_func = display.plot_ppi

elif radar.scan_type == 'rhi':
    limx = [0, 120]
    limy = [0, 15]
    sweep = 0
    display_func = display.plot


ncols = np.floor(np.sqrt(radar.nsweeps))
nrows = np.ceil(radar.nsweeps/ncols)
im_scale = 2.5

fig = plt.figure(figsize=(im_scale*ncols*1.25, im_scale*nrows))


axs = []

for ia in range(radar.nsweeps):


    axs.append(fig.add_subplot(nrows, ncols, ia+1))
    display_func(plot_var, sweep=ia, vmin=cfg.plot_params[plot_var]['lims'][0], vmax=cfg.plot_params[plot_var]['lims'][1], cmap=cfg.plot_params[plot_var]['cmap'], 
                mask_outside=True, colorbar_label='')
    display.set_limits(xlim=limx, ylim=limy)
    #display.generate_title(plot_var, ia)




if radar.scan_type == 'ppi':

    for j, a in enumerate(axs):
        a.set_xlabel('')
        a.set_ylabel('')

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

    #angle_title = np.average(radar.get_elevation(sweep))
        angle_title = np.average(radar.get_elevation(j))
    #cur_ax = plt.gca()
        a.set_title('%.1f$^{\circ}$'%(angle_title))


elif radar.scan_type == 'rhi':
    for j, a in enumerate(axs):
        angle_title = np.average(radar.get_azimuth(j))
        #cur_ax = plt.gca()
        a.set_title('%.1f$^{\circ}$'%(angle_title))





plt.tight_layout()

fig.suptitle('%s %s %s %s volume'%(radar_title, file_dt_string, radar.scan_type.upper(), plot_var), fontsize=14)
fig.subplots_adjust(top=0.90)

print 'saving %s volume figure'%(radar.scan_type.upper())
print '%s/volume/%s'%(fig_path, date_string)


ptools.save_figure(fig, '%s/volume/%s'%(fig_path, date_string), '%s%s_%s.png'%(prefix, dt_string, radar.scan_type), dpi=DPI)







