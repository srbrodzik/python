# Brody Fuchs, Sept 2017, brfuchs@atmos.colostate.edu


# BF: 10/19/18: This code will read in a QCed cfradial file and 
# just choose the variables that are passed into the command line

# I tried to poke around in subsampling the data in both range/rays
# and it's gonna be quite a pain in the butt. 
# Would either have to go in and modify the pyart source code to 
# modify things when the file is first being read in
# or would basically have to make my own netcdf from the data 
# that is read in. Not being able to alter netcdf files after they've been
# created is quite annoying. I think I could do this in the field
# if we're in a pinch, but I'd rather not!



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
parser.add_argument('--config', action="store", dest="config", default=None)
parser.add_argument('--vars', nargs='+', type=str, dest='vars', default=[])

pargs = parser.parse_args()

if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v

if 'downsample_path' in ycfg.keys():
    outpath = ycfg['downsample_path']
else:
    outpath = '.'


prefix = ycfg['qc_prefix']

rb = pargs.file

radarfile = deepcopy(rb)

rbase = os.path.basename(radarfile)

a_index = rbase.rfind('A')

if '.' in rbase:
    last_dot = rbase.rfind('.')
    dt_string = rbase[a_index+1:last_dot-4]
else:
    last_dot = None
    dt_string = rbase[a_index+1:-4]

#print dt_string
file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
file_dt_string = file_time.strftime('%Y-%m-%d %H:%M:%S')

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

print 'radar.fields.keys: {}'.format(radar.fields.keys())


phase_name = ycfg['orig_phase_name']
vel_name = ycfg['orig_vel_name']
sq_name = ycfg['orig_sq_name']

rscantype = radar.scan_type
nsweeps = radar.nsweeps


print 'Radar scan type: {}, number of sweeps: {}'.format(rscantype, nsweeps)


if True:


    radar_lat = radar.latitude['data'][0]
    radar_lon = radar.longitude['data'][0]


    # keys that need to be removed from the radar object
    #bad_keys = list(set(radar.fields.keys()) & set(lst2))

    for k in radar.fields.keys():
        if k not in pargs.vars:
            radar.fields.pop(k, None)




    #print 'radar.fields as about to save: {}'.format(radar.fields)
    print 'output radar fields: {}'.format(radar.fields.keys())

    radar_base = os.path.basename(radarfile)
    dot_loc = radar_base.rfind('.')
    new_radar_name = '%s%s_%s.nc'%(prefix, dt_string, radar.scan_type)
    #prefix, dt_string, radar.scan_type

    pyart.io.cfradial.write_cfradial('%s/%s'%(outpath, new_radar_name), radar, 
                            format='NETCDF4', time_reference=False)
    print 'saved cfradial %s/%s'%(outpath, new_radar_name)


    # ************* Now adding a UF file for DROPS purposes *******************
#        new_uf_name = new_radar_name.replace('nc', 'uf')
#        uf_path = outpath.replace('cfradial', 'uf')
#
#        print 'before uf radar fields: {}'.format(radar.fields.keys())
#
#        pyart.io.write_uf('%s/%s/%s'%(uf_path, radar.scan_type, new_uf_name), radar)


else:
        print 'File is not valid, not plotting or saving'











