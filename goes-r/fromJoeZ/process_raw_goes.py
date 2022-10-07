'''
process_raw_goes.py

Processes raw goes files downloaded from NOAA's CLASS archive
Saves brightness temperature w/ scene temp correction to netCDF

http://rammb.cira.colostate.edu/research/calibration_validation_and_visualization/
goes_image_display/conversion_ir_10_to_8.asp
'''

import numpy as np
import pandas as pd
import xarray as xr
import glob
import netCDF4 as nc4
import pdb
from goes_util import get_goes_constants
from datetime import datetime, timedelta

indir = '/home/disk/anvil2/ensembles/goes/'
outdir = '/home/disk/meso-home/jzagrod/Satellite/GOES/IR/'
domain = 'ne_pacific'

def find_files(indir):
    return glob.glob('{}goes*'.format(indir))

def convert_to_tb(counts,goes_num,channel):
    #convert data to brightness temp
    counts = (counts/256.) #reduce to 8-bit
    #counts = (counts/256).astype('uint8') 
    gcs = get_goes_constants(goes_num,channel)
    radiance = (counts - gcs['scaling_bias'])/gcs['scaling_gain']
    teff = (gcs['c2']*gcs['nu'])/np.log(1+((gcs['c1']*gcs['nu']**3)/radiance))
    tb_tmp = gcs['a']+(gcs['b']*teff)
    #tscene = np.zeros((np.shape(tb)))

    #tscene[np.where((counts >= 176))] = 418 - tb[np.where((counts >= 176))]
    #tscene[np.where((counts < 176))] = 660 - 2*tb[np.where((counts < 176))]
    tb = 418 - tb_tmp
    #tb = tb - 273.15
    return tb

files = find_files(indir)

for infile in files:
    f = nc4.Dataset(infile,'r') 
    goes_num = infile.split('/')[6].split('.')[0]
    year = int(infile.split('/')[6].split('.')[1])
    doy = int(infile.split('/')[6].split('.')[2])
    hour = int(infile.split('/')[6].split('.')[3][2:4])
    minute = int(infile.split('/')[6].split('.')[3][4:])
    goes_time = datetime(year,1,1,hour,minute)
    goes_time += timedelta(days = doy-1)
    timestr = goes_time.strftime('%Y%m%d_%H%M')

    channel = int(infile.split('/')[6].split('_')[1].split('.')[0])
    if channel == 4:
        goes_type = 'ir'

    lat_ir = f['lat'][:]
    lon_ir = f['lon'][:]
    counts = f['data'][0,:,:]
    tb = convert_to_tb(counts,goes_num,channel)
    lat_ir[np.where((np.isnan(tb)==True))]= float('nan')
    lon_ir[np.where((np.isnan(tb)==True))]= float('nan')

    #create xarray of data
    ds = xr.Dataset({'tb': (['x', 'y'],  tb)},
                coords={'lon': (['x', 'y'], lon_ir),
                        'lat': (['x', 'y'], lat_ir)})
    ds['tb'].attrs = {'satellite': goes_num,
                                    'channel': channel,
                      'time' : timestr}
    ds.to_netcdf('{}goes_{}_{}_{}.nc'.format(outdir,goes_type,domain,timestr), 'w')
    pdb.set_trace()

    
