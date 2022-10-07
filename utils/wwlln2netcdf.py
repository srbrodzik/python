# =============================== wwlln2netcdf =============================== #
# Convert World Wide Lightning Location Network (WWLLN) data from ASCII to
# NetCDF with compression.
# 
# Initial code by Joe Finlon, University of Washington. Copyright 2021.
# 
# Code History
# ----------
# 09/24/2021: Initial code commit [JF]
# 
# Parameters
# ----------
# start_date: start date in time range in YYYY-mm-dd format [str]
# end_date: end date in time range in YYYY-mm-dd format, inclusive [str]
#
# Execution
# ----------
# python wwlln2netcdf.py 2019-01-01 2019-12-31
# 
# Note: Dates that span multiple calendar years not currently supported.
# See the Directory Paths section (line 45) if you wish to change the base
#     input and output directories.
# ============================================================================ #
import os
import sys
import numpy as np
import xarray as xr
from datetime import datetime
import pandas as pd

# --- SUBROUTINES --- #
def date_parser(d_bytes):
    s = d_bytes.decode('utf-8')
    return np.datetime64(datetime.strptime(s, '%Y/%m/%d'))

def tis_parser(d_bytes):
    s = d_bytes.decode('utf-8')
    return (3600. * float(s[:2]) + 60. * float(s[3:5]) + float(s[6:8])
            + float(s[9:]) / 10.**6)

def dt(time_info):
    return (time_info[0] + np.timedelta64(int(np.fix(time_info[1])), 's')
            + np.timedelta64(int(10.**6 * np.modf(time_info[1])[0]), 'us'))

# --- DIRECTORY PATHS --- #
indir = '/home/disk/meso-home/jfinlon/lightning/AEfiles/global/'
outdir = '/home/disk/radar/data/wwlln/GlobalWithAdditionalFields/'

# --- READ USER INPUT DATE RANGE --- #
start_date = sys.argv[1]
end_date = sys.argv[2]
date_range = pd.date_range(start_date, end_date).to_pydatetime()

# --- LOOP THROUGH DATES IN DATE RANGE --- #
for i in range(len(date_range)):
    date_str = datetime.strftime(date_range[i], '%Y%m%d')
    infile = (indir + start_date[:4] + '/AE' + date_str + '.loc')
    outfile = (outdir + start_date[:4] + '/AE' + date_str + '.nc')
    if os.path.isfile(infile):
        time_data = np.genfromtxt(
            infile, delimiter=',', usecols=[0, 1],
            converters={0:date_parser, 1:tis_parser})
        time = np.array(
            [dt(time_data[i]) for i in range(len(time_data))],
            dtype='datetime64[us]')
        strike_data = np.genfromtxt(
            infile, delimiter=',', usecols=[2, 3, 4, 5, 6, 7, 8])

        lat = strike_data[:, 0]
        lon = strike_data[:, 1]
        res = strike_data[:, 2]
        nstn = strike_data[:, 3].astype(np.short)
        energy = strike_data[:, 4]
        energy_uncert = strike_data[:, 5]
        nstn_energy = strike_data[:, 6].astype(np.short)

        ds = xr.Dataset()

        ds.attrs['title'] = 'World Wide Lightning Location Network (WWLLN) Global Data'
        ds.attrs['institution'] = 'University of Washington'
        ds.attrs['source'] = 'ascii global lightning strikes'
        ds.attrs['start_time'] = np.datetime_as_string(time[0])
        ds.attrs['references'] = 'https://wwlln.net/'
        ds.attrs['comment'] = 'Raw data converted to NetCDF using wwlln2netcdf.py'

        ds.coords['time'] = ('time', time)

        ds['lon'] = ('time', lon)
        ds['lon'].attrs['valid_min'] = -180.
        ds['lon'].attrs['valid_max'] = 180.
        ds['lon'].attrs['standard_name'] = 'longitude'
        ds['lon'].attrs['units'] = 'degrees_east'

        ds['lat'] = ('time', lat)
        ds['lat'].attrs['valid_min'] = -90.
        ds['lat'].attrs['valid_max'] = 90.
        ds['lat'].attrs['standard_name'] = 'latitude'
        ds['lat'].attrs['units'] = 'degrees_north'

        ds['resid'] = ('time', res)
        ds['resid'].attrs['valid_min'] = 0.
        ds['resid'].attrs['valid_max'] = 30.
        ds['resid'].attrs['standard_name'] = 'resid'
        ds['resid'].attrs['long_name'] = 'residual_fit_error'
        ds['resid'].attrs['units'] = 'microseconds'

        ds['nstn'] = ('time', nstn)
        ds['nstn'].attrs['valid_min'] = 0.
        ds['nstn'].attrs['valid_max'] = 32767.
        ds['nstn'].attrs['standard_name'] = 'nstn'
        ds['nstn'].attrs['long_name'] = 'number_stations'
        ds['nstn'].attrs['description'] = (
            'number of WWLLN stations participating in the location fit')
        ds['nstn'].attrs['units'] = ''

        ds['E'] = ('time', energy)
        ds['E'].attrs['standard_name'] = 'E'
        ds['E'].attrs['long_name'] = 'rms_energy'
        ds['E'].attrs['description'] = (
            'root mean square energy of the stroke from 1.3 ms waveform '
            + 'sampling between 7 and 18 kHz')
        ds['E'].attrs['units'] = 'J'

        ds['E_error'] = ('time', energy_uncert)
        ds['E_error'].attrs['standard_name'] = 'E_error'
        ds['E_error'].attrs['long_name'] = 'energy_uncertainty'
        ds['E_error'].attrs['description'] = 'energy error of the location fit'
        ds['E_error'].attrs['units'] = 'J'

        ds['nstn_sub'] = ('time', nstn_energy)
        ds['nstn_sub'].attrs['valid_min'] = 0.
        ds['nstn_sub'].attrs['valid_max'] = 32767.
        ds['nstn_sub'].attrs['standard_name'] = 'nstn_sub'
        ds['nstn_sub'].attrs['long_name'] = 'number_stations_energy'
        ds['nstn_sub'].attrs['description'] = (
            'subset of nstn stations between 1000 and 8000 km from the stroke '
            + 'whose energy data were used in the energy estimate')
        ds['nstn_sub'].attrs['units'] = ''

        fillvals = {}
        ds.fillna(value=fillvals)

        comp = dict(zlib=True, complevel=5, _FillValue= -9999)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(outfile, encoding=encoding)
    else:
        print('File {} does not exist. Skipping {}.'.format(infile, date_str))