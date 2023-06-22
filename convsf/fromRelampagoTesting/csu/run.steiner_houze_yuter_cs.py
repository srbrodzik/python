import netCDF4 as nc4
import numpy as np
import sys
import os
from steiner_houze_yuter_cs import assign_radius
from steiner_houze_yuter_cs import haversine
from steiner_houze_yuter_cs import conv_strat_latlon
from netcdf_io import writeCFnetcdf

in_dir = '/home/storm/brodzik/data/radarCart/RMA1/latlon'
out_dir = '/home/storm/brodzik/data/raintype/RMA1/latlon/steiner_csu'
dates = ['20181111']
refl_name = 'DBZ'
lat_name = 'y0'
lon_name = 'x0'
alt_name = 'z0'
refl_level = 4
lat_origin = -31.44133
lon_origin = -64.19192

# Input params
core_thresh = 42
method = 'SYH'
a = 8
b = 64
tune_thresh = 42.0
sm_rad = 15
fill_dbz = -10  # brody used 25
bg_diff = 10

debug = 1

# Information about output
institution = 'University of Washington'
source = 'Cordoba AR radar data'
title = 'Rain type classifications'
references = 'CSU coding of Steiner algorithm: steiner_houze_yuter_cs.py'
comment = 'Based on 2.5km level of interpolated reflectivity data'

for date in dates:
    for fname in os.listdir(in_dir+'/'+date):
        if fname.endswith('nc'):

            if debug:
                print >>sys.stderr, 'fname = ', fname

            # Filename for output
            fname_out = str(out_dir+'/'+fname.replace('ncf','raintype') )
            if debug:
                print >>sys.stderr, 'fname_out = ', fname_out

            #If output dir does not exist, create it
            dir = os.path.dirname(fname_out)
            if not os.path.exists(dir):
                os.makedirs(dir)

            #Open input file
            ncid = nc4.Dataset(str(in_dir+'/'+date+'/'+fname),'r')

            #Make sure refl_name exists.
            try:
                ncid.variables[refl_name]
            except Exception,e:
                print >>sys.stderr, 'field name', e, 'does not exist'
                ncid.close()
                continue

            # Read in variables
            refl = np.array(np.squeeze(ncid.variables[refl_name][:,refl_level-1,:,:]))
            missing_value = ncid.variables[refl_name]._FillValue
            time = ncid.variables['time'][:]
            lat = ncid.variables[lat_name][:]
            lon = ncid.variables[lon_name][:]
            alt = ncid.variables[alt_name][:]

            # Create 2d lat and lon variables
            (lon2d,lat2d) = np.meshgrid(lon,lat)
            #if debug:
            #    print >>sys.stderr, 'lons = ', lon2d[0,:]
            #    print >>sys.stderr, 'lats = ', lat2d[:,0]

            # Turn refl missing_vals into NaN's
            #refl[refl==missing_value] = np.nan
            
            # Close file
            ncid.close()

            # Determine raintype
            (rain_type,conv_cores,bkgrnd) = conv_strat_latlon(refl, missing_value,lat2d,
                                                              lon2d, core_thresh,
                                                              method, a, b, tune_thresh,
                                                              sm_rad, fill_dbz, bg_diff)
            rain_type[np.isnan(rain_type)] = missing_value
            conv_cores[conv_cores==-1] = missing_value
            bkgrnd[np.isnan(bkgrnd)] = missing_value

            # Simple raintype
            rain_type_basic = np.copy(rain_type)
            rain_type_basic[rain_type==0] = 1
            rain_type_basic[rain_type>0] = 2

            # Output netcdf file
            writeCFnetcdf(fname_out,core_thresh,method,a,b,tune_thresh,sm_rad,
                          fill_dbz,bg_diff,institution,source,title,
                          references,comment,time,lat,lon,lat_origin,lon_origin,
                          rain_type,rain_type_basic,conv_cores,bkgrnd,missing_value)

