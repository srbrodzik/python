import os
import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt
import pytz

indir  = '/home/disk/radar/india_main/zeb-india-data/faam_trk/raw'
outdir = '/home/disk/radar/india_main/zeb-india-data/faam_trk'

for fname in os.listdir(indir):
    print fname
    if fname.endswith('nc'):

        # Read flight date from fname (FAAM_bxxx_YYYYMMDD.nc)
        tmp = fname.split('_')
        tmp = tmp[2]
        tmp = tmp.split('.')
        date = tmp[0]
        year = date[0:4]
        month = date[4:6]
        day = date[6:8]

        # This assumes t is in local time and converts it to UTC so time is 8 hours too much
        t = dt.datetime(int(year),int(month),int(day),tzinfo=pytz.utc )
        tempBase = int(tm.mktime( t.timetuple() ) )
        eightHrsInSecs = 60 * 60 * 8
        dateInSecsSince1970 = tempBase - eightHrsInSecs
        
        # open netcdf file for reading
        ncid = nc4.Dataset(indir+'/'+fname,'r')

        # read vars
        times = ncid.variables['time'][:]
        ntimes = times.size
        lats = ncid.variables['lat'][:]
        nlats = lats.size
        lons = ncid.variables['lon'][:]
        nlons = lons.size

        # close input file
        ncid.close()

        # create new netcdf file for zeb analysis
        ncname = outdir+'/'+fname
        # NOTES:
        #   NETCDF3_CLASSIC - limited to file sizes less than 2GB
        #   NETCDF3_64BIT_OFFSET - allows for file sizes greater than 2GB
        #   NETCDF3_64BIT_DATA - extends NETCDF3_64BIT_OFFSET to allow for unsigned/64 bit int data types & 64-bit dim sizes
        #   NETCDF4 - default
        ncid = nc4.Dataset(ncname,'w',format='NETCDF3_CLASSIC')

        # create dims
        time = ncid.createDimension('time',None) # None implies UNLIMITED

        # create vars
        base_time_var = ncid.createVariable('base_time',np.int32)
        time_offset_var = ncid.createVariable('time_offset',np.float32,('time') )
        y0_var = ncid.createVariable('latitude',np.float32,('time') )
        x0_var = ncid.createVariable('longitude',np.float32,('time') )

        # define attributes
        base_time_var.standard_name = 'time'
        base_time_var.units = 'seconds since 1970-01-01 00:00:00 +0000'

        time_offset_var.standard_name = 'time offset'
        time_offset_var.units = 'seconds since base_time'

        #y0_var.standard_name = 'latitude'
        y0_var.long_name = 'Latitude coordinate'
        y0_var.units = 'degrees_north'
        #y0_var.axis = 'Y'

        #x0_var.standard_name = 'longitude'
        x0_var.long_name = 'Longitude coordinate'
        x0_var.units = 'degrees_west'
        #x0_var.axis = 'X'

        # define global attributes
        #ncid.title = title
        #ncid.orbit = orbit
        #ncid.lat_min = lat_min
        #ncid.lat_max = lat_max
        #ncid.lon_min = lon_min
        #ncid.lon_max = lon_max

        # write vars to file
        base_time_var[:] = dateInSecsSince1970
        time_offset_var[:] = times
        y0_var[:] = lats
        x0_var[:] = lons

        # close file
        ncid.close()


