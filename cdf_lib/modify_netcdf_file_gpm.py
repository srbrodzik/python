import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

# open a compressed netcdf4 file for reading
#ncid = nc4.Dataset('/home/disk/bob/gpm/sam_ku/classify/ex_data/2014/03/GPM2Ku5_uw3_20140319.041420_to_20140319.043609_000303_SAM.nc','r')
#ncid = nc4.Dataset('/home/disk/bob/gpm/nam_ku/classify/ex_data/2015/04/test/GPM-2Ku.006621.20150429.053449.nc','r')
ncid = nc4.Dataset('/home/disk/bob/gpm/sind_ku/classify/ex_data/2016/06/GPM2Ku5_uw3_20160623.163247_to_20160623.163518_013178_SIND.nc','r')

# read vars
times = ncid.variables['time'][:]
ntimes = times.size
lats = ncid.variables['lat'][:]
nlats = lats.size
lons = ncid.variables['lon'][:]
nlons = lons.size
alts = ncid.variables['alt'][:]
nalts = alts.size
nalts_125 = 125
ntimes = 1
time = ncid.variables['time'][:]
refl = ncid.variables['refl'][:]
refl_missing = ncid.variables['refl']._FillValue
rtype = ncid.variables['rain_type'][:]
rtype_missing = ncid.variables['rain_type']._FillValue
#rtype_uw = ncid.variables['rain_type_uw'][:]
#rtype_uw_missing = ncid.variables['rain_type_uw']._FillValue
precip = ncid.variables['near_surf_rain'][:]
precip_missing = ncid.variables['near_surf_rain']._FillValue
swath = ncid.variables['swath'][:]
swath_missing = ncid.variables['swath']._FillValue

# read global attributes
title = ncid.title
orbit = ncid.orbit
lat_min = ncid.lat_min
lat_max = ncid.lat_max
lon_min = ncid.lon_min
lon_max = ncid.lon_max

# close input file
ncid.close()

# create new uncompressed netcdf file with subset of fields for zeb analysis
#ncname = '/home/disk/bob/gpm/sam_ku/classify/ex_data/zebData/GPM2Ku5_uw3_20140319.041420_to_20140319.043609_000303_SAM.nc'
#ncname = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2015/04/test/GPM-2Ku-zeb.006621.20150429.053449.nc'
ncname = '/home/disk/radar/india_main/zeb-india-data/gpm_ku/GPM2Ku5_uw3_20160623.163247_SIND.nc'
ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

# create dims
time = ncid.createDimension('time',None) # None implies UNLIMITED
z0 = ncid.createDimension('altitude',nalts_125)
y0 = ncid.createDimension('latitude',nlats)
x0 = ncid.createDimension('longitude',nlons)

# create vars
time_var = ncid.createVariable('time',np.float64,('time'),zlib=True )
#jday_var = ncid.createVariable('julian_date',np.float64,('time'),zlib=True )
z0_var = ncid.createVariable('altitude',np.float32,('altitude'),zlib=True )
y0_var = ncid.createVariable('latitude',np.float32,('latitude'),zlib=True )
x0_var = ncid.createVariable('longitude',np.float32,('longitude'),zlib=True )
refl_var = ncid.createVariable('refl',np.float32,('time','altitude','latitude','longitude'),
                               fill_value=refl_missing,zlib=True )
rtype_var = ncid.createVariable('rain_type',np.float32,('time','latitude','longitude'),
                               fill_value=rtype_missing,zlib=True )
#rtype_uw_var = ncid.createVariable('rain_type_uw',np.float32,('time','latitude','longitude'),
#                                   fill_value=rtype_uw_missing,zlib=True )
precip_var = ncid.createVariable('near_surf_rain',np.float32,('time','latitude','longitude'),
                                 fill_value=precip_missing,zlib=True )
swath_var = ncid.createVariable('swath',np.float32,('time','latitude','longitude'),
                                fill_value=swath_missing,zlib=True )

# define attributes
time_var.standard_name = 'time'
time_var.long_name = 'Data time'
time_var.units = 'seconds since 1970-01-01T00:00:00Z'
time_var.calendar = 'standard'
time_var.axis = 'T'

#jday_var.long_name = 'Julian Date Number'
#jday_var.units = 'Julian number after October 4, 1582-10-4 12:00:00'

z0_var.units = 'km'
z0_var.positive = 'up'
z0_var.axis = 'Z'

y0_var.standard_name = 'latitude'
y0_var.long_name = 'Latitude coordinate'
y0_var.units = 'degrees_north'
y0_var.axis = 'Y'

x0_var.standard_name = 'longitude'
x0_var.long_name = 'Longitude coordinate'
x0_var.units = 'degrees_west'
x0_var.axis = 'X'

refl_var.long_name = 'GPM-Ku Reflectivity'
refl_var.units = 'dBZ'

rtype_var.long_name = 'Rain Type'
rtype_var.units = 'none'
rtype_var.stratiform = 1
rtype_var.convective = 2
rtype_var.other = 3

#rtype_uw_var.long_name = 'Rain Type UW'
#rtype_uw_var.units = 'none'
#rtype_uw_var.stratiform = 1
#rtype_uw_var.convective = 2
#rtype_uw_var.other = 3

precip_var.long_name = 'Near Surface Rain'
precip_var.units = 'mm/hr'

swath_var.long_name = 'GPM-Ku coverage area'
swath_var.units = 'none'
swath_var.swath = 0
swath_var.no_swath = 1

# define global attributes
ncid.title = title
ncid.orbit = orbit
ncid.lat_min = lat_min
ncid.lat_max = lat_max
ncid.lon_min = lon_min
ncid.lon_max = lon_max

# write vars to file
time_var[:] = times
#jday_var[:] = jdate
z0_var[:] = alts[0:125]
y0_var[:] = lats
x0_var[:] = lons
refl_var[:,:,:,:] = refl[:,0:125,:,:]
rtype_var[:,:,:] = rtype
#rtype_uw_var[:,:,:] = rtype_uw
precip_var[:,:,:] = precip
swath_var[:,:,:] = swath

# close file
ncid.close()


