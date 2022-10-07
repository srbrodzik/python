import netCDF4 as nc4
import numpy as np
import time as tm
import datetime as dt

refl_missing = -999.0
precip_missing = -99.0

# open a compressed netcdf4 file for reading
ncid = nc4.Dataset('/home/disk/mjo/nmq/data.server/raw/nexrad/20110521/regrid_nmq_20110521-000000.nc','r')

# read vars
lats2d = ncid.variables['latitude'][:]
lats = lats2d[:,0]
nlats = lats.size
lons2d = ncid.variables['longitude'][:]
lons = lons2d[0,:]
nlons = lons.size
hts = ncid.variables['height'][:]
nhts = hts.size
ntimes = 1
refl_mean = ncid.variables['mrefl_mean'][:]
precip = ncid.variables['prec_hsr'][:]

# read global attributes
title = ncid.title
time_in_secs = ncid.Time

# close input file
ncid.close()

# create new netcdf file with renamed dims, resized lat and lon and a new time variable
ncname = '/home/disk/mjo/nmq/data.server/mdv/cart/nex_mosaic/20110521/regrid_nmq_20110521-000000_new.nc'
ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

# create dims
time = ncid.createDimension('time',None) # None implies UNLIMITED
z0 = ncid.createDimension('z0',nhts)
y0 = ncid.createDimension('y0',nlats)
x0 = ncid.createDimension('x0',nlons)

# create vars
time_var = ncid.createVariable('time',np.float64,('time'),zlib=True )
z0_var = ncid.createVariable('z0',np.float32,('z0'),zlib=True )
y0_var = ncid.createVariable('y0',np.float32,('y0'),zlib=True )
x0_var = ncid.createVariable('x0',np.float32,('x0'),zlib=True )
gm_var = ncid.createVariable('grid_mapping_0',np.int32)
refl_var = ncid.createVariable('mrefl_mean',np.float32,('time','z0','y0','x0'),
                               fill_value=refl_missing,zlib=True )
precip_var = ncid.createVariable('prec_hsr',np.float32,('time','y0','x0'),
                                 fill_value=precip_missing,zlib=True )

# define attributes
time_var.standard_name = 'time'
time_var.long_name = 'Data time'
time_var.units = 'seconds since 1970-01-01T00:00:00Z'
time_var.calendar = 'standard'
time_var.axis = 'T'
#time_var.comment = datetime

z0_var.standard_name = 'altitude'
z0_var.long_name = 'altitude MSL'
z0_var.units = 'km'
z0_var.positive = 'up'
z0_var.axis = 'Z'

y0_var.standard_name = 'latitude'
y0_var.long_name = 'latitude'
y0_var.units = 'degrees_north'
y0_var.axis = 'Y'

x0_var.standard_name = 'longitude'
x0_var.long_name = 'longitude'
x0_var.units = 'degrees_east'
x0_var.axis = 'X'

gm_var.grid_mapping_name = 'latitude_longitude'

refl_var.long_name = 'Regridded mosaic reflectivity (mean)'
refl_var.units = 'dBZ'

precip_var.long_name = 'preciprate_hsr'
precip_var.units = 'mm/hr'

# define global attributes
ncid.title = title
ncid.Time = time_in_secs

# write vars to file
time_var[:] = time_in_secs
z0_var[:] = hts
y0_var[:] = lats
x0_var[:] = lons
refl_var[:,:,:,:] = refl_mean
precip_var[:,:,:] = precip

# close file
ncid.close()


