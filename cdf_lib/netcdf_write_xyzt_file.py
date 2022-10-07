#from netCDF4 import Dataset
import netCDF4
#from datetime import datetime, timedelta
#import time
import numpy as np

nz = 2
ny = 10
nx = 10
btVal = np.array(1318728939)
toVal = np.arange(0,1)
latVal = np.array(-1.977621)
lonVal = np.array(71.75462)
altVal = np.array(0.5)
xspVal = np.array(2)
yspVal = np.array(2)
zspVal = np.array(0.5)
reflVal = np.zeros( (1,nz,ny,nx) )

# open a new netcdf file for writing (default is format = 'NETCDF4', not 'NETCDF4_CLASSIC')
#ncid = Dataset('simple_xyzt.nc','w',format='NETCDF4')
ncid = netCDF4.Dataset('simple_xyzt.nc','w',format='NETCDF4')
#print ncid.data_model

# create dimensions
time = ncid.createDimension('time',None) # None implies UNLIMITED
z = ncid.createDimension('z',nz)
y = ncid.createDimension('y',ny)
x = ncid.createDimension('x',nx)
#print ncid.dimensions
#print len(z)
#print len.is_unlimited()  # returns false
#print time.is_unlimited()  # returns true
#for dimobj in ncid.dimensions.values():
#    print dimobj

# create variables
bt = ncid.createVariable('base_time',np.float64 )
to = ncid.createVariable('time_offset',np.float32,('time',) )
lat = ncid.createVariable('lat',np.float32 )
lon = ncid.createVariable('lon',np.float32 )
alt = ncid.createVariable('alt',np.float32 )
xsp = ncid.createVariable('x_spacing',np.float32 )
ysp = ncid.createVariable('y_spacing',np.float32 )
zsp = ncid.createVariable('z_spacing',np.float32 )
refl = ncid.createVariable('refl',np.float32,('time','z','y','x',) )
#print ncid.variables
#print ncid.variables['lat']

# create global attributes
ncid.description = 'test file'
#ncid.history = 'Created ' + time.ctime(time.time())

# create var attributes
bt.units = 'seconds since 1970-01-01 00:00:00 +0000'
to.units = 'seconds since base_time'
lat.units = 'degrees_north'
lon.units = 'degrees_east'
alt.units = 'km'
xsp.units = 'km'
ysp.units = 'km'
zsp.units = 'km'
refl.units = 'dBZ'
refl.missing_value = -9999
#for name in ncid.ncattrs():
#    print 'Global attr', name '=', getattr(ncid,name)
#print ncid.__dict__

# write vars to file
bt[:] = btVal
to[:] = toVal
lat[:] = latVal
lon[:] = lonVal
alt[:] = altVal
xsp[:] = xspVal
ysp[:] = yspVal
zsp[:] = zspVal
refl[:] = reflVal

#close the file
ncid.close()
