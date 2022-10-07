from netCDF4 import Dataset
#import netCDF4
#from scipy.io import netcdf
import numpy as np

# open a netcdf file for reading
ncid = Dataset('/home/disk/mjo/dynamo/data.server/zebra/moments/sband/sur_2km/20111016/spolka.20111016.120032.nc','r')
#ncid = Dataset('simple_xyzt.nc','r')
#ncid = netCDF4.Dataset('simple_xyzt.nc','r')
#ncid = netcdf.netcdf_file('simple_xyzt.nc','r')

# to get dim info:
# print ncid.dimensions.keys()
# ncid.dimensions['time']
# ncid.dimensions['z']

# to get dim info about a var:
# ncid.variables['refl'].dimensions

# to get a dictionary of variables:
# for v in ncid.variables:
#       print(v)
# OR
# print(ncid)
# which gives dim and var info
# OR
# ncid.variables.keys()

# print variable metadata from netcdf file
print ncid.variables['DBZ']
# read variable from netcdf file
refl = ncid.variables['DBZ'][:]
#refl = ncid.variables['refl'][:]
nt,nz,ny,nx = refl.shape

# get individual global attribute info
# print ncid.Conventions
# OR
# find all global attributes
# for attr in ncid.ncattrs():
#     print attr, '=', getattr(ncid,attr)

# get variable attributes
# print refl.units
# OR
# find all var attributes
# for attr in 

# close the netcdf file
ncid.close()







# read dimensions
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
