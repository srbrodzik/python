from netCDF4 import Dataset
#import netCDF4
#from scipy.io import netcdf
import numpy as np

def read_nc_xyzt_file_func(ncid,):

    # get dim info:
    dimDict = {}
    for dim in ncid.dimensions:
        dimDict.update({dim:ncid.dimensions[dim].size})

    # get var dim and attr info
    varDict = {}
    for var in ncid.variables:
        print(var)
        print(ncid.variables[var].dimensions)
        varDict.update({var:ncid.variables[var].dimensions})
    
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

