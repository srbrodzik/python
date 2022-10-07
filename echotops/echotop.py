# Echo top heights
import numpy as np
import netCDF4 as nc4
import os

#infile='FILE_NAME'
infile='/data/hein/new/dynamo/20111024/RevelleTOGA.20111024_230001.nc'
outfile='echotop_'+os.path.basename(infile)
nx=300
ny=300
nz=40
dz=0.5

# read in reflectivity data
ncidin = nc4.Dataset(infile,'r')
refl = np.array(np.squeeze(ncidin.variables['REFL']))
ncidin.close()
print refl.shape
refl_level = [0,10,20,30,40,50]
nmax = len(refl_level)
echotop = np.zeros((nmax,ny,nx))-1.

# Find echo tops
for ii in range(nmax):
    maxht = np.zeros((ny,nx))-1.
    for nl in range(nz-1):
        #print nl, refl[nl,:,:].shape
        maxht = np.where(np.logical_and((refl[nl,:,:] >= refl_level[ii]),(refl[nl+1,:,:] < refl_level[ii])),nl*dz, maxht)
        #print maxht
    echotop[ii,:,:] = maxht
    print np.max(maxht)

    print outfile
# Write echo top values to a file
    ncidout = nc4.Dataset(outfile,'w',format='NETCDF4')
    # create dimensions
    refl_levels_dim = ncidout.createDimension('refl_levels_dim',nmax) 
    y = ncidout.createDimension('y',ny)
    x = ncidout.createDimension('x',nx)

    # create variables
    refl_levels = ncidout.createVariable('refl_level',np.float32,('refl_levels_dim'),zlib=True )
    etop = ncidout.createVariable('echo_top',np.float32,('refl_levels_dim','y','x'),zlib=True )

    # create variable attributes
    # refl_levels
    refl_levels.long_name = 'reflectivity_levels'
    refl_levels.units = 'dBZ'
    # etop
    etop.long_name = 'reflectivity_maximum_height'
    etop.units = 'km'

    # write vars to file
    refl_levels[:] = refl_level
    etop[:,:,:] = echotop

    ncidout.close()
