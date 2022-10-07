import numpy as np
import netCDF4 as nc4
import os

#inDir = '/home/disk/bob/olympex/cfradial/moments/xband_qc/rhi/20151203/output/20151203'
#inFile = 'cfrad.20151203_105607.000_to_20151203_105638.000_CAX1_RHI_RHI.nc'
inDir = '/home/disk/bob/olympex/cfradial/moments/xband_qc/rhi/20151203'
inFile = 'cfrad.20151203_105607.000_to_20151203_105638.000_CAX1_albert_head_v132_RHI.nc'

ncid = nc4.Dataset(inDir+'/'+inFile,'r')
dbz = ncid.variables['DBZ'][:]
scale_factor = ncid.variables['DBZ'].scale_factor
add_offset = ncid.variables['DBZ'].add_offset
ncid.close()

min_dbz = np.amin(dbz)
max_dbz = np.amax(dbz)
