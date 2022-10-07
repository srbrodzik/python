import numpy as np
import netCDF4 as nc4
import os

infile='/home/disk/bob/gpm/analysis/rainClim/data/2A_GPM_Ku_counts_2016_08_v1_save.nc'

ncidin = nc4.Dataset(infile,'r')
refl = np.array(np.squeeze(ncidin.variables['reflect25']))
ncidin.close()

refl_new = refl*100.

