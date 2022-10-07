import os
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma

indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/06'

for fname in os.listdir(indir):
    if fname.endswith('nc'):
        ncid = nc4.Dataset(indir+'/'+fname,'r')
        srt = ncid.variables['shallow_rain_type'][:]
        vals = np.unique(srt)
        if vals.size > 3:
            print fname
            print vals.data
        ncid.close()
