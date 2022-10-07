import netCDF4 as nc4
import numpy as np
import os
import re

indir = "/home/disk/bob/gpm/asia_ku/classify/ex_data/2014/03"
outFile = '/home/disk/bob/gpm/asia_ku/classify/ex_data/2014/201403_maxLat.txt'

totalMaxGT65 = 0

fid = open(outFile,'w')

for file in os.listdir(indir):
    if file.endswith('nc'):
        print(file)
        ncid = nc4.Dataset(indir+'/'+file,'r')
        lats = ncid.variables['lat'][:]
        ncid.close()
        maxLat = np.amax(lats)
        if maxLat > 65.0:
            totalMaxGT65 = totalMaxGT65 + 1
            fid.write("{} {}\n".format(file,maxLat))

fid.write("{} {}\n".format("Total max lats > 65N = ",totalMaxGT65))
fid.close()
