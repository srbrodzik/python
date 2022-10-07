import netCDF4 as nc4
import numpy as np
import os
import re

indir = "/home/disk/bob/gpm/sam_ku/classify/ex_data/2014/03"
outFile = '/home/disk/bob/gpm/sam_ku/classify/ex_data/2014/201403_minLat.txt'

totalMinsLT65 = 0

fid = open(outFile,'w')

for file in os.listdir(indir):
    if file.endswith('nc'):
        print(file)
        ncid = nc4.Dataset(indir+'/'+file,'r')
        lats = ncid.variables['lat'][:]
        ncid.close()
        minLat = lats[0]
        if minLat < -60.0:
            totalMinsLT65 = totalMinsLT65 + 1
            fid.write("{} {}\n".format(file,minLat))

fid.write("{} {}\n".format("Total min lats < 65S = ",totalMinsLT65))
fid.close()
