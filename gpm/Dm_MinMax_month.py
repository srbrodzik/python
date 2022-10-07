#!/usr/bin/python3

import os
import netCDF4 as nc
import numpy as np

indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v07/2014/07'
max_month = -99.
min_month = 99.

fout = open(indir+'/dm_info','w')

for file in os.listdir(indir):
    if file.endswith('nc'):
        print(file)
        fout.write('file = {}\t'.format(file))
        ncid = nc.Dataset(indir+'/'+file,'r')
        dm_id = ncid.variables['dsd_dm']
        dm = dm_id[:]
        ncid.close()
        max = np.max(dm)
        min = np.min(dm)
        fout.write('min = {:.2f}\tmax = {:.2f}\n'.format(min,max))
        if min < min_month:
            min_month = min
        if max > max_month:
            max_month = max

fout.write('min_month = {:.2f} and max_month = {:.2f}\n'.format(min_month,max_month))

fout.close()

        
        
