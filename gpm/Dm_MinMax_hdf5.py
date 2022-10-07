#!/usr/bin/python3

import os
import h5py
import numpy as np

indir = '/home/disk/bob/gpm/hdf_subsetted/2Ku/2022/01'
Dm_max = -99.
Dm_min = 99.

fout = open(indir+'/dm_info','w')

for file in os.listdir(indir):
    if file.endswith('HDF5'):
        #Dm_max_file = -99.
        #Dm_min_file = 99.
        print(file)
        fout.write('file = {}\t'.format(file))
        f = h5py.File(indir+'/'+file,'r')
        #list(f.keys())
        dset = f['FS']
        #list(dset.keys())
        dgrp = dset['SLV']
        #list(dgrp.keys())
        dvar = dgrp['paramDSD']
        # Can do this in one step:
        dvar = f['FS/SLV/paramDSD']
        
        (nscan,nray,nbin,nvar) = dvar.shape
        #nw = dvar[:,:,:,0]
        dm = dvar[:,:,:,1]
        f.close()

        min = np.min(dm[dm>=0])
        max = np.max(dm)
        fout.write('min = {:.2f}\tmax = {:.2f}\n'.format(min,max))
        
        #if min < Dm_min_file:
        #    Dm_min_file = min
        #if max > Dm_max_file:
        #    Dm_max_file = max

        if min < Dm_min:
            Dm_min = min
        if max > Dm_max:
            Dm_max = max

fout.write('SUMMARY: min = {:.2f}\tmax = {:.2f}\n'.format(Dm_min,Dm_max))


        
