#!/usr/bin/python3

# Read Nw and Dm values from 2AKu HDF5 files

import os
import h5py
import numpy as np

indir = '/home/disk/bob/gpm/hdf_subsetted/2Ku/2021/01_SP'

outfile = 'paramDSD_values.txt'
fout = open(indir+'/'+outfile,'w')

for file in os.listdir(indir):

    if file.startswith('2A-SP'):
    
        print(file)
        fout.write('{}\n'.format(file))

        f = h5py.File(file,"r")

        #for group in f.keys():
        #    if group == 'NS':
        #        print('NS group')
        #        for subgroup in f[group].keys():
        #            if subgroup == 'SLV':
        #                print('SLV subgroup')
        #                for dset in f[group][subgroup].keys():
        #                    print(dset)
        
        dsd = f['NS']['SLV']['paramDSD']
        FillValue = dsd.attrs['_FillValue']
        param1 = np.array(dsd[:,:,:,0])
        param2 = np.array(dsd[:,:,:,1])
        fout.write('  param1 min/max values = {min:6.2f}/{max:6.2f}\n'.format(min=np.min(param1[param1 != FillValue]),
                                                                              max=np.max(param1[param1 != FillValue])))
        #fout.write('  % non-missing param1 values = ',
        #           np.count_nonzero(param1 != FillValue)/dsd.size*100)
        fout.write('  param2 min/max values = {min:6.2f}/{max:6.2f}\n'.format(min=np.min(param2[param2 != FillValue]),
                                                                              max=np.max(param2[param2 != FillValue])))
        #fout.write('  num non-missing param2 values = ',
        #           np.count_nonzero(param2 != FillValue)/dsd.size*100)
        
        f.close()

fout.close()
