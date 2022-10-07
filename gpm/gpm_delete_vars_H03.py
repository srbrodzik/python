# Imports and functions
import os
import netCDF4 as nc4
import numpy as np

# User inputs
indir = '/home/disk/bob/gpm/h03_ku/classify/ex_data_v05'
months = ['01','02','06','07','08','12']
years = ['2015','2016','2017']

for iyear in years:
    for imonth in months:
        indir_full = indir+'/'+iyear+'/'+imonth
        os.chdir(indir_full)
        for fname in os.listdir(indir_full):
            if fname.endswith('nc'):

                print fname

                fname_parts = fname.split('.')
                outfile = fname_parts[0]+'.'+fname_parts[1]+'.'+fname_parts[2]+'.new.nc'

                cmd = 'ncks -a -x -v bsr_mask_str,dcc_mask_str,dwc_mask_str,wcc_mask_str,storm_mask_str,bsr_mask_mod,dcc_mask_mod,dwc_mask_mod,wcc_mask_mod,storm_mask_mod '+fname+' '+outfile
                os.system(cmd)

                cmd = '/bin/rm '+fname
                os.system(cmd)

                cmd = '/bin/mv '+outfile+' '+fname
                os.system(cmd)
                
