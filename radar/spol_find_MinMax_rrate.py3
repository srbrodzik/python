#!/usr/bin/python3

# Imports and functions
import os
import netCDF4 as nc4
import numpy as np

# User inputs
indir = '/home/disk/monsoon/precip/cfradial/spol_ncar/cfradial/rate/sur/20220715'

minRrateHybrid = 99
maxRrateHybrid = -99
minRrateZH = 99
maxRrateZH = -99

for fname in os.listdir(indir):
    
    if fname.endswith('nc'):

        print(fname)

        # open and read lat and lon from input file
        ncid = nc4.Dataset(indir+'/'+fname,'r')
        rr_hybrid = np.array(ncid.variables['RATE_HYBRID'])
        rr_zh = np.array(ncid.variables['RATE_ZH'])
        ncid.close()

        rr_hybrid[rr_hybrid == -9999.] = np.nan
        rr_zh[rr_zh == -9999.] = np.nan
        
        # find min and max lat and lon
        minRRhybrid = np.nanmin(rr_hybrid)
        minRRzh = np.nanmin(rr_zh)
        maxRRhybrid = np.nanmax(rr_hybrid)
        maxRRzh = np.nanmax(rr_zh)

        print('   min/max rr_hybrid = ',minRRhybrid,'/',maxRRhybrid)
        print('   min/max rr_zh     = ',minRRzh,'/',maxRRzh)
                
        if minRRhybrid < minRrateHybrid:
            minRrateHybrid = minRRhybrid
        if minRRzh < minRrateZH:
            minRrateZH = minRRzh
        if maxRRhybrid > maxRrateHybrid:
            maxRrateHybrid = maxRRhybrid
        if maxRRzh > maxRrateZH:
            maxRrateZH = maxRRzh

print('min/max rr_hybrid = ',minRrateHybrid,'/',maxRrateHybrid)
print('min/max rr_zh     = ',minRrateZH,'/',maxRrateZH)

