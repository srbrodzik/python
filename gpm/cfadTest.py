#!/usr/bin/python3

import os
import netCDF4 as nc

indir = '/home/disk/bob/gpm/nam_ku/classify/class_data_v07/stats_class_v12s/Cfad_Dm/06'
#infile = 'infoCfad_Conv_EchoCores_06_2014_NAM_v12s.nc'
infile = 'infoCfad_Dm_Conv_EchoCores_06_2014_NAM_v12s.nc'
typeDict = {0:'DCC',
            1:'WCC',
            2:'DWC'}

# read netcdf file
ncid = nc.Dataset(indir+'/'+infile,'r')
full_id = ncid.variables['Dm_CFAD_Full']
full = full_id[:][:]
fullCnts = full[:,:]
core_id = ncid.variables['Dm_CFAD_Core']
core = core_id[:][:][:]
ncid.close

# compare DC, WC, DWC counts to full count
(types,alts,refls) = core.shape
for type in range(0,types):
    print('cfad type = {}'.format(typeDict[type]))
    coreCnts = core[type,:,:]
    for alt in range(0,alts):
        for refl in range(0,refls):
            if coreCnts[alt,refl] > fullCnts[alt,refl]:
                print('at alt = {} and refl = {}:'.format(alt,refl))
                print('   coreCnts = {} > fullCnts = {}'.format(coreCnts[alt,refl],fullCnts[alt,refl]))
                
"""
# compare total of DC, WC, DWC counts to full count
(types,alts,refls) = core.shape
for alt in range(0,alts):
    for refl in range(0,refls):
        if sum(core[0:3,alt,refl]) > fullCnts[alt,refl]:
            print('at alt = {} and refl = {}:'.format(alt,refl))
            print('   sum(core) = {} > fullCnts = {}'.format(sum(core[0:3,alt,refl]),fullCnts[alt,refl]))
"""                
   
    
    
