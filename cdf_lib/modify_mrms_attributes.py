#!/usr/bin/python3

import os
import sys
import netCDF4 as nc

if len(sys.argv) != 2:
    print('Usage: {} [kdp/rhohv/sw/zdr]'.format(sys.argv[0]))
    sys.exit()
else:
    var = sys.argv[1]

if var not in ['kdp','rhohv','sw','zdr']:
    print('Input param must be either kdp, rhohv, sw or zdr')
    sys.exit()
    
inDirBase = '/home/disk/bob/impacts/mdv/mrms'
varDict = {'kdp':{'subdir':'3DKdp','units':'deg/km','long_name':'specific_differential_phase'},
           'rhohv':{'subdir':'3DRhoHV','units':'-','long_name':'correlation_coefficient'},
           'sw':{'subdir':'3DSpw','units':'m/s','long_name':'spectrum_width'},
           'zdr':{'subdir':'3DZdr','units':'dB','long_name':'differential_reflectivity'}}

inDirVar = inDirBase+'/'+varDict[var]['subdir']
for date in os.listdir(inDirVar):
    if date.startswith('2022'):
        print('date = {}'.format(date))
        inDir = inDirVar+'/'+date
        for file in os.listdir(inDir):
            if file.endswith('.nc'):
                print('file = {}'.format(file))
                ncid = nc.Dataset(inDir+'/'+file,'a')
                ncid.variables[var.upper()].long_name = varDict[var]['long_name']
                ncid.variables[var.upper()].units = varDict[var]['units']
                ncid.close()
                
