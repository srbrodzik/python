#!/usr/bin/python3

"""
Extracts reflectivity CFAD output from allStorms_xxx_v12x_xxx.pro saved in 
 .../classify/class_data/stats_class_v12x/Cfad/infoCfad_*.zip
netcdf files and create new netcdf4 files.

Modified from GPM code
Stacy Brodzik, 13 Oct 2022
Univ of Washington, Dept of Atmos Sciences
"""

import sys
import os
import glob
from zipfile import ZipFile
from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from datetime import datetime
from datetime import timedelta

if len(sys.argv) != 5:
    errMsg = 'Usage: '+sys.argv[0]+' region[AFC|CIO|EPO|H01|H02|H05|NAM|SAM|SAS|WMP] year[YYYY] month[MM] thresLevel[str|mod]'
    sys.exit(errMsg)
else:
    region = sys.argv[1]
    iyear = sys.argv[2]
    imonth = sys.argv[3]
    thresLevel = sys.argv[4]

verbose = True
codeVersion = 'TRMM2APR7_uw4'
trmmVersion = 'v07'
if thresLevel == 'str':
    features = {'Conv':'conv',
                'ShallowIsol': 'shi',
                'Stra':'stra'}
elif thresLevel == 'mod':
    features = {'Conv':'conv',
                'Stra':'stra'} 
data_type = 'grd'
path_tmp='/virtual'
convCoreTypes = ['DCC','WCC','DWC']
types_out = np.array(convCoreTypes, dtype='object')

missing_flt = -9999.
missing_int = -9999

regionDict = {'AFC':{'subdir':'afc_pr','long_name':'Africa','limits':[-40.,-30.,40.,60.]},
              'CIO':{'subdir':'cio_pr','long_name':'Central Indian Ocean','limits':[-40.,55.,10.,110.]},           
              'EPO':{'subdir':'epo_pr','long_name':'East Pacific Ocean','limits':[-40.,-178.,40.,-130.]},
              'H01':{'subdir':'h01_pr','long_name':'Hole 01 (West Pacific Ocean)','limits':[-40.,-140.,25.,-85.]},
              'H02':{'subdir':'h02_pr','long_name':'Hole 02 (North Atlantic Ocean)','limits':[15.,-65.,40.,-10.]},
              'H05':{'subdir':'h05_pr','long_name':'Hole 05 (Western Pacific)','limits':[5.,125.,40.,178.]},
              'NAM':{'subdir':'nam_pr','long_name':'North America','limits':[15.,-140.,40.,-55.]},
              'SAM':{'subdir':'sam_pr','long_name':'South America','limits':[-40.,-95.,20.,-25.]},
              'SAS':{'subdir':'sas_pr','long_name':'South Asia','limits':[5.,55.,40.,130.]},
              'WMP':{'subdir':'wmp_pr','long_name':'Warm Pool','limits':[-40.,105.,10.,178.]}}

if thresLevel == 'str':
    inDir = '/home/disk/bob/trmm/'+regionDict[region]['subdir']+'/classify/class_data_v07/stats_class_v12s/Cfad'
elif thresLevel == 'mod':
    inDir = '/home/disk/bob/trmm/'+regionDict[region]['subdir']+'/classify/class_data_v07/stats_class_v12m/Cfad'

outDirBase = '/home/disk/archive3/trmm/'+trmmVersion+'/'+region+'/cfad'
if not os.path.isdir(outDirBase):
    os.makedirs(outDirBase)

os.chdir(inDir+'/'+imonth)

# Get time in seconds since 1-Jan-1970
time_str = iyear+imonth
time_obj = datetime.strptime(time_str,'%Y%m')
secs_since_1970 = int( (time_obj - datetime(1970,1,1)).total_seconds() )

for ifeature in features.keys():
        
    for file in glob.glob('infoCfad_'+ifeature+'_EchoCores_*_'+iyear+'_*.zip'):

        print('file = ',file)

        prefix = os.path.splitext(file)[0]
        #strParts = prefix.split('_')
        [junk1,feature_long,junk2,month,year,junk3,version] = prefix.split('_')

        # Open original netcdf file and read var info & global atts
        with ZipFile(file, 'r') as zipObj:
            zipObj.extractall(path_tmp)
            
        cdfFile = path_tmp+'/'+prefix+'.nc'
        in_id =  Dataset(cdfFile,'r')

        cfad_full = in_id.variables['CFAD_Full'][:]
        cfad_core = in_id.variables['CFAD_Core'][:]
        refl      = in_id.variables['refl'][:]
        numRefls  = len(refl)
        alt       = in_id.variables['altitude'][:]
        numAlts   = len(alt)

        title     = in_id.getncattr('Title')
        source    = in_id.getncattr('source')
        if ifeature == 'Conv':
            etype = in_id.getncattr('echo_type')

        in_id.close()
        os.remove(cdfFile)

        if verbose:
            print('Done reading input file.')
            print('Creating new nc file')

        # Create new netcdf file
        outDir = outDirBase+'/'+features[ifeature]+'/'+month
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        cdfFile = outDir+'/'+codeVersion+'_reflCfad_'+features[ifeature]+'_'+thresLevel+'_'+data_type+'_'+year+month+'_'+region+'.nc'
        id = Dataset(cdfFile,'w',format='NETCDF4')

        # Make dimensions
        id.createDimension('time',None)
        if ifeature == 'Conv':
            id.createDimension('type_dim',len(convCoreTypes))
            #id.createDimension('echo_type',len(convCoreTypes))
        id.createDimension('alt_bin',numAlts)
        id.createDimension('refl_bin',numRefls)
        if verbose:
            print('Done with dims')

        # Define variables:
        time_var_id = id.createVariable('time','i4',('time',),zlib=True,complevel=5)
        time_var_id.units = 'seconds'
        time_var_id.long_name = 'seconds since 1970-01-01'

        if ifeature == 'Conv':
            #type_var_id = id.createVariable('echo_type','i4',('echo_type',))
            #type_var_id.units = 'none'
            #type_var_id.long_name = 'convective core type'
            type_var_id = id.createVariable('echo_type','str',('type_dim',))
            type_var_id.units = 'none'
            type_var_id.long_name = 'convective core type'
            
        alt_var_id = id.createVariable('altitude',np.float32,('alt_bin',),zlib=True,complevel=5)
        alt_var_id.units = 'km'
        alt_var_id.long_name = 'altitude bin value'

        refl_var_id = id.createVariable('refl','f4',('refl_bin',),zlib=True,complevel=5)
        refl_var_id.units = 'dBZ'
        refl_var_id.long_name = 'reflectivity bin value'

        full_var_id = id.createVariable('cfad_full','i4',('time','alt_bin','refl_bin',),
                                                fill_value=missing_int,zlib=True,complevel=5)
        full_var_id.units = 'none'
        full_var_id.long_name = 'refl counts for full storm'

        if ifeature == 'Conv':
            #core_var_id = id.createVariable('cfad_core','i4',('time','echo_type','alt_bin','refl_bin',),
            core_var_id = id.createVariable('cfad_core','i4',('time','type_dim','alt_bin','refl_bin',),
                                            fill_value=missing_int,zlib=True,complevel=5)
            core_var_id.units = 'none'
            core_var_id.long_name = 'refl counts for core regions only'
        else:
            core_var_id = id.createVariable('cfad_core','i4',('time','alt_bin','refl_bin',),
                                            fill_value=missing_int,zlib=True,complevel=5)
            core_var_id.units = 'none'
            core_var_id.long_name = 'refl counts for core regions only'
                    
        if verbose:
            print('Done with vars')
            
        # define global variables
        id.Title = title
        id.source = source
        if ifeature == 'Conv':
            id.echo_type = etype

        if verbose:
            print('Done with global attr')

        # write data to file
        time_var_id[:] = np.asarray(secs_since_1970)
        alt_var_id[:]  = alt
        refl_var_id[:] = refl
        if ifeature == 'Conv':
            #type_var_id[:] = np.asarray([0,1,2])
            type_var_id[:] = types_out
        full_var_id[:] = np.reshape(cfad_full,(1,numAlts,numRefls))
        if ifeature == 'Conv':
            core_var_id[:] = np.reshape(cfad_core,(1,len(convCoreTypes),numAlts,numRefls))
        else:
            core_var_id[:] = np.reshape(cfad_core,(1,numAlts,numRefls))

        if verbose:
            print('Done writing data')

        # close netcdf file
        id.close()
               

                
                
