#!/usr/bin/python3

"""
Converted from IDL code 
/home/disk/meso-home/brodzik/IDL/gpm/for_web_site/v06/convert_info_rain_nsr_bsr_grd_cdf_to_cdf.pro

Extracts stratiform rainrate sum output from allStorms_xxx_v12x_xxx.pro saved in 
 .../classify/class_data/stats_class_v12x/NearSurfRain/MM/infoRainNSR_Stra*zip
netcdf files and create new netcdf4 files.  The original files include both coarse and 
fine grids.  We only want to put fine grid rainrate sums on public website.

Stacy Brodzik, 01 Dec 2020 (converted to python)
Univ of Washington, Dept of Atmos Sciences

Stacy Brodzik, 11 Apr 2022 (converted v06 version for v07)

Stacy Brodzik, 13 Oct 2022 (modified GPM v07 code)
"""

import sys
import os
import glob
from zipfile import ZipFile
from netCDF4 import Dataset
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

print('region = ',region)
print('iyear = ',iyear)
print('imonth = ',imonth)
print('thresLevel = ',thresLevel)

verbose = False
codeVersion = 'TRMM2APR7_uw4'
trmmVersion = 'v07'
feature = 'BSR'
data_type = 'grd'
path_tmp='/virtual'

missing_flt = -9999.
missing_int = -9999

if region == 'AFC':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/afc_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/afc_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'Africa'
    limits = [-40.,-30.,40.,60.]   # [minLat,minLon,maxLat,maxLon]
elif region == 'CIO':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/cio_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/cio_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'Central Indian Ocean'
    limits = [-40.,55.,10.,110.]   # [minLat,minLon,maxLat,maxLon]
elif region == 'EPO':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/epo_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/epo_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'East Pacific Ocean'
    limits = [-40.,-178.,40.,-130.]# [minLat,minLon,maxLat,maxLon]
elif region == 'H01':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/h01_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/h01_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain' 
    region_long = 'Hole 01 (West Pacific Ocean)'
    limits = [-40.,-140.,25.,-85.] # [minLat,minLon,maxLat,maxLon]
elif region == 'H02':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/h02_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/h02_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'Hole 02 (North Atlantic Ocean)'
    limits = [15.,-65.,40.,-10.]   # [minLat,minLon,maxLat,maxLon]
elif region == 'H05':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/h05_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/h05_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'Hole 05 (Western Pacific)'
    limits = [5.,125.,40.,178.]  # [minLat,minLon,maxLat,maxLon]
elif region == 'NAM':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/nam_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/nam_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'North America'
    limits = [15.,-140.,40.,-55.]  # [minLat,minLon,maxLat,maxLon]
elif region == 'SAM':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/sam_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/sam_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'South America'
    limits = [-40.,-95.,20.,-25.]  # [minLat,minLon,maxLat,maxLon]
elif region == 'SAS':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/sas_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/sas_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'South Asia'
    limits = [5.,55.,40.,130.]     # [minLat,minLon,maxLat,maxLon]
elif region == 'WMP':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/wmp_pr/classify/class_data_v07/stats_class_v12s/NearSurfRain'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/wmp_pr/classify/class_data_v07/stats_class_v12m/NearSurfRain'
    region_long = 'Warm Pool'
    limits = [-40.,105.,10.,178.]  # [minLat,minLon,maxLat,maxLon] 

outDirBase = '/home/disk/archive3/trmm/'+trmmVersion+'/'+region 
#outDirBase = '/home/disk/archive3/trmm/'+trmmVersion+'/'+region+'_test'

if verbose:
    print('inDir = ',inDir)
    print('outDirBase = ',outDirBase)
    print('region = ',region)
    print('region_long = ',region_long)
    print('limits = ',limits) 

os.chdir(inDir+'/'+imonth)

# Get time in seconds since 1-Jan-1970
time_str = iyear+imonth
time_obj = datetime.strptime(time_str,'%Y%m')
secs_since_1970 = int( (time_obj - datetime(1970,1,1)).total_seconds() )
                
for file in glob.glob('infoRainNSR_Stra_EchoCores_*_'+iyear+'_*.zip'):

    print('file = ',file)

    prefix = os.path.splitext(file)[0]
    strParts = prefix.split('_')
    [junk1,feature_long,junk2,month,year,junk3,version] = prefix.split('_')

    # Open original netcdf file and read var info
    with ZipFile(file, 'r') as zipObj:
        zipObj.extractall(path_tmp)
            
    cdfFile = path_tmp+'/'+prefix+'.nc'
    in_id =  Dataset(cdfFile,'r')

    core_rain_rate_sum    = in_id.variables['rain_Core_BS'][:]
    core_rain_rate_count  = in_id.variables['nRain_Core_BS'][:]
    storm_rain_rate_sum   = in_id.variables['rain_Full_BS'][:]
    storm_rain_rate_count = in_id.variables['nRain_Full_BS'][:]
    lon = in_id.variables['lonsF'][:]
    numLons = len(lon)
    lat = in_id.variables['latsF'][:]
    numLats = len(lat)
        
    in_id.close()
    os.remove(cdfFile)

    if verbose:
        print('Done reading input file.')
        print('Creating new nc file')
        
    # Create new netcdf file
    outDir = outDirBase+'/'+feature+'/'+month
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    cdfFile = outDir+'/'+codeVersion+'_'+feature+'_'+thresLevel+'_'+data_type+'_'+year+month+'_'+region+'.nc'
    id = Dataset(cdfFile,'w',format='NETCDF4')

    # Make dimensions
    id.createDimension('time',None)
    id.createDimension('longitude',numLons)
    id.createDimension('latitude',numLats)
    if verbose:
        print('Done with dims')

    # Define variables:
    time_var_id = id.createVariable('time','i4',('time',),zlib=True,complevel=5)
    time_var_id.units = 'seconds'
    time_var_id.long_name = 'seconds since 1970-01-01'
            
    lon_var_id = id.createVariable('longitude',np.float32,('longitude',),zlib=True,complevel=5)
    lon_var_id.units = 'degrees'
    lon_var_id.long_name = 'longitude'

    lat_var_id = id.createVariable('latitude','f4',('latitude',),zlib=True,complevel=5)
    lat_var_id.units = 'degrees'
    lat_var_id.long_name = 'latitude'

    crr_var_id = id.createVariable('bsr_core_rain_rate_sum','f4',('time','latitude','longitude',),
                                   fill_value=missing_flt,zlib=True,complevel=5)
    crr_var_id.units = 'mm/hr'
    crr_var_id.long_name = 'broad stratiform core rainrate sum based on TRMM PR /FS/SLV/precipRateNearSurface'

    cct_var_id = id.createVariable('bsr_core_rain_rate_count','i4',('time','latitude','longitude',),
                                   fill_value=missing_int,zlib=True,complevel=5)
    cct_var_id.units = 'none'
    cct_var_id.long_name = 'number of broad stratiform core pixels with TRMM PR /FS/SLV/precipRateNearSurface greater than zero'

    srr_var_id = id.createVariable('bsr_storm_rain_rate_sum','f4',('time','latitude','longitude',),
                                   fill_value=missing_flt,zlib=True,complevel=5)
    srr_var_id.units = 'mm/hr'
    srr_var_id.long_name = 'broad stratiform storm rainrate sum based on TRMM PR /FS/SLV/precipRateNearSurface'

    sct_var_id = id.createVariable('bsr_storm_rain_rate_count','i4',('time','latitude','longitude',),
                                   fill_value=missing_int,zlib=True,complevel=5)
    sct_var_id.units = 'none'
    sct_var_id.long_name = 'number of broad stratiform storm pixels with TRMM PR /FS/SLV/precipRateNearSurface greater than zero'

    if verbose:
        print('Done with vars')
            
    # define global variables
    id.Title = 'broad stratiform monthly rain rate sums based on TRMM PR /FS/SLV/precipRateNearSurface and /FS/CSF/typePrecip'
    if thresLevel == 'str':
        id.BroadStratiform_Criterion = 'contiguous stratiform echo >= 50,000km^2'
    elif thresLevel == 'mod':
        id.BroadStratiform_Criterion = 'contiguous stratiform echo >= 40,000km^2'
    id.source_uw = 'Created by the Mesoscale Group, University of Washington'
    id.source_nasa = 'Original data obtained from NASA Goddard Earth Sciences http://pmm.nasa.gov/'
    id.reference = 'Methodology described in Houze et al. 2015 (Reviews of Geophysics)'
    id.data_location = 'http://trmm.atmos.washington.edu'
    id.region = region+' '+region_long
    id.lon_min = limits[1]
    id.lon_max = limits[3]
    id.lat_min = limits[0]
    id.lat_max = limits[2]
    id.time_range = 'Month starting at date defined by time variable'

    if verbose:
        print('Done with global attr')

    # write the data into file
    time_var_id[:] = np.asarray(secs_since_1970)
    lon_var_id[:] = lon
    lat_var_id[:] = lat
    crr_var_id[:] = np.reshape(core_rain_rate_sum,(1,numLats,numLons))
    cct_var_id[:] = np.reshape(core_rain_rate_count,(1,numLats,numLons))
    srr_var_id[:] = np.reshape(storm_rain_rate_sum,(1,numLats,numLons))
    sct_var_id[:] = np.reshape(storm_rain_rate_count,(1,numLats,numLons))

    if verbose:
        print('Done writing data')

    # close netcdf file
    id.close()
