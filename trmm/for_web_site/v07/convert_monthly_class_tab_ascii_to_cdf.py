#!/usr/bin/python3

"""
Converted from IDL code 
/home/disk/shear2/brodzik/IDL/gpm/for_web_site/v06/convert_monthly_class_tab_ascii_to_cdf.pro

Converts tabular output from allStorms_xxx_v12x.pro and saved in 
..../classify/class_data/monthly_class_v12s and
..../classify/class_data/monthly_class_v12m 
to netcdf format.

Stacy Brodzik, 01 Dec 2020 (converted to python)
Univ of Washington, Dept of Atmos Sciences

Stacy Brodzik, 11 Apr 2022 (converted v06 version for v07)

Stacy Brodzik, 13 Oct 2022 (modified GPM v07 code)
"""

import sys
import os
import glob
from netCDF4 import Dataset
import numpy as np
import pandas
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

verbose = False
codeVersion = 'TRMM2APR7_uw4'
trmmVersion = 'v07'
data_type = 'tab'
min_parts_in_infile_name = 6

features = {'BroadStratiform':'BSR',
            'DeepWideConvective':'DWC',
            'DeepConvective':'DCC',
            'ShallowIsolated':'SHI',
            'WideConvective':'WCC'}

if region == 'AFC':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/afc_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/afc_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'Africa'
    limits = [-40.,-30.,40.,60.]   # [minLat,minLon,maxLat,maxLon]
elif region == 'CIO':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/cio_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/cio_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'Central Indian Ocean'
    limits = [-40.,55.,10.,110.]   # [minLat,minLon,maxLat,maxLon]
elif region == 'EPO':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/epo_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/epo_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'East Pacific Ocean'
    limits = [-40.,-178.,40.,-130.]# [minLat,minLon,maxLat,maxLon]
elif region == 'H01':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/h01_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/h01_pr/classify/class_data_v07/monthly_class_v12m' 
    region_long = 'Hole 01 (West Pacific Ocean)'
    limits = [-40.,-140.,25.,-85.] # [minLat,minLon,maxLat,maxLon]
elif region == 'H02':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/h02_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/h02_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'Hole 02 (North Atlantic Ocean)'
    limits = [15.,-65.,40.,-10.]   # [minLat,minLon,maxLat,maxLon]
elif region == 'H05':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/h05_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/h05_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'Hole 05 (Western Pacific)'
    limits = [5.,125.,40.,178.]  # [minLat,minLon,maxLat,maxLon]
elif region == 'NAM':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/nam_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/nam_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'North America'
    limits = [15.,-140.,40.,-55.]  # [minLat,minLon,maxLat,maxLon]
elif region == 'SAM':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/sam_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/sam_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'South America'
    limits = [-40.,-95.,20.,-25.]  # [minLat,minLon,maxLat,maxLon]
elif region == 'SAS':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/sas_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/sas_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'South Asia'
    limits = [5.,55.,40.,130.]     # [minLat,minLon,maxLat,maxLon]
elif region == 'WMP':
    if thresLevel == 'str':
        inDir = '/home/disk/bob/trmm/wmp_pr/classify/class_data_v07/monthly_class_v12s'
    elif thresLevel == 'mod':
        inDir = '/home/disk/bob/trmm/wmp_pr/classify/class_data_v07/monthly_class_v12m'
    region_long = 'Warm Pool'
    limits = [-40.,105.,10.,178.]  # [minLat,minLon,maxLat,maxLon] 

outDirBase = '/home/disk/archive3/trmm/'+trmmVersion+'/'+region

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
                
for file_info in glob.glob('*_'+iyear+'_*.info'):

    print('file_info = ',file_info)

    #----------------------------
    # Read and prepare ascii data
    #----------------------------
    # Get matching info and dat filenames
    prefix = os.path.splitext(file_info)[0]
    file_dat = prefix+'.dat'

    # Get feature and thresh
    strParts = prefix.split('_')
    numParts = len(strParts)
    version = strParts[numParts-1]
    #year = strParts[numParts-3]
    #month = strParts[numParts-4]
    feature_long = strParts[0]
    if numParts >= min_parts_in_infile_name:
        for ipart in range(1,numParts-(min_parts_in_infile_name-1)+1):
            feature_long = feature_long+strParts[ipart]
    feature = features[feature_long]
    if feature == 'SHI':
        thresh = 'xxx'
    else:
        thresh = thresLevel

    # Read info file
    f_in = open(file_info,'r')
    orbits = []
    dates = []
    times = []
    core_ids = []
    for line in f_in:
        line = line.strip()
        (orbit,date,time,core_id) = line.split('.')
        orbits.append(orbit)
        dates.append(date)
        times.append(time)
        core_ids.append(int(core_id))  # orig code used float
    f_in.close()
    ncases = len(orbits)

    if ncases > 0:
            
        # Read corresponding dat file (array(n_cases,numCols))
        dat = np.genfromtxt(file_dat)   # np.loadtxt(file) is faster
        if ncases == 1:
            cols = len(dat)
            dat = np.reshape(dat,(1,cols))
        #alternative using pandas
        #frame = pd.DataFrame(data=dat,
        #                     columns=['coreLon','coreLat','coreArea','coreTopHt',
        #                              'coreBotHt','coreDimX','coreDimY','coreTerHt',
        #                              'coreOceanLand','coreRainMean','coreRainStdDev',
        #                              'coreRainMax','coreRainMin',
        #                              'coreRainPixAll','coreRainPixStra','coreRainPixConv',
        #                              'coreRainMeanAll','coreRainMeanStra','coreRainMeanConv',
        #                              'coreRainVolAll','coreRainVolStra','coreRainVolConv',
        #                              'stormLon','stormLat','stormArea','stormTopHt',
        #                              'stormBotHt','stormDimX','stormDimY','stormTerHt',
        #                              'stormOceanLand','stormRainMean','stormRainStdDev',
        #                              'stormRainMax','stormRainMin',
        #                              'stormRainPixAll','stormRainPixStra','stormRainPixConv',
        #                              'stormRainMeanAll','stormRainMeanStra','stormRainMeanConv',
        #                              'stormRainVolAll','stormRainVolStra','stormRainVolConv',
                                          
        coreLon = dat[:,0]
        coreLat = dat[:,1]           
        coreArea = dat[:,2]
        coreTopHt = dat[:,3]
        coreBotHt = dat[:,4]
        coreDimX = dat[:,5]
        coreDimY = dat[:,6]
        coreTerHt = dat[:,7]
        coreOceanLand = dat[:,8]
        coreRainMean = dat[:,9]
        coreRainStdDev = dat[:,10]
        coreRainMax = dat[:,11]
        coreRainMin = dat[:,12]
        coreRainPixAll = dat[:,13]
        coreRainPixStra = dat[:,14]
        coreRainPixConv = dat[:,15]
        coreRainMeanAll = dat[:,16]
        coreRainMeanStra = dat[:,17]
        coreRainMeanConv = dat[:,18]
        coreRainVolAll = dat[:,19]
        coreRainVolStra = dat[:,20]
        coreRainVolConv = dat[:,21]
        coreDmMin = dat[:,22]
        coreDmMax = dat[:,23]
        coreDmMean = dat[:,24]
        coreDmStdDev = dat[:,25]
        
        stormLon = dat[:,26]
        stormLat = dat[:,27]
        stormArea = dat[:,28]
        stormTopHt = dat[:,29]
        stormBotHt = dat[:,30]
        stormDimX = dat[:,31]
        stormDimY = dat[:,32]
        stormTerHt = dat[:,33]
        stormOceanLand = dat[:,34]
        stormRainMean = dat[:,35]
        stormRainStdDev = dat[:,36]
        stormRainMax = dat[:,37]
        stormRainMin = dat[:,38]
        stormRainPixAll = dat[:,39]
        stormRainPixStra = dat[:,40]
        stormRainPixConv = dat[:,41]
        stormRainMeanAll = dat[:,42]
        stormRainMeanStra = dat[:,43]
        stormRainMeanConv = dat[:,44]
        stormRainVolAll = dat[:,45]
        stormRainVolStra = dat[:,46]
        stormRainVolConv = dat[:,47]
        stormDmMin = dat[:,48]
        stormDmMax = dat[:,49]
        stormDmMean = dat[:,50]
        stormDmStdDev = dat[:,51]

    #-------------------
    # Create netcdf file
    #-------------------
    outDir = outDirBase+'/'+feature+'/'+imonth
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    cdfFile = outDir+'/'+codeVersion+'_'+feature+'_'+thresh+'_'+data_type+'_'+iyear+imonth+'_'+region+'.nc'
    id = Dataset(cdfFile,'w',format='NETCDF4')

    # Make dimensions
    if ncases > 0:
        id.createDimension('case',ncases)
        id.createDimension('date_len',len(dates[0]))
        id.createDimension('time_len',len(times[0]))
        id.createDimension('orbit_len',len(orbits[0]))
    else:
        id.createDimension('case',ncases)
        id.createDimension('date_len',8)
        id.createDimension('time_len',6)
        id.createDimension('orbit_len',6)
    if verbose:
        print('Done with dims')

    # Define variables    
    orbit_id = id.createVariable('orbit','c',('case','orbit_len',),zlib=True,complevel=5) 
    orbit_id.units = 'none'
    orbit_id.long_name = 'orbit number'
        
    date_id = id.createVariable('date','c',('case','date_len',),zlib=True,complevel=5) 
    date_id.units = 'YYYYMMDD'
    date_id.long_name = 'date'
        
    utc_id = id.createVariable('time','c',('case','time_len',),zlib=True,complevel=5) 
    utc_id.units = 'hhmmss'
    utc_id.long_name = 'UTC time'
        
    core_id_id = id.createVariable('core_id','i4',('case',),zlib=True,complevel=5) 
    core_id_id.units = 'none'
    core_id_id.long_name = 'core id'
        
    coreLon_id = id.createVariable('core_lon','f4',('case',),zlib=True,complevel=5) 
    coreLon_id.units = 'degrees'
    coreLon_id.long_name = 'longitude of core centroid'
        
    coreLat_id = id.createVariable('core_lat','f4',('case',),zlib=True,complevel=5) 
    coreLat_id.units = 'degrees'
    coreLat_id.long_name = 'latitude of core centroid'
        
    coreArea_id = id.createVariable('core_area','f4',('case',),zlib=True,complevel=5) 
    coreArea_id.units = 'km^2'
    coreArea_id.long_name = 'core area'
        
    coreTopHt_id = id.createVariable('core_top_ht','f4',('case',),zlib=True,complevel=5) 
    coreTopHt_id.units = 'km'
    coreTopHt_id.long_name = 'top height of core'

    coreBotHt_id = id.createVariable('core_bot_ht','f4',('case',),zlib=True,complevel=5) 
    coreBotHt_id.units = 'km'
    coreBotHt_id.long_name = 'bottom height of core'

    coreDimX_id = id.createVariable('core_dim_x','f4',('case',),zlib=True,complevel=5) 
    coreDimX_id.units = 'degrees'
    coreDimX_id.long_name = 'max x-dimension of core'

    coreDimY_id = id.createVariable('core_dim_y','f4',('case',),zlib=True,complevel=5) 
    coreDimY_id.units = 'degrees'
    coreDimY_id.long_name = 'max y-dimension of core'

    coreTerHt_id = id.createVariable('core_ter_ht','f4',('case',),zlib=True,complevel=5) 
    coreTerHt_id.units = 'm'
    coreTerHt_id.long_name = 'terrain height at [core_lon,core_lat]'

    coreOL_id = id.createVariable('core_ocean_land_flag','i4',('case',),zlib=True,complevel=5) 
    coreOL_id.units = 'none'
    coreOL_id.long_name = 'terrain type at [core_lon,core_lat]'
    coreOL_id.ocean = '0'
    coreOL_id.land = '1'

    coreRainMean_id = id.createVariable('core_rain_rate_mean','f4',('case',),zlib=True,complevel=5) 
    coreRainMean_id.units = 'mm/hr'
    coreRainMean_id.long_name = 'mean rain rate of core'

    coreRainStdDev_id = id.createVariable('core_rain_rate_std_dev','f4',('case',),zlib=True,complevel=5) 
    coreRainStdDev_id.units = 'mm/hr'
    coreRainStdDev_id.long_name = 'standard deviation of core rain rate'

    coreRainMax_id = id.createVariable('core_rain_rate_max','f4',('case',),zlib=True,complevel=5) 
    coreRainMax_id.units = 'mm/hr'
    coreRainMax_id.long_name = 'max core rain rate'

    coreRainMin_id = id.createVariable('core_rain_rate_min','f4',('case',),zlib=True,complevel=5) 
    coreRainMin_id.units = 'mm/hr'
    coreRainMin_id.long_name = 'min core rain rate'

    coreRainPixAll_id = id.createVariable('core_rain_pix_all','i4',('case',),zlib=True,complevel=5) 
    coreRainPixAll_id.units = 'none'
    coreRainPixAll_id.long_name = 'total number of pixels in core'

    coreRainPixStra_id = id.createVariable('core_rain_pix_stra','i4',('case',),zlib=True,complevel=5) 
    coreRainPixStra_id.units = 'none'
    coreRainPixStra_id.long_name = 'number of stratiform pixels in core'

    coreRainPixConv_id = id.createVariable('core_rain_pix_conv','i4',('case',),zlib=True,complevel=5) 
    coreRainPixConv_id.units = 'none'
    coreRainPixConv_id.long_name = 'number of convective pixels in core'

    coreRainMeanAll_id = id.createVariable('core_rain_rate_mean_all','f4',('case',),zlib=True,complevel=5) 
    coreRainMeanAll_id.units = 'mm/hr'
    coreRainMeanAll_id.long_name = 'mean rain rate of all core precip'

    coreRainMeanStra_id = id.createVariable('core_rain_rate_mean_stra','f4',('case',),zlib=True,complevel=5) 
    coreRainMeanStra_id.units = 'mm/hr'
    coreRainMeanStra_id.long_name = 'mean rain rate of stratiform core precip'

    coreRainMeanConv_id = id.createVariable('core_rain_rate_mean_conv','f4',('case',),zlib=True,complevel=5) 
    coreRainMeanConv_id.units = 'mm/hr'
    coreRainMeanConv_id.long_name = 'mean rain rate of convective core precip'

    coreRainVolAll_id = id.createVariable('core_rain_vol_all','f4',('case',),zlib=True,complevel=5) 
    coreRainVolAll_id.units = '10^6 kg/s'
    coreRainVolAll_id.long_name = 'volume of all core precip'

    coreRainVolStra_id = id.createVariable('core_rain_vol_stra','f4',('case',),zlib=True,complevel=5) 
    coreRainVolStra_id.units = '10^6 kg/s'
    coreRainVolStra_id.long_name = 'volume of stratiform core precip'

    coreRainVolConv_id = id.createVariable('core_rain_vol_conv','f4',('case',),zlib=True,complevel=5) 
    coreRainVolConv_id.units = '10^6 kg/s'
    coreRainVolConv_id.long_name = 'volume of convective core precip'

    coreDmMin_id = id.createVariable('core_dm_min','f4',('case',),zlib=True,complevel=5) 
    coreDmMin_id.units = 'mm'
    coreDmMin_id.long_name = 'minimum Dm'

    coreDmMax_id = id.createVariable('core_dm_max','f4',('case',),zlib=True,complevel=5) 
    coreDmMax_id.units = 'mm'
    coreDmMax_id.long_name = 'maximum Dm'

    coreDmMean_id = id.createVariable('core_dm_mean','f4',('case',),zlib=True,complevel=5) 
    coreDmMean_id.units = 'mm'
    coreDmMean_id.long_name = 'mean Dm'

    coreDmStdDev_id = id.createVariable('core_dm_std_dev','f4',('case',),zlib=True,complevel=5) 
    coreDmStdDev_id.units = 'mm'
    coreDmStdDev_id.long_name = 'Dm standard deviation'

    stormLon_id = id.createVariable('storm_lon','f4',('case',),zlib=True,complevel=5) 
    stormLon_id.units = 'degrees'
    stormLon_id.long_name = 'longitude of storm centroid'
        
    stormLat_id = id.createVariable('storm_lat','f4',('case',),zlib=True,complevel=5) 
    stormLat_id.units = 'degrees'
    stormLat_id.long_name = 'latitude of storm centroid'
        
    stormArea_id = id.createVariable('storm_area','f4',('case',),zlib=True,complevel=5) 
    stormArea_id.units = 'km^2'
    stormArea_id.long_name = 'storm area'
        
    stormTopHt_id = id.createVariable('storm_top_ht','f4',('case',),zlib=True,complevel=5) 
    stormTopHt_id.units = 'km'
    stormTopHt_id.long_name = 'top height of storm'

    stormBotHt_id = id.createVariable('storm_bot_ht','f4',('case',),zlib=True,complevel=5) 
    stormBotHt_id.units = 'km'
    stormBotHt_id.long_name = 'bottom height of storm'

    stormDimX_id = id.createVariable('storm_dim_x','f4',('case',),zlib=True,complevel=5) 
    stormDimX_id.units = 'degrees'
    stormDimX_id.long_name = 'max x-dimension of storm'

    stormDimY_id = id.createVariable('storm_dim_y','f4',('case',),zlib=True,complevel=5) 
    stormDimY_id.units = 'degrees'
    stormDimY_id.long_name = 'max y-dimension of storm'

    stormTerHt_id = id.createVariable('storm_ter_ht','f4',('case',),zlib=True,complevel=5) 
    stormTerHt_id.units = 'm'
    stormTerHt_id.long_name = 'terrain height at [storm_lon,storm_lat]'

    stormOL_id = id.createVariable('storm_ocean_land_flag','i4',('case',),zlib=True,complevel=5) 
    stormOL_id.units = 'none'
    stormOL_id.long_name = 'terrain type at [storm_lon,storm_lat]'
    stormOL_id.ocean = '0'
    stormOL_id.land = '1'

    stormRainMean_id = id.createVariable('storm_rain_mean','f4',('case',),zlib=True,complevel=5) 
    stormRainMean_id.units = 'mm/hr'
    stormRainMean_id.long_name = 'mean rain rate of storm'

    stormRainStdDev_id = id.createVariable('storm_rain_rate_std_dev','f4',('case',),zlib=True,complevel=5) 
    stormRainStdDev_id.units = 'mm/hr'
    stormRainStdDev_id.long_name = 'standard deviation of storm rain rate'

    stormRainMax_id = id.createVariable('storm_rain_rate_max','f4',('case',),zlib=True,complevel=5) 
    stormRainMax_id.units = 'mm/hr'
    stormRainMax_id.long_name = 'max storm rain rate'
            
    stormRainMin_id = id.createVariable('storm_rain_rate_min','f4',('case',),zlib=True,complevel=5) 
    stormRainMin_id.units = 'mm/hr'
    stormRainMin_id.long_name = 'min storm rain rate'

    stormRainPixAll_id = id.createVariable('storm_rain_pix_all','i4',('case',),zlib=True,complevel=5) 
    stormRainPixAll_id.units = 'none'
    stormRainPixAll_id.long_name = 'total number of pixels in storm (stra, conv, other, anvil, etc)'

    stormRainPixStra_id = id.createVariable('storm_rain_pix_stra','i4',('case',),zlib=True,complevel=5) 
    stormRainPixStra_id.units = 'none'
    stormRainPixStra_id.long_name = 'total number of stratiform pixels in storm'

    stormRainPixConv_id = id.createVariable('storm_rain_pix_conv','i4',('case',),zlib=True,complevel=5) 
    stormRainPixConv_id.units = 'none'
    stormRainPixConv_id.long_name = 'total number of convective pixels in storm'

    stormRainMeanAll_id = id.createVariable('storm_rain_rate_mean_all','f4',('case',),zlib=True,complevel=5) 
    stormRainMeanAll_id.units = 'mm/hr'
    stormRainMeanAll_id.long_name = 'mean rain rate of all storm precip'

    stormRainMeanStra_id = id.createVariable('storm_rain_rate_mean_stra','f4',('case',),zlib=True,complevel=5) 
    stormRainMeanStra_id.units = 'mm/hr'
    stormRainMeanStra_id.long_name = 'mean rain rate of stratiform storm precip'

    stormRainMeanConv_id = id.createVariable('storm_rain_rate_mean_conv','f4',('case',),zlib=True,complevel=5) 
    stormRainMeanConv_id.units = 'mm/hr'
    stormRainMeanConv_id.long_name = 'mean rain rate of convective storm precip'

    stormRainVolAll_id = id.createVariable('storm_rain_vol_all','f4',('case',),zlib=True,complevel=5) 
    stormRainVolAll_id.units = '10^6 kg/s'
    stormRainVolAll_id.long_name = 'volume of all storm precip'

    stormRainVolStra_id = id.createVariable('storm_rain_vol_stra','f4',('case',),zlib=True,complevel=5) 
    stormRainVolStra_id.units = '10^6 kg/s'
    stormRainVolStra_id.long_name = 'volume of stratiform storm precip'

    stormRainVolConv_id = id.createVariable('storm_rain_vol_conv','f4',('case',),zlib=True,complevel=5) 
    stormRainVolConv_id.units = '10^6 kg/s'
    stormRainVolConv_id.long_name = 'volume of convective storm precip'

    stormDmMin_id = id.createVariable('storm_dm_min','f4',('case',),zlib=True,complevel=5) 
    stormDmMin_id.units = 'mm'
    stormDmMin_id.long_name = 'minimum Dm'

    stormDmMax_id = id.createVariable('storm_dm_max','f4',('case',),zlib=True,complevel=5) 
    stormDmMax_id.units = 'mm'
    stormDmMax_id.long_name = 'maximum Dm'

    stormDmMean_id = id.createVariable('storm_dm_mean','f4',('case',),zlib=True,complevel=5) 
    stormDmMean_id.units = 'mm'
    stormDmMean_id.long_name = 'mean Dm'

    stormDmStdDev_id = id.createVariable('storm_dm_std_dev','f4',('case',),zlib=True,complevel=5) 
    stormDmStdDev_id.units = 'mm'
    stormDmStdDev_id.long_name = 'Dm standard deviation'

    if thresh == 'xxx':
        id.Title = feature_long+' Cases Based on TRMM PR /NS/SLV/precipRateNearSurface'
        id.ShallowIsolated_Criterion = 'TRMM PR /NS/CSF/flagShallowRain = [10,11]'
    elif thresh == 'str':
        id.Title = feature_long+' Cases Based on TRMM PR /NS/SLV/precipRateNearSurface Using Strong (Land) Thresholds'
        if feature == 'BSR':
            id.BroadStratiform_Criterion = 'contiguous stratiform echo >= 50,000km^2'
        elif feature == 'DCC':
            id.DeepConvective_Criterion = '40dBZ echos exceed 10km in height'
        elif feature == 'WCC':
            id.WideConvective_Criterion = '40dBZ echos with horiz dim at some alt > 1000km^2'
        elif feature == 'DWC':
            id.DeepWideConvective_Criterion = 'meet both DeepConvective and WideConvective criteria'
        else:
            id.Title = feature_long+' Cases Based on TRMM PR /NS/SLV/precipRateNearSurface Using Moderate (Ocean) Thresholds'
            if feature == 'BSR':
                id.BroadStratiform_Criterion = 'contiguous stratiform echo >= 40,000km^2'
            elif feature == 'DCC':
                id.DeepConvective_Criterion = '30dBZ echos exceed 8km in height'
            elif feature == 'WCC':
                id.WideConvective_Criterion = '30dBZ echos with horiz dim at some alt > 800km^2'
            elif feature == 'DWC':
                id.DeepWideConvective_Criterion = 'meet both DeepConvective and WideConvective criteria'

        id.source_uw = 'Created by the Mesoscale Group, University of Washington'
        id.source_nasa = 'Original data obtained from NASA Goddard Earth Sciences http://pmm.nasa.gov/'
        id.reference = 'Methodology described in Houze et al. 2015 (Reviews of Geophysics)'
        id.data_location = 'http://trmm.atmos.washington.edu'
        id.region = region+' '+region_long
        id.lon_min = limits[1]
        id.lon_max = limits[3]
        id.lat_min = limits[0]
        id.lat_max = limits[2]

        # Write the data into file
        if ncases > 0:
            orbit_id[:] = orbits
            date_id[:] = dates
            utc_id[:] = times
            core_id_id[:] = core_ids
            coreLon_id[:] = coreLon
            coreLat_id[:] = coreLat
            coreArea_id[:] = coreArea
            coreTopHt_id[:] = coreTopHt
            coreBotHt_id[:] = coreBotHt
            coreDimX_id[:] = coreDimX
            coreDimY_id[:] = coreDimY
            coreTerHt_id[:] = coreTerHt
            coreOL_id[:] = coreOceanLand
            coreRainMean_id[:] = coreRainMean
            coreRainStdDev_id[:] = coreRainStdDev
            coreRainMax_id[:] = coreRainMax
            coreRainMin_id[:] = coreRainMin
            coreRainPixAll_id[:] = coreRainPixAll
            coreRainPixStra_id[:] = coreRainPixStra
            coreRainPixConv_id[:] = coreRainPixConv
            coreRainMeanAll_id[:] = coreRainMeanAll
            coreRainMeanStra_id[:] = coreRainMeanStra
            coreRainMeanConv_id[:] = coreRainMeanConv
            coreRainVolAll_id[:] = coreRainVolAll
            coreRainVolStra_id[:] = coreRainVolStra
            coreRainVolConv_id[:] = coreRainVolConv
            coreDmMin_id[:] = coreDmMin
            coreDmMax_id[:] = coreDmMax
            coreDmMean_id[:] = coreDmMean
            coreDmStdDev_id[:] = coreDmStdDev
                
            stormLon_id[:] = stormLon
            stormLat_id[:] = stormLat
            stormArea_id[:] = stormArea
            stormTopHt_id[:] = stormTopHt
            stormBotHt_id[:] = stormBotHt
            stormDimX_id[:] = stormDimX
            stormDimY_id[:] = stormDimY
            stormTerHt_id[:] = stormTerHt
            stormOL_id[:] = stormOceanLand
            stormRainMean_id[:] = stormRainMean
            stormRainStdDev_id[:] = stormRainStdDev
            stormRainMax_id[:] = stormRainMax
            stormRainMin_id[:] = stormRainMin
            stormRainPixAll_id[:] = stormRainPixAll
            stormRainPixStra_id[:] = stormRainPixStra
            stormRainPixConv_id[:] = stormRainPixConv
            stormRainMeanAll_id[:] = stormRainMeanAll
            stormRainMeanStra_id[:] = stormRainMeanStra
            stormRainMeanConv_id[:] = stormRainMeanConv
            stormRainVolAll_id[:] = stormRainVolAll
            stormRainVolStra_id[:] = stormRainVolStra
            stormRainVolConv_id[:] = stormRainVolConv
            stormDmMin_id[:] = stormDmMin
            stormDmMax_id[:] = stormDmMax
            stormDmMean_id[:] = stormDmMean
            stormDmStdDev_id[:] = stormDmStdDev

            # close netcdf file
            id.close()         


