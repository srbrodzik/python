# This function is used to remove files from ex_data_v05 for H03 that do not fall in proper domain

# Imports and functions
import os
import netCDF4 as nc4
import numpy as np
import shutil

# User inputs
indir = '/home/disk/bob/gpm/h03_ku/classify/ex_data_v05'
#months = ['01','02','03','04','05','06','07','08','09','10','11','12']
#years = ['2014','2015','2016','2017']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
years = ['2014','2015','2016','2017','2018']

# H03 domain
minLatitude = -67
maxLatitude = -35
minLongitude = -30
maxLongitude = 75

for iyear in years:
    for imonth in months:
        # Create outOfBounds_partial directory if it doesn't exist
        if not os.path.exists(indir+'/'+iyear+'/'+imonth+'/outOfBounds_partial'):
            os.makedirs(indir+'/'+iyear+'/'+imonth+'/outOfBounds_partial')
        os.chdir(indir+'/'+iyear+'/'+imonth)
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
            if fname.endswith('nc'):

                print fname

                # open and read latfrom input file
                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'r')
                lat = np.array(ncid.variables['lat'])
                ncid.close()

                # find max lat
                maxLat = np.amax(lat)
                print maxLat
                
                if maxLat > maxLatitude:
                    fname_parts = fname.split('.')
                    outfile = fname_parts[0]+'.'+fname_parts[1]+'.'+fname_parts[2]+'.new.nc'
                    command = 'ncks -a -v lon,lat,rain_type,rain_type_raw,phase_type,phase_type_raw,shallow_rain_type,near_surf_rain,width_bb,height_bb,swath,flag_echo,refl,rain_type_uw,bsr_mask_str,dcc_mask_str,dwc_mask_str,wcc_mask_str,storm_mask_str,bsr_mask_mod,dcc_mask_mod,dwc_mask_mod,wcc_mask_mod,storm_mask_mod -d lat,-67.,-35. -d lon,-30.,75. '+fname+' -O '+outfile
                    #print command
                    os.system(command)
                    command = 'mv '+fname+' outOfBounds_partial'
                    os.system(command)
                    
