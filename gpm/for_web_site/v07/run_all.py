#!/usr/bin/python3

import os
import sys

if len(sys.argv) != 2:
    print('Usage: {} region[AFC|AKA|CIO|EPO|EUR|H01|H02|H03|H04|H05|NAM|NAS|SAM|SAS|WMP]'.format(sys.argv[0]))
    sys.exit()
else:
    region = sys.argv[1]

binDir = '/home/disk/meso-home/brodzik/python/gpm/for_web_site/v07'
    
years = ['2023']
months = ['05']
#years = ['2017']
#months = ['03','04','05','06','07','08','09','10','11','12']
#months = ['01','02','03','04','05','06','07','08','09','10','11','12']

for year in years:
    for month in months:
        cmd = binDir+'/convert_info_rain_nsr_bsr_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' str'
        os.system(cmd)
        cmd = binDir+'/convert_info_rain_nsr_bsr_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' mod'
        os.system(cmd)
        cmd = binDir+'/convert_info_rain_nsr_conv_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' str'
        os.system(cmd)
        cmd = binDir+'/convert_info_rain_nsr_conv_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' mod'
        os.system(cmd)
        cmd = binDir+'/convert_info_rain_nsr_shi_grd_cdf_to_cdf.py '+region+' '+year+' '+ month
        os.system(cmd)
        cmd = binDir+'/convert_monthly_class_tab_ascii_to_cdf.py '+region+' '+year+' '+ month+' str'
        os.system(cmd)
        cmd = binDir+'/convert_monthly_class_tab_ascii_to_cdf.py '+region+' '+year+' '+ month+' mod'
        os.system(cmd)
        cmd = binDir+'/convert_cfad_dbz_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' str'
        os.system(cmd)
        cmd = binDir+'/convert_cfad_dbz_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' mod'
        os.system(cmd)
        cmd = binDir+'/convert_cfad_dm_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' str'
        os.system(cmd)
        cmd = binDir+'/convert_cfad_dm_grd_cdf_to_cdf.py '+region+' '+year+' '+ month+' mod'
        os.system(cmd)
