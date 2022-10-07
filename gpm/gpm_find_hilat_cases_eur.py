# Imports and functions
import os
import netCDF4 as nc4
import numpy as np
from shutil import copyfile
import sys

# User inputs
indir = '/home/disk/bob/gpm/eur_ku/classify/ex_data_v05'
outdir = '/home/disk/bob/gpm/eur_ku/classify/ex_data_v05_hilat'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
years = ['2014','2015','2016','2017','2018']

minLatitude = 55.0

for iyear in years:
    for imonth in months:
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
            if fname.endswith('nc'):

                print fname

                # open and read lat from input file
                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'r')
                lat = np.array(ncid.variables['lat'])
                ncid.close()

                # find max lat
                maxLat = np.amax(lat)
                
                if maxLat > minLatitude:
                    infile = indir+'/'+iyear+'/'+imonth+'/'+fname
                    outfile = outdir+'/'+iyear+'/'+imonth+'/'+fname
                    outfileNew = outdir+'/'+iyear+'/'+imonth+'/'+fname+'.new'
                    print 'Copying '+infile+' to '+outfile
                    copyfile(infile,outfile)


                    command = "ncatted -a 'BroadStratiform_Criteria_Strong',global,d,, -h "+outfile
                    os.system(command)
                    command = "ncatted -a 'DeepConvective_Criteria_Strong',global,d,, -h "+outfile
                    os.system(command)
                    command = "ncatted -a 'WideConvective_Criteria_Strong',global,d,, -h "+outfile
                    os.system(command)
                    command = "ncatted -a 'DeepWideConvective_Criteria_Strong',global,d,, -h "+outfile
                    os.system(command)
                    command = "ncatted -a 'BroadStratiform_Criteria_Moderate',global,d,, -h "+outfile
                    os.system(command)
                    command = "ncatted -a 'DeepConvective_Criteria_Moderate',global,d,, -h "+outfile
                    os.system(command)
                    command = "ncatted -a 'WideConvective_Criteria_Moderate',global,d,, -h "+outfile
                    os.system(command)
                    command = "ncatted -a 'DeepWideConvective_Criteria_Moderate',global,d,, -h "+outfile
                    os.system(command)
                    command = 'ncks -x -v bsr_mask_str,dcc_mask_str,dwc_mask_str,wcc_mask_str,storm_mask_str,bsr_mask_mod,dcc_mask_mod,dwc_mask_mod,wcc_mask_mod,storm_mask_mod '+outfile+' '+outfileNew
                    os.system(command)
                    command = '/bin/mv '+outfileNew+' '+outfile
                    os.system(command)
                    command = 'ncks -d lon,-180.0,180.0 -d lat,53.0,67.0 '+outfile+' '+outfileNew
                    os.system(command)
                    command = '/bin/mv '+outfileNew+' '+outfile
                    os.system(command)
