#!/usr/bin/python3

#Wrapper script for running allStorms

#from deg2km import deg2km
#from findStormNew import findStorm
#from allStorms import allStorms_v11_v06

# define inputs
years = ['2014','2015'] 
months = ['01','02']
region = 'WMP'
path_in_base = '/home/disk/bob/gpm'
path_out_base = '/home/disk/bob/gpm'
path_topo = '/home/disk/shear2/brodzik/GTOPO30'
data_version = 'v06'

if region == 'WMP':
    limits = [-40.0,105.0,10.0,178.0]
    path_in = path_in_base+'/'+region.lower()+'_ku/classify/ex_data_'+data_version
    path_out = path_out_base+'/'+region.lower()+'_ku/classify/class_data_'+data_version
    topo_file = path_topo+'/gtopo5km_WMP_gpm_v06.hdf'
elif region == 'SAS':
    limits = [5.0,55.0,40.0,130.0]
    path_in = path_in_base+'/'+region.lower()+'_ku/classify/ex_data_'+data_version
    path_out = path_out_base+'/'+region.lower()+'_ku/classify/class_data_'+data_version
    topo_file = path_topo+'/gtopo5km_SAS_gpm.hdf'

for year in years:

    #if year == '2014':
    #    months = ['12']
    #elif year == '2015':
    #    months = ['01','02']
    #else:
    #    months = ['01','02','03','04','05','06','07','08','09','10','11','12']  

     for month in months:
         typeFile = month+'_'+year+'_'+region
         print('typeFile = ',typeFile)
         #allStorms_v11_v06,region,limits,path_in,path_out,topo_file,typeFile



