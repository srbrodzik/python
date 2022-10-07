import os
import netCDF4 as nc4
import numpy as np
import fnmatch as fnm

region = 'ASIA'
indir = '/home/disk/archive3/gpm/'+region
seasons = ['DJF','MAM','JJA','SON']
#region = 'NAM'
#indir = '/home/disk/archive3/gpm/'+region
#seasons = ['DJF','MAM','JJA','SON']
#events = ['DWC']
events = ['DWC']
thresholds = ['str','mod']
outdir = '/home/disk/shear2/brodzik/IDL/ERAi_analysis/data'

for season in seasons:
    print season
    if season == 'DJF':
        months = ['12','01','02']
    elif season == 'MAM':
        months = ['03','04','05']
    elif season == 'JJA':
        months = ['06','07','08']
    elif season == 'SON':
        months = ['09','10','11']
    else:
        print 'season='+season+' is invalid . . . exiting'
        exit()
    for event in events:
        print event
        for threshold in thresholds:
            print threshold
            orbit_all=np.zeros(1)
            date_all=np.zeros(1)
            time_all=np.zeros(1)
            core_lon_all=np.zeros(1)
            core_lat_all=np.zeros(1)
            core_area_all=np.zeros(1)
            core_top_ht_all=np.zeros(1)
            core_bot_ht_all=np.zeros(1)
            core_dim_x_all=np.zeros(1)
            core_dim_y_all=np.zeros(1)
            core_ter_ht_all=np.zeros(1)
            core_ocean_land_flag_all=np.zeros(1)
            for month in months:
                print month
                for fname in os.listdir(indir+'/'+event+'/'+month):
                    if fnm.fnmatch(fname,'*'+threshold+'_tab*'):
                        print fname
                        ncid = nc4.Dataset(indir+'/'+event+'/'+month+'/'+fname,'r')
                        orbit = nc4.chartostring(ncid.variables['orbit'][:,:])
                        date = nc4.chartostring(ncid.variables['date'][:,:])
                        time = nc4.chartostring(ncid.variables['time'][:,:])
                        core_lon = np.array(ncid.variables['core_lon'])
                        core_lat = np.array(ncid.variables['core_lat'])
                        core_area = np.array(ncid.variables['core_area'])
                        core_top_ht = np.array(ncid.variables['core_top_ht'])
                        core_bot_ht = np.array(ncid.variables['core_bot_ht'])
                        core_dim_x = np.array(ncid.variables['core_dim_x'])
                        core_dim_y = np.array(ncid.variables['core_dim_y'])
                        core_ter_ht = np.array(ncid.variables['core_ter_ht'])
                        core_ocean_land_flag = np.array(ncid.variables['core_ocean_land_flag'])
                        ncid.close()
                        orbit_all = np.concatenate([orbit_all,orbit])
                        date_all = np.concatenate([date_all,date])
                        time_all = np.concatenate([time_all,time])
                        core_lon_all = np.concatenate([core_lon_all,core_lon])
                        core_lat_all = np.concatenate([core_lat_all,core_lat])
                        core_area_all = np.concatenate([core_area_all,core_area])
                        core_top_ht_all = np.concatenate([core_top_ht_all,core_top_ht])
                        core_bot_ht_all = np.concatenate([core_bot_ht_all,core_bot_ht])
                        core_dim_x_all = np.concatenate([core_dim_x_all,core_dim_x])
                        core_dim_y_all = np.concatenate([core_dim_y_all,core_dim_y])
                        core_ter_ht_all = np.concatenate([core_ter_ht_all,core_ter_ht])
                        core_ocean_land_flag_all = np.concatenate([core_ocean_land_flag_all,core_ocean_land_flag])
            print 'Done with months for '+threshold
            
            # remove first element of arrays
            orbit_all = orbit_all[1:]
            date_all = date_all[1:]
            time_all =  time_all[1:]
            core_lon_all = core_lon_all[1:]
            core_lat_all = core_lat_all[1:]
            core_area_all = core_area_all[1:]
            core_top_ht_all = core_top_ht_all[1:]
            core_bot_ht_all = core_bot_ht_all[1:]
            core_dim_x_all = core_dim_x_all[1:]
            core_dim_y_all = core_dim_y_all[1:]
            core_ter_ht_all = core_ter_ht_all[1:]
            core_ocean_land_flag_all = core_ocean_land_flag_all[1:]

            # write file for threshold, event and season
            n_cases = len(core_lon_all)
            fid = open(outdir+'/'+region+'_'+event+'_'+threshold+'_'+season+'_uw.txt','w')
            fid.write('%5d\n' % (n_cases) )
            for i in range(0,n_cases):
                fid.write('%6s %8s %6s %7.2f %6.2f %8.2f %5.2f %5.2f %5.2f %5.2f %5d %1d\n' % (orbit_all[i],date_all[i],time_all[i],core_lon_all[i],core_lat_all[i],core_area_all[i],core_top_ht_all[i],core_bot_ht_all[i],core_dim_x_all[i],core_dim_y_all[i],core_ter_ht_all[i],core_ocean_land_flag_all[i]))
            fid.close()
    
    
