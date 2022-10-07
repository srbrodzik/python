import os
import netCDF4 as nc4
import numpy as np

region = 'H02'
indir = '/home/disk/archive3/gpm'
eventType = 'WCC'
eventThres = 'str'
years = ['2014','2015']
#years = ['2014','2015','2016','2017']
#months = ['03','04','05','06','07','08','09','10','11','12']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
LAT_THRES = 60.0
outdir = '/home/disk/archive3/gpm'
outfile = 'highLat_'+eventType+'_'+region+'_'+eventThres+'.txt'

num_years = len(years)
num_months = len(months)

# open output file
fout = open(outdir+'/'+outfile, 'w')
fout.write("orbit\tdate\t\ttime\tid\tlon\tlat\tarea\ttop_ht\tbot_ht\n")

for iyear in range(0,num_years):
    for imonth in range(0,num_months):
        print 'year = ',years[iyear],' month = ',months[imonth]
        year_month = years[iyear]+months[imonth]
        substring = eventThres+'_tab_'+year_month
        for fname in os.listdir(indir+'/'+region+'/'+eventType+'/'+months[imonth]):
            if substring in fname:
                ncid = nc4.Dataset(indir+'/'+region+'/'+eventType+'/'+months[imonth]+'/'+fname,'r')
                core_lat = np.array(ncid.variables['core_lat'])
                if np.sum(core_lat >= LAT_THRES) > 0:
                    core_lat_sub = core_lat[core_lat >= LAT_THRES]
                
                    orbit = nc4.chartostring(ncid.variables['orbit'][:,:])
                    orbit_sub = orbit[core_lat >= LAT_THRES]
                
                    date = nc4.chartostring(ncid.variables['date'][:,:])
                    date_sub = date[core_lat >= LAT_THRES]
                
                    time = nc4.chartostring(ncid.variables['time'][:,:])
                    time_sub = time[core_lat >= LAT_THRES]
                
                    core_id = np.array(ncid.variables['core_id'])
                    core_id_sub = core_id[core_lat >= LAT_THRES]
                
                    core_lon = np.array(ncid.variables['core_lon'])
                    core_lon_sub = core_lon[core_lat >= LAT_THRES]
                
                    core_area = np.array(ncid.variables['core_area'])
                    core_area_sub = core_area[core_lat >= LAT_THRES]
                
                    core_top_ht = np.array(ncid.variables['core_top_ht'])
                    core_top_ht_sub = core_top_ht[core_lat >= LAT_THRES]
                
                    core_bot_ht = np.array(ncid.variables['core_bot_ht'])
                    core_bot_ht_sub = core_bot_ht[core_lat >= LAT_THRES]

                    for i in range(0,len(core_lat_sub)):
                        fout.write("%6s\t%8s\t%6s\t%2d\t%7.2f\t%5.2f\t%7.2f\t%5.2f\t%5.2f\n" %
                                   (orbit_sub[i],date_sub[i],time_sub[i],core_id_sub[i],core_lon_sub[i],
                                    core_lat_sub[i],core_area_sub[i],core_top_ht_sub[i],core_bot_ht_sub[i]))

                ncid.close()

fout.close()
