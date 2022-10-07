import os
import glob

season = 'MAM'
subregion = 'ChinaNC'
region = 'ASIA'

flist_in = '/home/disk/bob/gpm/asia_ku/classify/class_data_v05_uw/monthly_class_v10s/Wide_Convective_'+subregion+'_1418_'+season+'_v10s.info'
indir = '/home/disk/bob/gpm/asia_ku/classify/ex_data_v05'
flist_out_net4 = '/home/disk/bob/gpm/asia_ku/classify/class_data_v05_uw/monthly_class_v10s/Wide_Convective_'+subregion+'_1418_'+season+'_v10s.fnames_net4'
flist_out_zeb  = '/home/disk/bob/gpm/asia_ku/classify/class_data_v05_uw/monthly_class_v10s/Wide_Convective_'+subregion+'_1418_'+season+'_v10s.fnames_zeb'

f_net4 = open(flist_out_net4,'w')
f_zeb = open(flist_out_zeb,'w')

with open(flist_in) as f:
    for line in f:
        line = line.strip('\n')
        print line
 
        parts = line.split('\t')
        orbit = parts[0]
        date = parts[1]
        time = parts[2]
        #print orbit, date, time

        year = date[0:4]
        month = date[4:6]
        #print year, month

        fname1 = 'GPM2Ku5_uw3_'+date+'.'+time+'_to_*_'+orbit+'_'+region+'.nc'
        #print fname1
        fname2 = glob.glob(indir+'/'+year+'/'+month+'/'+fname1)
        #print fname2
        fname_full = fname2[0]
        fname_net4 = os.path.basename(fname_full)
        f_net4.write(fname_net4+"\n")
        fname_zeb = 'GPM2Ku5_uw3_'+date+'.'+time+'_'+orbit+'_'+region+'.nc'
        f_zeb.write(fname_zeb+"\n")

f_net4.close()
f_zeb.close()


