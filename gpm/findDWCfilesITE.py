import os
import glob
import shutil

region = 'NHW'

#flist = '/home/disk/bob/gpm/nht_ku/classify/class_data_ite132/monthly_class_v10s/DeepWide_Convective_'+region+'_v10s.info'
flist = '/home/disk/bob/gpm/nht_ku/classify/class_data_v05_uw/monthly_class_v10s/Wide_Convective_'+region+'_v10s.info'
indir = '/home/disk/bob/gpm/nht_ku/classify/ex_data_ite132'
outdir = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_ite/netcdf4_wcc_str'

#flist = '/home/disk/bob/gpm/nht_ku/classify/class_data_v05_uw/monthly_class_v10s/DeepWide_Convective_'+region+'_v10s.info'
#indir = '/home/disk/bob/gpm/nht_ku/classify/ex_data_v05'
#outdir = '/home/disk/bob/gpm/nht_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dwc_str'

with open(flist) as f:
    for line in f:
        line = line.strip('\n')
        print line
 
        parts = line.split('.')
        orbit = parts[0]
        date = parts[1]
        time = parts[2]

        year = date[0:4]
        month = date[4:6]

        fname1 = 'GPM2Ku5_uw3_'+date+'.'+time+'_to_*_'+orbit+'_'+region+'.nc'
        fname2 = glob.glob(indir+'/'+year+'/'+month+'/'+fname1)
        fname = fname2[0]
        fname_new = outdir+'/GPM2Ku5_uw3_'+date+'.'+time+'_'+orbit+'_'+region+'.nc'

        shutil.copy2(fname,fname_new)
