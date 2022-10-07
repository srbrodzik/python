import os
import glob
import shutil

#flist = '/home/disk/bob/gpm/nam_ku/classify/class_data_ite132/monthly_class_v9/BroadStratiform_SON_all_v9.info'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_ite132'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_ite/netcdf4_bsr_str'

#flist = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw/monthly_class_v9/BroadStratiform_SON_all_v9.info'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_bsr_str'

#flist = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw_testing/monthly_class_v10/BroadStratiform_JJA_1416J_v10.info'
#flist = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw/monthly_class_v10s/BroadStratiform_JJA_1417_NAM_v10s.info'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_bsr_str_mask'

#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v06'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v06/netcdf4_bsr_str'
#indir = '/home/disk/archive3/gpm/NAM/interp_data'
#outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str_mask'
#flist = outdir+'/nam_wcc_str_files.csv'

indir = '/home/disk/archive3/gpm/NAM/interp_data/v05'
outdir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str_mask'
flist = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw/monthly_class_v10s/Wide_Convective_Gulf_JJA_2018_v10s.csv'

with open(flist) as f:
    for line in f:
        line = line.strip('\n')
        print line
 
#        parts = line.split('.')
        parts = line.split(',')
        orbit = parts[0]
        date = parts[1]
        time = parts[2]

        year = date[0:4]
        month = date[4:6]

        fname1 = 'GPM2Ku5_uw3_'+date+'.'+time+'_to_*_'+orbit+'_NAM.nc'
        fname2 = glob.glob(indir+'/'+year+'/'+month+'/'+fname1)
        fname = fname2[0]
        fname_new = outdir+'/GPM2Ku5_uw3_'+date+'.'+time+'_'+orbit+'_NAM.nc'
        print fname
        print fname_new
        
        shutil.copy2(fname,fname_new)
