import os
import glob

#flist_in = '/home/disk/bob/gpm/nam_ku/classify/class_data_ite132/monthly_class_v9/DeepWide_Convective_MAM_all_v9.info'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_ite132'
#flist_out_net4 = '/home/disk/bob/gpm/nam_ku/classify/class_data_ite132/monthly_class_v9/DeepWide_Convective_MAM_v9.fnames_net4'
#flist_out_zeb  = '/home/disk/bob/gpm/nam_ku/classify/class_data_ite132/monthly_class_v9/DeepWide_Convective_MAM_v9.fnames_zeb'

#flist_in = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05/monthly_class_v9/DeepWide_Convective_SON_all_v9.info'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
#flist_out_net4 = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05/monthly_class_v9/DeepWide_Convective_SON_v9.fnames_net4'
#flist_out_zeb  = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05/monthly_class_v9/DeepWide_Convective_SON_v9.fnames_zeb'

#flist_in = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw/monthly_class_v9/DeepWide_Convective_MAM_all_v9.info'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
#flist_out_net4 = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw/monthly_class_v9/DeepWide_Convective_MAM_v9.fnames_net4'
#flist_out_zeb  = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw/monthly_class_v9/DeepWide_Convective_MAM_v9.fnames_zeb'

flist_in = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw_testing/monthly_class_v10/DeepWide_Convective_JJA_1415_v10.info'
indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
flist_out_net4 = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw_testing/monthly_class_v10/DeepWide_Convective_JJA_v10.fnames_net4'
flist_out_zeb  = '/home/disk/bob/gpm/nam_ku/classify/class_data_v05_uw_testing/monthly_class_v10/DeepWide_Convective_JJA_v10.fnames_zeb'

f_net4 = open(flist_out_net4,'w')
f_zeb = open(flist_out_zeb,'w')

with open(flist_in) as f:
    for line in f:
        line = line.strip('\n')
        print line
 
        parts = line.split('.')
        orbit = parts[0]
        date = parts[1]
        time = parts[2]

        year = date[0:4]
        month = date[4:6]

        fname1 = 'GPM2Ku5_uw3_'+date+'.'+time+'_to_*_'+orbit+'_NAM.nc'
        fname2 = glob.glob(indir+'/'+year+'/'+month+'/'+fname1)
        fname_full = fname2[0]
        fname_net4 = os.path.basename(fname_full)
        f_net4.write(fname_net4+"\n")
        fname_zeb = 'GPM2Ku5_uw3_'+date+'.'+time+'_'+orbit+'_NAM.nc'
        f_zeb.write(fname_zeb+"\n")

f_net4.close()
f_zeb.close()


