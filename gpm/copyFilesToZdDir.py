import os

ncdf_dir = '/home/disk/bob/gpm/asia_ku/classify/ex_data_v05'
zeb_dir = '/home/disk/bob/gpm/asia_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str'
flist = 'fileList'

#ncdf_dir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
#zeb_dir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str_mask'
#flist = 'fileList'

#ncdf_dir = '/home/disk/bob/gpm/asia_ku/classify/ex_data_v05'
#zeb_dir = '/home/disk/bob/gpm/asia_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dwc_str'
#flist = 'fileList'

#ncdf_dir = '/home/disk/bob/gpm/h02_ku/classify/ex_data_v05'
#zeb_dir = '/home/disk/bob/gpm/h02_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_mod'
#flist = 'fileList'

#ncdf_dir = '/home/disk/bob/gpm/aka_ku/classify/ex_data_v05'
#zeb_dir = '/home/disk/bob/gpm/aka_ku/zeb-data/gpm_ku_v05_uw/netcdf4_dcc_str'
#flist = 'fileList'

#ncdf_dir = '/home/disk/bob/gpm/eur_ku/classify/ex_data_v05'
#zeb_dir = '/home/disk/bob/gpm/eur_ku/zeb-data/gpm_ku_v05_uw/netcdf4_wcc_str'
#flist = 'fileList'

os.chdir(zeb_dir)
fid = open(flist,'r')
for line in fid:
    line = line.strip()
    print(repr(line))
    tmp = line.split('_')
    datetime = tmp[2]
    year = datetime[0:4]
    month = datetime[4:6]
    short_filename = tmp[0]+'_'+tmp[1]+'_'+tmp[2]+'_'+tmp[5]+'_'+tmp[6]
    command = 'cp '+ncdf_dir+'/'+year+'/'+month+'/'+line+' '+short_filename
    os.system(command)
fid.close()


