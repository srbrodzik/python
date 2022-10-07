import numpy as np
import os

ex_data_dir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
zeb_dir = '/home/disk/bob/gpm/nam_ku/zeb-data/gpm_ku_v5'

flist = zeb_dir+'/wcc_nasa.txt'

f = open(flist, 'r')
for line in f:
    line = line.strip()
    print(repr(line))
    tmp = line.split('.')
    orbit = tmp[0]
    datetime = tmp[1]+'.'+tmp[2]
    year = datetime[0:4]
    month = datetime[4:6]
    #print(orbit,datetime,year,month)
    name_in = ex_data_dir+'/'+year+'/'+month+'/GPM2Ku5_uw3_'+datetime+'_to_*_'+orbit+'_NAM.nc'
    name_out = zeb_dir+'/GPM2Ku5_uw3_'+datetime+'_'+orbit+'_NAM.nc'
    command = '/bin/cp '+name_in+' '+name_out
    os.system(command)

f.close()


