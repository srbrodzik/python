#!/usr/bin/python3

import pyart
from chivo_qc_more import radar_qc_more
import os
import fnmatch

indir = '/home/disk/monsoon/relampago/cfradial/csu_moments_v1bf/sur/'
outdir = '/home/disk/monsoon/relampago/cfradial/csu_moments_v1bf_more/sur/'

for path,dirs,files in sorted(os.walk(indir)): 
        for file in sorted(files):
            if fnmatch.fnmatch(file,'*.nc'):
                fullname = os.path.join(path,file)
                new_radar = radar_qc_more(fullname)
                pyart.io.write_cfradial(outdir+file, new_radar)
                print ('done with ' + file)
