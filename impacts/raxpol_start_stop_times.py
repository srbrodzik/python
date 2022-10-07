#!/usr/bin/python3

import os
import glob

rawdir = '/home/disk/bob/impacts/raw/radar/raxpol/20220224'
indir = '/home/disk/bob/impacts/cfradial/raxpol/sur/20220224'
outdir = indir
outfile = 'times.txt'

f = open(outdir+'/'+outfile,'w')

for file in os.listdir(indir):
    if file.endswith('nc'):
        print('file = {}'.format(file))
        (cfrad,start,stop,junk1,junk2,ext) = file.split('.')
        (start_date,start_time) = start.split('_')
        (junk1,junk2,stop_date,stop_time) = stop.split('_')
        rawfiles = glob.glob(rawdir+'/RAXPOL-'+start_date+'-'+start_time+'*Z.nc')
        if len(rawfiles) == 1:
            (path,rawfile) = os.path.split(rawfiles[0])
        else:
            rawfile = 'unknown'
        f.write('rawfile = {}\tcfradial: start = {} and end = {}\n'.format(rawfile,start_time,stop_time))

f.close()
    
