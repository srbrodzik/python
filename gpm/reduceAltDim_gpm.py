#!/usr/bin/python

import os
import sys

indir = '/home/disk/bob/gpm/nam_imp_ku/classify/ex_data_v07/2023/02'
outdir = '/home/disk/bob/impacts/raw/gpm_ku/202302'
#indir = '/home/disk/shear2/brodzik/gpm/PMMproposal'
#outdir = '/home/disk/shear2/brodzik/gpm/PMMproposal/reducedAlts'

for file in os.listdir(indir):
    if file.endswith('.nc'):
        print file
        cmd = 'ncks -d alt,0,124 '+indir+'/'+file+' '+outdir+'/'+file
        os.system(cmd)
