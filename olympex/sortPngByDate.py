import numpy as np
import os
import logging as log

inDir = '/home/disk/funnel/olympex/archive/research/ts_snds'

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

os.chdir(inDir)
for fname in os.listdir(inDir):
    if fname.endswith('png'):
        log.info( "fname = {}".format(fname) )
        date = fname[17:25]
        if not os.path.exists(inDir+'/'+date):
            os.mkdir(inDir+'/'+date)
        os.rename(inDir+'/'+fname, inDir+'/'+date+'/'+fname)
