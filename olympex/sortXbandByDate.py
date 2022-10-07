import numpy as np
import os
import logging as log

inDir = '/home/disk/bob/olympex/xband/rainbowFiles/fullDataSet/RHI_A'
fields = ['KDP','PhiDP','RhoHV','SNR','SQI','V','W','ZDR','dBZ','dBuZ','uKDP','uPhiDP']
#suffix = 'azi'
suffix = 'ele'

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

#for field in fields:
#log.info( "field = {}".format(field) )
#os.chdir(inDir+'/'+field)
os.chdir(inDir)
#for fname in os.listdir(inDir+'/'+field):
for fname in os.listdir(inDir):
    if fname.startswith('2016') & fname.endswith(suffix):
    #if fname.endswith(suffix):
        log.info( "fname = {}".format(fname) )
        date = fname[:8]
        if not os.path.exists(inDir+'/'+date):
            os.mkdir(inDir+'/'+date)
        #os.rename(inDir+'/'+field+'/'+fname, inDir+'/'+date+'/'+fname)
        os.rename(inDir+'/'+fname, inDir+'/'+date+'/'+fname)
