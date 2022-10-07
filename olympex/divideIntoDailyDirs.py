import os
import numpy as np

baseDir = '/home/disk/funnel/olympex/archive/ops/text_sounding'
qcDir = baseDir+'/QC'

for fname in os.listdir(qcDir):
    if fname.endswith('txt'):
        print fname
        [cat,type,dateTime,site,imgType] = str.split(fname,'.')
        date = dateTime[0:8]
        newFname = cat+'.'+type+'.'+dateTime+'.'+site+'.'+imgType
        os.rename(qcDir+'/'+fname, baseDir+'/'+date+'/'+newFname)
        
