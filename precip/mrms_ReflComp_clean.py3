#!/usr/bin/python3

import os
import shutil
import glob

baseDir = '/home/disk/monsoon/precip/raw/radar/2021/mrms/MergedReflectivityQCComposite'

for date in os.listdir(baseDir):
    if date.startswith('20210') and os.path.isdir(baseDir+'/'+date):
        print('date =',date)
        os.chdir(baseDir+'/'+date)
        os.mkdir('10')
        flist = glob.glob('MRMS_MergedReflectivityQCComposite_*_20210???-???6??.*')
        for file in flist:
            shutil.move(file,'10/'+file)
        if len(os.listdir(baseDir+'/'+date+'/10')) == 144:
            extraFiles = glob.glob('MRMS*')
            for file in extraFiles:
                os.remove(file)
            goodFiles = glob.glob('10/*')
            for file in goodFiles:
                shutil.move(file,baseDir+'/'+date)
            os.rmdir('10')
        else:
            print('Need to check date',date)
