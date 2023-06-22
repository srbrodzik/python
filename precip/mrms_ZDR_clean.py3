#!/usr/bin/python3

import os
import shutil
import glob

baseDir = '/home/disk/monsoon/precip/raw/radar/2021/mrms/MergedZdr'

for date in os.listdir(baseDir):
    if date.startswith('20210') and os.path.isdir(baseDir+'/'+date):
        print('date =',date)
        os.chdir(baseDir+'/'+date)
        os.mkdir('10')
        flist = glob.glob('MRMS_MergedZdr_19.00_20210???-???5??.*')
        print('len(flist) =',len(flist))
        if len(flist) >= 144:
            flistAll = glob.glob('MRMS_MergedZdr_*_20210???-???5??.*')   
            for file in flistAll:
                shutil.move(file,'10/'+file)
            extraFiles = glob.glob('MRMS*')
            for file in extraFiles:
                os.remove(file)
            goodFiles = glob.glob('10/*')
            for file in goodFiles:
                shutil.move(file,baseDir+'/'+date)
            os.rmdir('10')
        else:
            flistAll = glob.glob('MRMS_MergedZdr_*_20210???-???5??.*')
            for file in flistAll:
                shutil.move(file,'10/'+file)
            print('Need to check date',date)
