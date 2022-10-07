#!/usr/bin/python3

import os

paramDir = '/home/disk/monsoon/relampago/git/lrose-projects-relampago/projDir/ingest/params'
paramFile = 'RadxBufr.RMA1_archive'
baseDir = '/home/disk/monsoon/relampago/raw/radar/RMA1'
dates = ['20180208','20181213','20181214','20201024']
hours = ['00','01','02','03','04','05','06','07',
         '08','09','10','11','12','13','14','15',
         '16','17','18','19','20','21','22','23']

for date in dates:
    for hour in hours:
        currDir = baseDir+'/'+date+'/'+hour
        for vol in os.listdir(currDir):
            if os.path.isdir(currDir+'/'+vol):
                os.chdir(currDir+'/'+vol)
                cmd = 'RadxBufr -v -params '+paramDir+'/'+paramFile+' -f *.BUFR'
                os.system(cmd)
