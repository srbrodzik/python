#!/usr/bin/python3

import sys
import os
import shutil

if len(sys.argv) != 2:
    print('Usage: {} [indir]'.format(sys.argv[0]))
    sys.exit()
else:
    indir = sys.argv[1]

#indir = '/home/wrf1/data/201504'

for dir in os.listdir(indir):
    if os.path.isdir(indir+'/'+dir):
        print(dir)
        runDir = indir+'/'+dir
        #os.chdir(runDir)
        if os.path.isdir(runDir+'/zipped'):
            zipDir = runDir+'/zipped'
            for file in os.listdir(zipDir):
                print(file)
                shutil.move(zipDir+'/'+file,
                            runDir+'/'+file)
            os.rmdir(zipDir)
