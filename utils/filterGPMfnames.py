import os
import re

indir = "/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/05"
outFile = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/201405_bad_times.txt'

fid = open(outFile,'w')

for file in os.listdir(indir):
    print(file)
    phrase = file.replace('_',' ').replace('.',' ').split()
    startTime = phrase[3]
    endTime = phrase[6]
    if endTime < startTime:
        fid.write("{} {} {}\n".format(file,startTime,endTime))

fid.close()
