import os
import re

indir = "/home/disk/bob/gpm/asia_ku/classify/ex_data/2014/07"
outFile = '/home/disk/bob/gpm/asia_ku/classify/ex_data/2014/201407_bad_times.txt'

fid = open(outFile,'w')

for file in os.listdir(indir):
    print(file)
    phrase = file.replace('_',' ').replace('.',' ').split()
    startDay = phrase[2]
    endDay = phrase[5]
    startTime = phrase[3]
    endTime = phrase[6]
    if endTime < startTime and startDay == endDay:
        fid.write("{} {} {}\n".format(file,startTime,endTime))

fid.close()
