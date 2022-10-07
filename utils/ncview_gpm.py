import os

#inDir = '/home/disk/bob/gpm/sam_ku/classify/ex_data/2014'
inDir = '/home/disk/bob/gpm/asia_ku/classify/ex_data/2014'
month = '03'
listFile = '201403_maxLat_fnames.txt'

f = open(inDir+'/'+listFile,'r')
for line in f:
    filename = repr(line.strip())
    os.system("ncview "+inDir+'/'+month+'/'+filename)
