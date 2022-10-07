import os

indir = "/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/08"
outFile = '/home/disk/bob/gpm/nam_ku/classify/ex_data/2014/makeLinks.201408.csh'

fid = open(outFile,'w')

for file in os.listdir(indir):
    if file.startswith('GPM2Ku4_uw1'):
        print(file)
        linkName = file[0:27]+'.nc'
        fid.write("{}{} {}\n".format('ln -s ../2014/08/',file,linkName))

fid.close()
