import os
import shutil

indir = '/home/disk/bob/relampago/raw/soundings/cordoba'

for file in os.listdir(indir):
    if file.endswith('.nc'):
        print file
        year = file[12:16]
        month = file[17:19]
        if not os.path.exists(indir+'/'+year):
            os.makedirs(indir+'/'+year)
        shutil.move(indir+'/'+file,indir+'/'+year)

