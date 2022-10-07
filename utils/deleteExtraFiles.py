import os
import logging as log
import numpy as np

indir = "/home/disk/blitz/bin/TOREMOVE"
lrosedir = "/home/disk/meso-home/meso/lrose/bin"

#log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

for file in os.listdir(indir):
    print(file)
    if os.path.isfile(lrosedir+'/'+file):
        os.remove(indir+'/'+file)
