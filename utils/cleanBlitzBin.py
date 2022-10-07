import os
import logging as log
import numpy as np

indir = "/home/disk/blitz/bin"

flist = "flist.blitzbin"

log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

f = open(indir+'/'+flist, 'r')
for binary in f:
    binary = binary.strip()
    os.rename(indir+'/'+binary,indir+'/TOREMOVE/'+binary)

f.close()

