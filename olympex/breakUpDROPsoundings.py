import os
import numpy as np

indir = '/home/disk/funnel/olympex/archive/missions/dc8/sondes/QC'
infile = 'OLYMPEX_upaL4.0_drop'
outdir = indir

# read in all data & get indices for start of each sounding
in_id = open(indir+'/'+infile, 'r')
lineindex = -1
indices = np.array([],dtype=int)
data = np.array([],dtype=str)
for line in in_id:
    data = np.hstack((data, [line]))
    lineindex = lineindex + 1
    if line.startswith('STN',1,4):
        indices = np.hstack((indices, [lineindex]))
in_id.close()
numLines = data.size
numSoundings = indices.size

# write out all but last sounding to a separate file
for i in range(len(indices)-1):
    dateTime = str.split(data[indices[i]+1])
    outfile = 'drop.20'+dateTime[1]+'.'+dateTime[2]+'.txt'
    print outfile
    out_id = open(outdir+'/'+outfile, 'w')
    for j in range(indices[i],indices[i+1]):
        out_id.write(data[j])
    out_id.close()

# write out last sounding to a separate file
dateTime = str.split(data[indices[numSoundings-1]+1])
outfile = 'drop.20'+dateTime[1]+'.'+dateTime[2]+'.txt'
print outfile
out_id = open(outdir+'/'+outfile, 'w')
for j in range(indices[numSoundings-1],numLines):
    out_id.write(data[j])
out_id.close()
