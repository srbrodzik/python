import numpy as np
import os

indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05_testing/2016/09'

for fname in os.listdir(indir):

    if fname.endswith('NAM.nc'):

        print 'fname = ',fname

        # get info from input fname
        parts = fname.split('_')
        fname_out = parts[0]+'_'+parts[1]+'_'+parts[2]+'_'+parts[5]+'_'+parts[6]
        print 'fname_out = ',fname_out

        # rename file
        os.rename(indir+'/'+fname,indir+'/'+fname_out)

