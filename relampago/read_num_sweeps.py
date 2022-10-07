import os
import netCDF4 as nc4
import numpy as np
import numpy.ma as ma

indir = '/home/disk/monsoon/relampago/cfradial/csu_chivo/sur'
dates = ['20181110','20181111','20181112','20181113',
         '20181115','20181116','20181117','20181118',
         '20181120','20181121','20181122','20181123','20181124',
         '20181125','20181126','20181127','20181128','20181129',
         '20181130','20181201','20181202','20181203','20181204',
         '20181205','20181206','20181208','20181209',
         '20181210','20181211','20181212','20181213','20181214']

for date in dates:

    outfile = indir+'/'+date+'/sweep_nums'
    f = open(outfile,'w')

    for fname in os.listdir(indir+'/'+date):
        if fname.endswith('nc'):
            f.write('fname = %s' % fname)
            ncid = nc4.Dataset(indir+'/'+date+'/'+fname,'r')
            num_sweeps = len(ncid.dimensions['sweep'])
            sweep_mode = ncid.variables['sweep_mode'][:]
            tmp = sweep_mode[0]
            mode = ''.join(ma.getdata(tmp))
            f.write('   num_sweeps = %d\t' % (num_sweeps) )
            f.write('   sweep_mode = %s\n' % (mode) )
            ncid.close()
            
    f.close()
