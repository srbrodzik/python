#!/usr/bin/python3

import os

indir = '/home/disk/flood2/mm5/home/mm5rt/data'

f = open('/home/disk/meso-home/brodzik/misc/mm5counts_test.txt','w')

for date in os.listdir(indir):
    print(date)
    if date.startswith('1999') or date.startswith('2000'):
        os.chdir(indir+'/'+date)
        #f.write(date,'   ',len([name for name in os.listdir('.') if os.path.isfile(name)]) )
        f.write('{0} {1}\n'.format(date,len([name for name in os.listdir('.') if os.path.isfile(name)]) ))
        
f.close()
