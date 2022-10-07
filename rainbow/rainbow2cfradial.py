# ***** rainbow to cfradial conversion code for Canadian xband data *****
# Author: Stacy Brodzik, University of Washington
# Date: October 14, 2016
# Description: 

#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function
import pyart
import logging as log
import os

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
inDir1 = '/home/disk/bob/olympex/xband/rainbowFiles/RAW_DOPVOL1_A/'
inDir2 = '/home/disk/bob/olympex/xband/rainbowFiles/RAW_DOPVOL1_B/'
inDir3 = '/home/disk/bob/olympex/xband/rainbowFiles/RAW_DOPVOL1_C/'

outDir = '/home/disk/bob/olympex/cfradial/moments/xband_dbz/sur_single'

dates = ['20151114']
#dates = ['20151114','20151115',
#         '20151116','20151117','20151118','20151119','20151120',
#         '20151121','20151122','20151123','20151124','20151125',
#         '20151126','20151127','20151128','20151129','20151130',
#         '20151201','20151202','20151203','20151204','20151205',
#         '20151206','20151207','20151208','20151209','20151210',
#         '20151211','20151212','20151213','20151214','20151215',
#         '20151216','20151217','20151218','20151219','20151220',
#         '20151221','20151222','20151223','20151224','20151225',
#         '20151226','20151227','20151228','20151229','20151230',
#         '20151231',
#         '20160101','20160102','20160103','20160104','20160105',
#         '20160106','20160107','20160108','20160109','20160110',
#         '20160111','20160112','20160113','20160114','20160115',
#         '20160116','20160117','20160118','20160119','20160120',
#`         '20160121','20160122','20160123','20160124','20160125']

field1 = 'dBZ'
field2 = 'V'
# ------------------------------------- END INPUTS -------------------------------------

# set up logging
#log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# process data for each date
for date in dates:
    print date
    for file in os.listdir(inDir1+field1):
        if file.startswith(date):
            if not os.path.exists(outDir+'/'+date):
                os.makedirs(outDir+'/'+date)
            print file
            time = file[8:14]
            print time
            prefix = date+time+'00'
            radar1 = pyart.aux_io.read_rainbow_wrl(inDir1+field1+'/'+prefix+field1+'.azi')
            radar2 = pyart.aux_io.read_rainbow_wrl(inDir2+field2+'/'+prefix+field2+'.azi')
            radar = pyart.util.join_radar(radar1,radar2)
            angle = radar.fixed_angle['data'][0]
            mode = radar.scan_type
            outFile = 'cfrad.'+date+'_'+time+'_xband_'+str(angle)+'_'+mode.upper()+'.nc'
            print outFile
            pyart.io.write_cfradial(outDir+'/'+date+'/'+outFile,radar)
            
    for file in os.listdir(inDir2):
        if file.startswith(date):
            if not os.path.exists(outDir+'/'+date):
                os.makedirs(outDir+'/'+date)
            print file
            time = file[8:14]
            print time
            radar = pyart.aux_io.read_rainbow_wrl(inDir2+'/'+file)
            angle = radar.fixed_angle['data'][0]
            mode = radar.scan_type
            outFile = 'cfrad.'+date+'_'+time+'_xband_'+field+'_'+str(angle)+'_'+mode.upper()+'.nc'
            print outFile
            pyart.io.write_cfradial(outDir+'/'+date+'/'+outFile,radar)
            
    for file in os.listdir(inDir3):
        if file.startswith(date):
            if not os.path.exists(outDir+'/'+date):
                os.makedirs(outDir+'/'+date)
            print file
            time = file[8:14]
            print time
            radar = pyart.aux_io.read_rainbow_wrl(inDir3+'/'+file)
            angle = radar.fixed_angle['data'][0]
            mode = radar.scan_type
            outFile = 'cfrad.'+date+'_'+time+'_xband_'+field+'_'+str(angle)+'_'+mode.upper()+'.nc'
            print outFile
            pyart.io.write_cfradial(outDir+'/'+date+'/'+outFile,radar)
            
