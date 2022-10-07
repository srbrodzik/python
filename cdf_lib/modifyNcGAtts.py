# Modify global attribute in existing netcdf4 file

import numpy as np
import netCDF4 as nc4
import os

#--------------------- START INPUTS ---------------------

#inDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/sur_1km_cf/zdr'
inDir = '/home/disk/mjo/dynamo/data.server/interp/QCed/spolka/rain_rate_sur'
gattName1 = 'source'
gattVal1 = 'SPolKa radar data'
gattName2 = 'title'
gattVal2 = 'Rainrate values: minimum, expected, maximum'
gattName3 = 'references'
gattVal3 = 'Code used https://github.com/CSU-Radarmet/CSU_RadarTools'
gattName4 = 'comment'
gattVal4 = 'Based on 2.5km level of interpolated radar data'

#dates = ['20111001','20111002','20111003','20111004','20111005',
#         '20111006','20111007','20111008','20111009','20111010',
#         '20111011','20111012','20111013','20111014','20111015',
#         '20111016','20111017','20111018','20111019','20111020',
#         '20111021','20111022','20111023','20111024','20111025',
#         '20111026','20111027','20111028','20111029','20111030',
#         '20111031',
#         '20111101','20111102','20111103','20111104','20111105',
#         '20111106','20111107','20111108','20111109','20111110',
#         '20111111','20111112','20111113','20111114','20111115',
#         '20111116','20111117','20111118','20111119','20111120',
#         '20111121','20111122','20111123','20111124','20111125',
#         '20111126','20111127','20111128','20111129','20111130',
#         '20111201','20111202','20111203','20111204','20111205',
#         '20111206','20111207','20111208','20111209','20111210',
#         '20111211','20111212','20111213','20111214','20111215',
#         '20111216','20111217','20111218','20111219','20111220',
#         '20111221','20111222','20111223','20111224','20111225',
#         '20111226','20111227','20111228','20111229','20111230',
#         '20111231',
#         '20120101','20120102','20120103','20120104','20120105',
#         '20120106','20120107','20120108','20120109','20120110',
#         '20120111','20120112','20120113','20120114','20120115',
#         '20120116']
dates = ['20111031']

#---------------------- END INPUTS ----------------------

for date in dates:
    print date

    for infile in os.listdir(inDir+'/'+date):
        if infile.endswith('nc'):
            print infile
            ncid = nc4.Dataset(inDir+'/'+date+'/'+infile,'a')
            setattr(ncid,gattName1,gattVal1)
            setattr(ncid,gattName2,gattVal2)
            setattr(ncid,gattName3,gattVal3)
            setattr(ncid,gattName4,gattVal4)
            ncid.close()
