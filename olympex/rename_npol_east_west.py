import os

dates = ['20151113','20151114','20151115',
         '20151116','20151117','20151118','20151119','20151120',
         '20151121','20151122','20151123','20151124','20151125',
         '20151126','20151130',
         '20151201','20151202','20151203','20151204','20151205',
         '20151206','20151207','20151208','20151209','20151210',
         '20151211','20151212','20151213','20151214','20151215',
         '20151216','20151217','20151218','20151219',
         '20160103','20160104','20160105',
         '20160106','20160108','20160110',
         '20160111','20160112','20160113','20160114','20160115']
#dates = ['20151112']

indir_moments = '/home/disk/bob/olympex/cfradial/moments/npol_qc2/rhi'
indir_pid     = '/home/disk/bob/olympex/cfradial/partrain/npol_qc2/rhi'

for date in dates:
    for fname in os.listdir(indir_moments+'/'+date):
        if fname.endswith('nc'):
            print 'fname = ',fname
            parts = fname.split('_')
            base_name = parts[0]+'_'+parts[1]+'_'+parts[2]+'_'+parts[3]+'_'+parts[4]+'_'+parts[5]+'_'+parts[6]
            suffix = parts[7]
            os.rename(indir_pid+'/'+date+'/'+base_name+'.nc',indir_pid+'/'+date+'/'+fname)
            
