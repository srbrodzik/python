import os, os.path

main_path_a = '/home/disk/bob/olympex/cfradial/moments/xband_all/rhi_a'
main_path_b = '/home/disk/bob/olympex/cfradial/moments/xband_all/rhi_b'
main_path_c = '/home/disk/bob/olympex/cfradial/moments/xband_all/rhi_c'
main_path_d = '/home/disk/bob/olympex/cfradial/moments/xband_all/rhi_d'
main_path_e = '/home/disk/bob/olympex/cfradial/moments/xband_all/rhi_e'

dates = ['20151113','20151114','20151115',
         '20151116','20151117','20151118','20151119','20151120',
         '20151121','20151122','20151123','20151124','20151125',
         '20151127','20151128','20151129','20151130',
         '20151201','20151202','20151203','20151204','20151205',
         '20151206','20151207','20151208','20151209','20151210',
         '20151211','20151212','20151213','20151214','20151215',
         '20151216','20151217','20151218','20151219','20151220',
         '20151221','20151222','20151223','20151224','20151225',
         '20151226','20151227','20151228','20151229','20151230',
         '20151231',
         '20160101','20160102','20160103','20160104','20160105',
         '20160106','20160107','20160108','20160109','20160110',
         '20160111','20160112','20160113','20160114','20160115',
         '20160116','20160117','20160118','20160119','20160120',
         '20160121','20160122','20160123','20160124','20160125',
         '20160126','20160127','20160128','20160129','20160130',
         '20160131',
         '20160201','20160202','20160203','20160204','20160205',
         '20160206','20160207','20160208','20160209','20160210',
         '20160211','20160212','20160213','20160214','20160215',
         '20160216','20160217','20160218','20160219','20160220',
         '20160221','20160222','20160223','20160224','20160225',
         '20160226','20160227','20160228','20160229',
         '20160301','20160302','20160303','20160304','20160305',
         '20160306','20160307','20160308','20160309','20160310',
         '20160311','20160312','20160313','20160314','20160315',
         '20160316','20160317','20160318','20160319','20160320',
         '20160321','20160322','20160323','20160324','20160325',
         '20160326','20160327','20160328','20160329','20160330',
         '20160331',
         '20160401','20160420']

for date in dates:
    #print date
    path_a = main_path_a + '/' + date
    path_b = main_path_b + '/' + date
    path_c = main_path_c + '/' + date
    
    nfiles_a = len(next(os.walk(path_a))[2])
    nfiles_b = len(next(os.walk(path_b))[2])
    nfiles_c = len(next(os.walk(path_c))[2])
    #print 'nfiles_a = '+str(nfiles_a)
    #print 'nfiles_b = '+str(nfiles_b)
    #print 'nfiles_c = '+str(nfiles_c)

    if nfiles_a != nfiles_b or nfiles_a != nfiles_c or nfiles_b != nfiles_c:
        print date+': File counts are different'
    else:
        print date+': All nfiles = '+str(nfiles_a)
