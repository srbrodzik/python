# ***** reorganize rainbow data for Radx processing *****
# Author: Stacy Brodzik, University of Washington
# Date: November 3, 2016
# Description: 

#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function
import logging as log
import os

# ------------------------------------ START INPUTS ------------------------------------
## input and output directories
inDir = '/home/disk/bob/olympex/xband/rainbowFiles/fullDataSet'

#dates = ['20151113','20151114','20151115',
#         '20151116','20151117','20151118','20151119','20151120',
#         '20151121','20151122','20151123','20151124','20151125',
#         '20151127','20151128','20151129','20151130',
#         '20151201','20151202','20151203','20151204','20151205',
#         '20151206','20151207','20151208','20151209','20151210',
#         '20151211','20151212','20151213','20151214','20151215',
#         '20151216','20151217','20151218','20151219','20151220',
#         '20151221','20151222','20151223','20151224','20151225',
#         '20151226','20151227','20151228','20151229','20151230',
#         '20151231',
#         '20160101','20160102','20160103','20160104','20160105',
#         '20160106','20160107','20160108','20160109','20160110',
#         '20160111','20160112','20160113','20160114']
dates = ['20160115',
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

# ------------------------------------- END INPUTS -------------------------------------

# set up logging
#log.basicConfig(format='%(levelname)s:%(message)s',level=log.INFO)

# process data for each date
for date in dates:
    print date
    
    os.chdir(inDir+'/'+date+'/BIRDBATH')
    for file in os.listdir(inDir+'/'+date+'/BIRDBATH'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/BIRDBATH'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/BIRDBATH/png/'+swpFile)
                elif swpFile.endswith('azi'):
                    if 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/dBZ/'+swpFile)
                    elif 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/dBuZ/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/RhoHV/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/SNR/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/SQI/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/ZDR/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/BIRDBATH/PhiDP/'+swpFile)
            os.remove(base[0])
            
    os.chdir(inDir+'/'+date+'/DOPVOL1_A')
    for file in os.listdir(inDir+'/'+date+'/DOPVOL1_A'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/DOPVOL1_A'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/DOPVOL1_A/png/'+swpFile)
                elif swpFile.endswith('azi'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_A/KDP/'+swpFile)
            os.remove(base[0])
    
    os.chdir(inDir+'/'+date+'/DOPVOL1_B')
    for file in os.listdir(inDir+'/'+date+'/DOPVOL1_B'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/DOPVOL1_B'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/DOPVOL1_B/png/'+swpFile)
                elif swpFile.endswith('azi'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_B/KDP/'+swpFile)
            os.remove(base[0])

    os.chdir(inDir+'/'+date+'/DOPVOL1_C')
    for file in os.listdir(inDir+'/'+date+'/DOPVOL1_C'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/DOPVOL1_C'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/DOPVOL1_C/png/'+swpFile)
                elif swpFile.endswith('azi'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/DOPVOL1_C/KDP/'+swpFile)
            os.remove(base[0])
        
    os.chdir(inDir+'/'+date+'/RHI_A')
    for file in os.listdir(inDir+'/'+date+'/RHI_A'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/RHI_A'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/RHI_A/png/'+swpFile)
                elif swpFile.endswith('ele'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_A/KDP/'+swpFile)
            os.remove(base[0])
        
    os.chdir(inDir+'/'+date+'/RHI_B')
    for file in os.listdir(inDir+'/'+date+'/RHI_B'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/RHI_B'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/RHI_B/png/'+swpFile)
                elif swpFile.endswith('ele'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_B/KDP/'+swpFile)
            os.remove(base[0])

    os.chdir(inDir+'/'+date+'/RHI_C')
    for file in os.listdir(inDir+'/'+date+'/RHI_C'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/RHI_C'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/RHI_C/png/'+swpFile)
                elif swpFile.endswith('ele'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_C/KDP/'+swpFile)
            os.remove(base[0])

    os.chdir(inDir+'/'+date+'/RHI_D')
    for file in os.listdir(inDir+'/'+date+'/RHI_D'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/RHI_D'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/RHI_D/png/'+swpFile)
                elif swpFile.endswith('ele'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_D/KDP/'+swpFile)
            os.remove(base[0])

    os.chdir(inDir+'/'+date+'/RHI_E')
    for file in os.listdir(inDir+'/'+date+'/RHI_E'):
        if file.endswith('gz'):
            print file
            os.system('gunzip '+file)
            base = os.path.splitext(file)
            os.system('tar xvf '+base[0])
            for swpFile in os.listdir(inDir+'/'+date+'/RHI_E'):
                if swpFile.endswith('png'):
                    os.rename(swpFile, inDir+'/RHI_E/png/'+swpFile)
                elif swpFile.endswith('ele'):
                    if 'dBuZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/dBuZ/'+swpFile)
                    elif 'dBZ' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/dBZ/'+swpFile)
                    elif 'uPhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/uPhiDP/'+swpFile)
                    elif 'PhiDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/PhiDP/'+swpFile)
                    elif 'RhoHV' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/RhoHV/'+swpFile)
                    elif 'SQI' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/SQI/'+swpFile)
                    elif 'SNR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/SNR/'+swpFile)
                    elif 'V' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/V/'+swpFile)
                    elif 'W' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/W/'+swpFile)
                    elif 'ZDR' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/ZDR/'+swpFile)
                    elif 'uKDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/uKDP/'+swpFile)
                    elif 'KDP' in swpFile:
                        os.rename(swpFile, inDir+'/RHI_E/KDP/'+swpFile)
            os.remove(base[0])

    
