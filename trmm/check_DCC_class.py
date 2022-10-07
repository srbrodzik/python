import os
import numpy as np

# ---------- Start of Inputs ----------
types = ['str','mod']

regions = ['afc','cio','epo','h01','h02','h03','h04','h05','h06','h07','h08','ind','msa','pac','sam','tra','usa','wmp']

main_dir = '/home/disk/bob/trmm_v7'
sub_dir = ['africa','cent_ind_ocean','east_pac','hole_regions/Hole01','hole_regions/Hole02','hole_regions/Hole03',
          'hole_regions/Hole04','hole_regions/Hole05','hole_regions/Hole06','hole_regions/Hole07','hole_regions/Hole08',
          'india_main_v7','mesoAmer','pacific_ocean','southAmer','trop_atl','usa','warm_pool']

outdir = '/home/disk/bob/trmm_v7/check_dcc'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
prefix = 'Deep_Convective'
# ---------- End of Inputs ----------

M_TO_KM = 0.001
THRES1 = 4.0
THRES2 = 5.0

num_types = len(types)
num_regions = len(regions)
num_months = len(months)

for tt in range(0,num_types):
    print 'type = ',types[tt]
    if types[tt] == 'str':
        rest_of_path = 'classify/class_data/monthly_class_v6'
    else:
        rest_of_path = 'classify/class_data/monthly_class_v61'

    for rr in range(0,num_regions):
        print 'region = ',regions[rr]

        f_out = open(outdir+'/check_DCC_class_'+regions[rr]+'_'+types[tt]+'.out','w')
        f_out.write('{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('orbitDateTime','lon','lat','top','bot','xdim','ydim','area','terrHt','bot-terrHt'))

        for mm in range(0,num_months):
            #print mm
            path = main_dir+'/'+sub_dir[rr]+'/'+rest_of_path+'/'+months[mm]
            print 'path = ',path
            #for fname in os.listdir(in_dir+'/'+months[mm]):
            for fname in os.listdir(path):
                if prefix in fname and fname.endswith('dat'):
                    print fname
                    fname_info = os.path.splitext(fname)[0]+'.info'

                    # read info file
                    f_info = open(path+'/'+fname_info,'r')
                    dtg = f_info.readlines()
                    f_info.close()
                    
                    # read ascii dat file
                    #f = open(indir+'/'+months[mm]+'/'+fname,'r')
                    f = open(path+'/'+fname,'r')
                    index = 0
                    for line in f:
                        #print(repr(line))
                        line = line.strip()
                        columns = line.split()
                        lon = float(columns[0])
                        lat = float(columns[1])
                        area = float(columns[2])
                        top_ht = float(columns[3])
                        bot_ht = float(columns[4])
                        x_dim = float(columns[5])+0.05
                        y_dim = float(columns[6])+0.05
                        terr_ht = float(columns[7]) * M_TO_KM
                        diff = bot_ht - terr_ht
                        if diff < THRES2:
                            diff_flag = 'XXXXX'
                        else:
                            diff_flag = ''
                        if diff > THRES1:
                            dateTimeOrbit = dtg[index].strip()
                            f_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(dateTimeOrbit,lon,lat,top_ht,bot_ht,x_dim,y_dim,area,terr_ht,diff,diff_flag))
                        index+=1
                    f.close
            f_out.write('\n')
        f_out.close()
    
                
