import os

indir = '/home/disk/bob/gpm/eur_ku/classify/class_data_v05_uw/monthly_class_v10m'
#months = ['01','02','03','04','05','06','07','08','09','10','11','12']
months = ['12']

num_months = len(months)
for i in range(0,num_months):
    month = months[i]
    print 'month =',month
    for fname in os.listdir(indir+'/'+month):
        if fname.endswith('dat'):
            infile = indir+'/'+month+'/'+fname
            f_in = open(infile,'r')
            outfile = indir+'/'+month+'/'+fname+'.NEW'
            f_out = open(outfile,'w')
            for index,line in enumerate(f_in):
                #print (index,line)
                if index%2 == 0:
                    a = line.strip('\n')
                else:
                    outline = a + line
                    f_out.write(outline)
            f_in.close()
            f_out.close()
            os.rename(outfile,infile)
