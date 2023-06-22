# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
# THIS IS CODE THAT WILL ATTRIBUTE BOTH NLDN AND LMA FLASHES TO YOUR STORMS
# AS LONG AS THEY ARE ORGANIZED BY LAT/LON


from lmatools_brody.flashsort.autosort.autorun import run_files_with_params, test_output, logger_setup
import numpy as np
from lmatools_brody.flashsort.autosort.autorun_sklearn import cluster
#from autorun_sklearn import cluster
import glob
import datetime 
import os
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import time as timemod
import lightning_tools as LT
import numpy.lib.recfunctions as rec
from Config import Config

base_path = '/data3/brfuchs/ALCOcases'

# first, let's try to open with yaml
radar_path = 'AL12mar2010'
conf = Config(path='%s/%s/lightningconfig.yaml' % (base_path, radar_path))


start_time = datetime.datetime.now()
print 'START: %s' %start_time

# since attributing flashes all the way up to next radar file, set to the ENTIRE time spacing
# between radar files

##### HERE ARE SOME TUNABLE OPTIONS FOR PROCESSING #########
mintime = 300. # need to define this!!
zipold = 0 # zip old files? only if not case study
untar = 1 # if need to untar older files first
start = 0 # testing purposes, which LMA file to start with in a day
# THESE ARE THINGS TO CHANGE IF YOU WANT TO
t_offset = 0. # flash/source attribution time offset in seconds-> default (0) is forward, negative for middle or backward
#year, month, date, n_days = 2012, 6, 6, 2
nldn_flag = 1 # if wanna attribute nldn data too
sdist = 10. # distance in km for attributing flashes to cells
################################################################




outpath = '/data2/brfuchs/research/tracking/countflash/'# just dummy directory for lma log output
f_dir = '/data2/brfuchs/research/radartest/lma_nc' # keep these as Brody's directories
l_dir = '/data2/brfuchs/research/lma/%s'%conf.v['region']
nldn_dir = '/net/longs/data2/brfuchs/research/nldn/%d'%conf.v['year'] # directory where NLDN flash files sit
logger_setup(outpath)
b_dir = '/net/longs/data2/brfuchs/research'
flashloop = 0 # LOOP THRU FLASHES AND MAKE PLOTS OF EACH ONE??

data = None

mintime_d = mintime/86400. # mintime in units of day
re = 6374.0
dlist = sorted(glob.glob('%s/%s/%d????'%(base_path, conf.v['radar_dir'],conf.v['year']))) 
# list of radar files (that have been cell tracked)
sdlist = [os.path.basename(sd) for sd in dlist] # list of directories with date names
t_data = None

reg_conf = Config(path='./%slightningconfig.yaml' % (conf.v['region']))


lma_params = {'stations':(reg_conf.v['s_min'],14), 'mask_length': reg_conf.v['ml'], 
          'chi2':(0,2.0),
          'ascii_flashes_out':'flashes_out.dat',
          'ctr_lat':reg_conf.v['ctrLat'], 'ctr_lon':reg_conf.v['ctrLon'],
          }


#os.system('tar -zvxf 20110804/080411_nldn.tar.gz -C 20110804')
print 'directory: %s'%conf.v['radar_dir']



# CREATE JULIAN DAY ARRAY OF RADAR FILE TIMES
if conf.v['file_type'] == 'chill':
    files = sorted(glob.glob('%s/%s/%s*/*.%s' %(base_path, conf.v['radar_dir'],conf.v['year'],conf.v['file_ext']))) 
    # ALL radar files FOR ALL DAYS
elif conf.v['file_type'] == 'nmq':
    files = sorted(glob.glob('%s/%s/%s/%s*/*.%s' %(base_path, conf.v['radar_dir'],conf.v['year'],conf.v['file_ext']))) 
    # ALL radar files FOR ALL DAYS
	# this looks like it might be the same as chill so can combine


rl = len(conf.v['radar_name'])



# GONNA HAVE TO GET RID OF FDATES STUFF BECAUSE NOT EVERYONE HAS IT     
jul_days = []
for i in range(len(files)):
    fname = os.path.basename(files[i])
    if conf.v['file_type'] == 'chill':
    	yyyy,mm,dd,hh,min,ss = int(fname[rl:rl+4]),int(fname[rl+4:rl+6]),int(fname[rl+6:rl+8]),\
		int(fname[rl+9:rl+11]),int(fname[rl+11:rl+13]),int(fname[rl+13:rl+15])
    elif conf.v['file_type'] == 'nmq':
	yyyy,mm,dd,hh,min,ss = int(fname[12:16]),int(fname[16:18]),\
		int(fname[18:20]),int(fname[21:23]),int(fname[23:25]),0
    to_date = datetime.datetime(yyyy,mm,dd,hh,min,ss) # convert to datetime object
    jul_day = datetime.datetime.timetuple(to_date) # gives julian day of year
    jul_days.append(jul_day[7]+jul_day[3]/24.0+jul_day[4]/1440.0+jul_day[5]/86400.0)
jul_days = np.array(jul_days)
rflag = np.zeros(len(files),np.int) # array that has write flag for each file in a day
unique_jul_days = np.unique(np.floor(jul_days))

try:
    if start == 0: start = input('What file do you want to start on? Files are %s seconds [0]: '%(conf.v['file_duration']))
except SyntaxError:
    start = 0



print 'Starting at %d'%start

for iday in range(len(sdlist[:conf.v['ndays']])): # master loop thru each day
    print 'day is', sdlist[iday]

    file_yr = sdlist[iday][0:4]
    file_mo = sdlist[iday][4:6]
    file_day = sdlist[iday][6:8]


    tar_file = '%s/%s/%s%s%s/tracked.tar.gz' %(base_path, conf.v['radar_dir'], file_yr, file_mo, file_day)


    if untar == 1:
        print 'INITIAL UNTARRING OF ALL RADAR FILES....'
            # CHANGE THIS IF NEEDED
        os.system('tar -zxvf %s -C %s/%s/%s%s%s'
               %(tar_file, base_path, conf.v['radar_dir'], file_yr, file_mo, file_day))

        # NOW NEED COPYING ROUTINE HERE
    #files = sorted(glob.glob('%s/%s/CHL*.cdf'%(radar_dir,sdlist[iday]))) # commented b/c defined files above over
    # multiple days

    # NOW TO SET UP TIME ARRAY FOR EACH DAY
    mint = 0.
    maxt = 86400.
    nd2 = len(files) # nd2 is now TOTAL number of files (over multiple days)
    rsec = np.zeros(nd2)
    mnr = np.zeros(nd2)
    tt = np.zeros(nd2)
    rstart = np.zeros(nd2)
    rend = np.zeros(nd2)

    """for i in range(nd2):
        fname = os.path.basename(files[i])
        hh = float(fname[12:14])
        mm = float(fname[14:16])
        ss = float(fname[16:18])
        rsec[i] = hh*3600.0+60.0*mm+ss"""
    for i in range(nd2):
        rsec[i] = jul_days[i] # simply grab the jul_days for rsec
    for i in range(1,nd2):
        mnr[i] = (rsec[i]+rsec[i-1])/2.0
    for i in range(1,nd2-1):
        tt[i] = mnr[i+1]-mnr[i]
        rstart[i] = mnr[i]
        rend[i] = mnr[i+1]
    tt[0] = mnr[1] - mint
    rstart[0] = mint
    rend[0] = mnr[1]
    rstart[nd2-1] = maxt-mnr[nd2-1]
    rend[nd2-1] = maxt
    rdur = tt



# NOW READ IN DATA TO BE FLASH PROCESSED

    # this is currently set up so LMA raw files for each day are in separate directories -
    # so can use these lma files presumably for any and all cases
    l_files = sorted(glob.glob('%s/%d/%s/LYLOUT*%s.dat*' %(l_dir, conf.v['year'], sdlist[iday], conf.v['file_duration']))) 
		# still only grabbing LMA files for one day


    for ifi, lfile in enumerate(l_files[start:]): # looping thru each LMA file in a day

        print lfile
	n_src = LT.n_sources(lfile) # check how many sources are in the file
	ml = LT.mask_length(lfile) # figure out the mask length in the file
	lma_params['mask_length'] = ml # set the mask length for flash counting purposes


	if n_src > conf.v['src_thresh']:
# now have n_src and mask length, can do flash counting 
    	    flash_points, src = LT.flash_count(lfile, lma_params, min_points=reg_conf.v['min_points'])
	
            n_fl = flash_points['lat'].shape[0] # number of flashes?
            #print 'NUMBER OF FLASHES?? (n_fl): %i' %len(n_fl)
            n_src = src['lat'].shape[0]
	    print 'source_points.shape: %d'%n_src


# JULIAN DAY CONVERSIONS HERE TO ATTRIBUTE INDIVIDUAL FLASHES AND SOURCES TO CORRECT RADAR FILES	

            # flash file attribution
            iradar = LT.file_attribute(flash_points['time'], sdlist[iday], jul_days, 
			time_offset=t_offset, mintime=mintime) 
            # t_offset (in seconds) to change how sources, flashes are attributed to files, i.e forward or center or backward

            good = iradar != -1 # finding good attributed sources

            iradar = iradar[good]

            flash = LT.filter_dict(flash_points, good) # flashes with enough points and attributed to radar files
            print 'flash.shape: ', flash['lat'].shape[0]

            # source file attribution
            s_iradar = LT.file_attribute(src['time'], sdlist[iday], jul_days, time_offset=t_offset)

            s_good = s_iradar != -1
            src = LT.filter_dict(src, s_good)
            s_iradar = s_iradar[s_good]

            array = iradar[iradar.sort()] # just sorting the array
            s_array = s_iradar[s_iradar.sort()]

            wh = np.unique(array) # indices representing unique radar files associated with LMA file
            s_wh = np.unique(s_array)
            cnt = wh.shape[0] # size of wh array
            s_cnt = s_wh.shape[0]

# loop thru each radar file
            for i in range(cnt):
                if wh[i] != -1: # only interested in iradar indices not equal to -1
                    print 'LOOPING THRU EACH RADAR FILE IN LMA FILE'
                    cin = np.where(iradar == wh[i]) # indices of all flashes to be attributed to a given radar file
                    s_cin = np.where(s_iradar == s_wh[i])
                    ccnt = cin[0].shape[0]
                    s_ccnt = s_cin[0].shape[0]
# FLASHES
                    this_fl = LT.filter_dict(flash, cin) # flashes in this radar file
                    print 'this_fl.shape: ', this_fl['lat'].shape[0]
# SOURCES
                    this_src = LT.filter_dict(src, s_cin)

                    r_file = files[wh[i]] # select the radar file to write to
# OPENING NETCDF FILE
                    print 'OPENING:', os.path.basename(r_file)
                    ncid = Dataset(r_file)
                    cell = ncid.variables['cell'][:]
                    cell = np.ma.getdata(cell)
                    print 'cell max ', cell.max()
                    if cell.max() > 0:	
			# THESE ARE THE ATTRIBUTES FROM NMQ FILES
			if conf.v['file_type'] == 'nmq':
                            lonfirst = getattr(ncid, 'Loncut')
                            latfirst = getattr(ncid, 'Latcut')
                            londelta = getattr(ncid, 'LonGridSpacing')
                            latdelta = getattr(ncid, 'LatGridSpacing')
                            nlon = cell.shape[1]
                            nlat = cell.shape[0]
			# THESE ARE FOR NEW CHL FILES
			elif conf.v['file_type'] == 'chill':
			    lat_array = ncid.variables['LAT'][:]
			    lon_array = ncid.variables['LON'][:]
			
			    latfirst, lonfirst = lat_array[0], lon_array[0]
			    latdelta = np.average(lat_array[1:] - lat_array[:-1])
			    londelta = np.average(lon_array[1:] - lon_array[:-1])
			    nlat, nlon = lat_array.shape[0], lon_array.shape[0]


                    # HERE TAKING OUT THE RAW CODE AND SUBBING IN THE FUNCTION TO SEE IF IT WORKS

                        this_fl['feat'] = LT.cell_attribute(this_fl['initlat'], this_fl['initlon'], cell, latfirst, \
                    		lonfirst, reg_conf.v['ctrLat'], reg_conf.v['ctrLon'], nlat, nlon, latdelta=latdelta, londelta=londelta, \
                    		latincreasing=conf.v['latinc'], sdist=10.0)

                        this_src['feat'] = LT.cell_attribute(this_src['lat'], this_src['lon'], cell, latfirst, \
                            lonfirst, reg_conf.v['ctrLat'], reg_conf.v['ctrLon'], nlat, nlon, latdelta=latdelta, londelta=londelta, \
                             latincreasing=conf.v['latinc'], sdist=10.0)

                    # END OF SOURCE ATTRIBUTING, NOW NEED TO DO SOME FILTERING ON EACH

                        fin = np.where(this_fl['feat'] > 0)
                        s_fin = np.where(this_src['feat'] > 0)
                        tj = fin[0].shape[0]
                        s_tj = s_fin[0].shape[0]
                        print 'NUM OF FLASHES ATTRIBUTED TO CELLS: ', tj
                        print 'NUM OF SOURCES ATTRIBUTED TO CELLS: ', s_tj

                        tk = tj
                        s_tk = s_tj
# FLASHES
                        #if (tj > 0):
                            # flash dictionary that has survived this far and will be written				
                        #    new_flash = LT.filter_dict(this_fl, fin)
                        #    print 'write_flash.shape: ', new_flash['lat'].shape[0]
			new_flash = LT.filter_dict(this_fl, fin)
			new_src = LT.filter_dict(this_src, s_fin)
# SOURCES
                        #if (s_tj > 0):
                        #    new_src = LT.filter_dict(this_src, s_fin)

                        if (rflag[wh][i] == 1):
                            print 'RADAR FILE WRITTEN TO BEFORE'
                        f_flag = 0
                        try:
    # FLASHES
                            old_flash = {}

                            for key in this_fl.keys():
                                old_flash[key] = ncid.variables['lma_flash_%s'%key][:]

                            f_flag = 1
                            max_id = np.amax(old_flash['id'])

                            # writing old nc variables to new ones being written out
                            if (tj > 0):
                                new_flash['id'] = new_flash['id'] + old_flash['id'].max() + 1 
                                write_flash = LT.dict_concatenate(old_flash, new_flash)
                                    # ADDING MAX_ID TO ALL SOURCES FOR NO OVERLAP

                        except KeyError:
                            max_id = 0
                            write_flash = new_flash # if no old flashes, just write the new flashes
				# but what if no flashes but there are sources, new_flash isnt a thing
                            pass
    # NOW SOURCES
                        try:
                            old_src = {}

                            for key in this_src.keys(): # grabbing any existing sources
                                old_src[key] = ncid.variables['lma_src_%s'%key][:]

                            if (s_tj > 0):
                                new_src['id'] = new_src['id'] + old_flash['id'].max() + 1 
                                write_src = LT.dict_concatenate(old_src, new_src)

    # IF SOMETHING THERE IN EITHER FLASHES OR SOURCES, THEN COPY THE FILE OVER AGAIN	

                        except KeyError:
                            write_src = new_src
#			    if (f_flag == 1):
                        if (rflag[wh][i] == 1):
                            #os.system('tar -xvf %s%s/tracked.tar.gz -C %s%s %s'
                                    #%(radar_dir,sdlist[iday],radar_dir,sdlist[iday],os.path.basename(files[wh[i]])))
                            file_base = os.path.basename(files[wh[i]])
                            #                                print 'RE-UNTARRING FILE %s FOR ADDITIONAL WRITING' %file_base
                            os.system('tar -xvf %s -C %s/%s/%s%s%s %s'
                                    %(tar_file, base_path, conf.v['radar_dir'], file_yr, file_mo, file_day, file_base))
                            #os.system('mv %s %s'%(files[wh[i]], sdlist[iday]))
                            #copy file over again to rewrite bigger array
                        ncid.close()

                       # NOW OPEN THE NEW NETCDF FOR WRITING OF SOURCES AND FLASHES 

                        ncwid = Dataset(r_file,'a')
                        if tj > 0:
                            ncwid.createDimension('f', write_flash['lat'].shape[0])
			    ncwid.createDimension('dummy', 1)
                            # WRITING THE FLASH INFO
                            # MAKE FUNCTION THAT WILL LOOP OVER THE DICTIONARY??

                            LT.write_nc_dict(write_flash, ncwid, 'f', 'd', prefix='lma_flash_')
			# now writing out minimum number of stations so info is preserved down the road
			    LT.write_nc_variable(ncwid, reg_conf.v['s_min'], 'min_stn', 'dummy', 'd')  

                            # NOW WRITING THE SOURCE INFO

                            if s_tj > 0:
                                ncwid.createDimension('s', write_src['lat'].shape[0])
                                LT.write_nc_dict(write_src, ncwid, 's', 'd', prefix='lma_src_')


                        ncwid.close() 
                        if ((s_tj != 0) | (tj !=0)):
                            rflag[wh[i]] = 1

                    else:
                        print 'NO ACTIVE CELLS IN FILE'


        else: # trying flash counting
            print 'NO FLASHES'
	    timemod.sleep(0.1)
            pass

    # reset start to 0 incase going to next day, probably need to start with first file
    start = 0

    # NOW A QUICK LITTLE NLDN ATTRIBUTION
    if nldn_flag:
        print 'NOW DOING NLDN ATTRIBUTION'
        nldn_day = unique_jul_days[iday]
        print 'OPENING %s'%nldn_day
        nldn_file = 'Nflash%d.%03d_daily_v3_lit.raw'%(conf.v['year'], nldn_day)

        nldn_lines = np.genfromtxt('%s/%s'%(nldn_dir, nldn_file), dtype = [('date', 'S10'), ('time', 'S15'), ('lat', 'f8'), \
                    ('lon', 'f8'), ('ss', 'f8'), ('cg', 'S5')], usecols = (0,1,2,3,4,13))

        nldn_dts = np.array([ datetime.datetime(2000+int(n[6:8]), int(n[0:2]), int(n[3:5]), int(t[0:2]), int(t[3:5]), \
                            int(t[6:8]), 1000*int(t[9:12])) for n, t in zip(nldn_lines['date'], nldn_lines['time']) ] )
            # the last entry in the datetime is microseconds which is why it's multiplied by 1000
        nldn_tuples = np.array( [datetime.datetime.timetuple(nd) for nd in nldn_dts] ) # gives julian day of year
        nldn_times = np.array( [j[3]*3600.0+j[4]*60.0+j[5]+ndt.microsecond/1.e6 for j, ndt in \
                            zip(nldn_tuples, nldn_dts)] ) # want this in day seconds
    # j[3] gives hours, j[4] gives minutes, j[5] seconds, had to get microseconds from datetime objects cuz not in timetuples

            # this file_attribute function does not need the day because it handles it by itself
        inldn = LT.file_attribute(nldn_times, sdlist[iday], jul_days, time_offset=t_offset, mintime=mintime)

        unldn = np.unique(inldn)
        if unldn.shape[0] > 1: # means there are actually some files to try attributing to
            for nl in unldn[1:]: # now loop thru the radar files

                this_time = np.where(inldn == nl) # which of the day's nldn flashes are in this time
                nldn_lines = rec.drop_fields(nldn_lines, ('date', 'time')) # dropping date and time fields because they're strings
                nldn_lines = rec.append_fields(nldn_lines, 'time', nldn_times, usemask = False)
                n_flashes = nldn_lines[this_time]
                nfile = files[nl]
                print 'OPENING:', os.path.basename(nfile)
                ncid = Dataset(nfile)
                cell = ncid.variables['cell'][:]
                cell = np.ma.getdata(cell)
                print 'cell max ', cell.max()
                if cell.max() > 0:	
                    # THESE ARE THE ATTRIBUTES FROM NMQ FILES
		    if conf.v['file_type'] == 'nmq':
                    	lonfirst = getattr(ncid, 'Loncut')
                    	latfirst = getattr(ncid, 'Latcut')
                    	londelta = getattr(ncid, 'LonGridSpacing')
                    	latdelta = getattr(ncid, 'LatGridSpacing')
                    	nlon = cell.shape[1]
                    	nlat = cell.shape[0]
                    # THESE ARE FOR NEW CHL FILES
		    elif conf.v['file_type'] == 'chill':
                    	lat_array = ncid.variables['LAT'][:]
                    	lon_array = ncid.variables['LON'][:]
                    	latfirst, lonfirst = lat_array[0], lon_array[0]
                    	latdelta = np.average(lat_array[1:] - lat_array[:-1])
                    	londelta = np.average(lon_array[1:] - lon_array[:-1])
                    	nlat, nlon = lat_array.shape[0], lon_array.shape[0]

                    ncid.close()

                # HERE TAKING OUT THE RAW CODE AND SUBBING IN THE FUNCTION TO SEE IF IT WORKS
                    latinc = 1 # latitude is increasing with CHL files

                    nldn_feat = LT.cell_attribute(n_flashes['lat'], n_flashes['lon'], cell, latfirst, \
                            lonfirst, reg_conf.v['ctrLat'], reg_conf.v['ctrLon'], nlat, nlon, latdelta=latdelta, londelta=londelta, \
                            latincreasing=conf.v['latinc'], sdist=10)
                    n_flashes = rec.append_fields(n_flashes, 'feat', nldn_feat, usemask = False) # adding in the feature field
                    # still need to add in a seconds of day field


                    fin = np.where( (n_flashes['feat'] > 0) & (n_flashes['feat'] < 1e10) )
                    tj = fin[0].shape[0]
                    print 'NUM OF FLASHES ATTRIBUTED TO CELLS: ', tj

                    names = n_flashes.dtype.names 

    #		nldn_dict = [dict(zip(names, record)) for record in nldn_lines]
                    nldn_dict = {}
                    for n in names:
                        nldn_dict[n] = n_flashes[n]

                    # figure this part out, need to convert C and G to 1 and 0
                    cg_dict = {'C': 1, 'G': 0}
    #		for dummy in nldn_dict['cg']:
    #		    nldn_dict['cg'][dummy] = cg_dict[dummy] # this can probably get sped up with list comprension	
                    nldn_dict['cg'] = np.array( [cg_dict[c] for c in nldn_dict['cg']] ) 

                    new_nldn = LT.filter_dict(nldn_dict, fin)

                    ncwid = Dataset(nfile,'a')
                    if tj > 0:
                        ncwid.createDimension('n', new_nldn['lat'].shape[0])
                        # WRITING THE FLASH INFO
                        # MAKE FUNCTION THAT WILL LOOP OVER THE DICTIONARY??

                        LT.write_nc_dict(new_nldn, ncwid, 'n', 'd', prefix = 'nldn_flash_')

                    ncwid.close()

                else: ncid.close()


# SEND ME A TEXT SAYING CODE IS DONE
    #os.system("echo It is done | mail -s LMA flash processing 3202668731@vtext.com")


    


# BF ABOVE WORKS JUST FOR .DAT FILES BUT BELOW NEEDS TO WORK FOR NC FILES AS WELL
end_time = datetime.datetime.now()
duration = end_time - start_time
print 'END: %s' %end_time 
print 'DURATION: %s' %duration
#os.system('python plot_flash_case.py')
