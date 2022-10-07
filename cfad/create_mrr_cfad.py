'''
create_mrr_cfad.py
Joe Zagrodnik 19-July-2016
Create an MRR CFAD
Using my .nc files from the "Ave" data
Start by trying for one day, later expand to BB detection/other days

Known bugs: (28-July)
-Not happy with some areas flagged as non-BB
-Returns BB height even if non-BB
'''

import numpy as np
import netCDF4 as nc
import pdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import calendar
from matplotlib.colors import LogNorm


import pyart
from matplotlib import cm, colors

#Colormaps
cmap3 = []
for value, color in zip([0,2./10,3./10,4./10,5./10,6./10,7./10,8./10,0.9,1],\
    ['#990000','#ffffff','#aec2ea','#2951a3','#152951','#732673','#c03fc0','#ff6600','#ffb380','#000066']):
    cmap3.append((value,color))
    velmap = colors.LinearSegmentedColormap.from_list("custom",cmap3)
cmap4 = []
for value, color in zip([0,1./10,2./10,3./10,4./10,5./10,6./10,7./10,8./10,9./10,1],\
    ['#ccd8ff','#3366ff','#9fdf9f','#00b300','#ffff00','#ffcc30','#e62e00','#ff6600','#fff0e5','#c03fc0','#602060']):
    cmap4.append((value,color))
dbzmap = colors.LinearSegmentedColormap.from_list("custom",cmap4)
carbone = matplotlib.cm.get_cmap('pyart_Carbone42')

#stuff for converting the timestamps
def unix2timestamp(unix):
  return datetime.datetime.utcfromtimestamp(unix).strftime("%Y%m%d%H%M")
def get_title_timestamp(unix):
  return datetime.datetime.utcfromtimestamp(unix).strftime("%Y-%b-%d")
def timestamp2unix(timestamp):
  return calendar.timegm(datetime.datetime(year = int(timestamp[0:4]), month = int(timestamp[4:6]), day = int(timestamp[6:8]), hour = 0, minute = 0, second = 0).timetuple())
def convert_to_datetime(timearr):
    dt_arr = []
    for t in timearr:
        ttmp = unix2timestamp(t)
        dt_arr.append(datetime.datetime(int(ttmp[0:4]),int(ttmp[4:6]),int(ttmp[6:8]),int(ttmp[8:10]),int(ttmp[10:12])))
    dt_arr = np.array(dt_arr)
    return dt_arr

#-------------BB detection algorithm-----------------
def eval_dbz10(dbz):
    #check if near-surface (3rd bin from surface) is greater than 10 dbz
    if dbz[2] <= 10 or np.isnan(dbz[2]) == True: return 0
    else: return 1

def eval_velocity(vel,height):
    #change of < -1.5 m/s over two range gates (200 m)
    velheights = height[1:-1] #middle height of differencing
    vel_diff = vel[2:]-vel[0:-2] #difference values
    if np.min(vel_diff) <= -1.5: velcheck = 1
    else: velcheck = 0
    if velcheck == 1:
        vel_good = vel_diff[np.where((vel_diff <= -1.5))[0]]
        gate_good = np.where((vel_diff <= -1.5))[0]
    else:
        vel_good = np.array([])
        gate_good = np.array([])    
    return velcheck, vel_good, gate_good

def eval_dbz(dbz,height):
     #change of > 2.5 dBZ over two range gates (200 m)
    dbzheights = height[1:-1] #middle height of differencing
    dbz_diff = dbz[2:]-dbz[0:-2] #difference values
    if np.max(dbz_diff) > 2.5: dbzcheck = 1
    else: dbzcheck = 0
    if dbzcheck == 1:
        dbz_good = dbz_diff[np.where((dbz_diff >= 1.0) & (dbz[1:-1] > 10))[0]]
        gate_good = np.where((dbz_diff >= 1.0) & (dbz[1:-1] > 10))[0]
    else:
        dbz_good = np.array([])
        gate_good = np.array([])   
    return dbzcheck, dbz_good, gate_good

def average_dbz(dbzarr):
    #Average an array of dbz values (top, bottom of BB)
    ztmp = 10**(dbzarr/10)
    zmean = np.mean(ztmp)
    dbz_avg = 10*np.log(zmean)
    return dbz_avg

def eval_hybrid_profile(dbz,height,bb_loc_bot):
    #check dbz gradient from BB bottom to surface
    height_diff_tmp2 = height[bb_loc_bot] - height[2] #bb_bot to 3rd bin height
    dbz_diff = dbz[bb_loc_bot] - dbz[2]
    dbz_gradient = dbz_diff/(height_diff_tmp2/1000)
    if dbz_gradient <= -0.1: hybridcheck = 1
    else: hybridcheck = 0
    return hybridcheck, dbz_gradient, height_diff_tmp2, dbz[2], dbz[bb_loc_bot]

def find_bb_height(dbz,height,bb_loc_bot):
    #check 500 m above BB bottom for main bright band peak
    bb_possibilities = dbz[bb_loc_bot-1:bb_loc_bot+6]
    bb_height = height[np.where((dbz > np.max(bb_possibilities - 0.001)) & (dbz < np.max(bb_possibilities + .001)))[0]]
    bb_dbz = dbz[np.where((dbz > np.max(bb_possibilities - 0.001)) & (dbz < np.max(bb_possibilities + .001)))[0]]
    return bb_height[0],bb_dbz[0]  

class mrr_raw(object):
    '''
    Class to hold MRR raw data
    '''
    def __init__(self):
        self.time = []
        self.rangegate = []
        self.Z = []
        self.rainrate = []
        self.lwc = []
        self.vel = []

    def load_mrr(self,infile):
        data = nc.Dataset(infile,'r')          
        self.time = data.variables['time'][:] #1440
        self.rangegate = data.variables['MRR rangegate'][:] #1440 by 31
        self.Z = data.variables['MRR_Capital_Z'][:] #1440 by 31
        self.rainrate = data.variables['MRR_RR'][:] #1440 by 31
        self.lwc = data.variables['MRR_LWC'][:] #1440 by 31
        self.vel = data.variables['MRR_W'][:] #1440 by 31

class mrr_cfad(object):
    '''
    Class to hold MRR CFAD data
    '''
    def __init__(self):
        self.time = []
        self.Z = []
        self.vel = []
        self.bb_height = []
        self.category = [] #category (Non-bb, bb, hybrid)

    def add_first_time(self,times,ranges,Zs,vels,bb_height,category):
        self.time = times
        self.rangegate = ranges
        self.Z = Zs
        self.vel = vels
        self.bb_height = bb_height
        self.category = category

    def add_additional_times(self,times,ranges,Zs,vels,bb_height,category):
        self.time = np.hstack((self.time,times))
        self.rangegate = np.vstack((self.rangegate,ranges))
        self.Z = np.vstack((self.Z,Zs))
        self.vel = np.vstack((self.vel,vels))
        self.bb_height = np.hstack((self.bb_height,bb_height))
        self.category = np.hstack((self.category,category))

    def pickle_cfad(self,outfile):
        import cPickle as pickle
        struct = (self.time,self.rangegate,self.Z,self.vel,self.bb_height,self.category)
        pickle.dump( struct, open( outfile, 'wb' ) )        

        
def brightband_algorithm(mrr):
    '''
    Brightband detection algorithm
    '''
    dt_time = convert_to_datetime(mrr.time)

    sttime = [] #start time of 30-min period
    num_prof_gt_10dbz = [] #num profiles with 3rd bin > 10 dBZ (need > 15)
    prof_type = [] #1 for nonBB, 2 for BB, 3 for hybrid
    bb_bottom = [] #-999 if no BB, otherwise range height where BB is detected
    slope_below_bb = [] #slope of profile below BB

    add_30m = datetime.timedelta(minutes=30) 
    starttime = dt_time[0]
    endtime = starttime + datetime.timedelta(days=1) #00:00 UTC next day
    starttime = dt_time[0].replace(hour=0,minute=0) #00:00 UTC current day

    while starttime < endtime:
        #1. Pull out 30-min time period
        goodt = np.where((dt_time >= starttime) & (dt_time < starttime+add_30m))[0]
        starttime += add_30m
        #define check lists
        dbz10_check = []
        vel_check = []
        dbz_check = []
        hybrid_check_ind_tmp = [] #hybrid check using >50% of individual profiles
        height_diff = []
        dbz_300m = []
        dbz_bottomofbb = []
        bb_height = []
        bb_dbz = []
        for t in goodt:
            dbz10_check.append(eval_dbz10(mrr.Z[t,:]))
            vel_check_tmp, vel_good, vel_gate_good = eval_velocity(mrr.vel[t,:],mrr.rangegate[t,:])
            vel_check.append(vel_check_tmp)
            dbz_check_tmp, dbz_good, dbz_gate_good = eval_dbz(mrr.Z[t,:],mrr.rangegate[t,:])
            dbz_check.append(dbz_check_tmp)
            vel_dbz_match = np.intersect1d(vel_gate_good,dbz_gate_good) #checks that vel and dbz gradient are at same level
            vel_dbz_match = vel_dbz_match[np.where((vel_dbz_match > 1))[0]] #removes matches in the lowest 2 levels
            if len(vel_good) > 0 and len(dbz_good) > 0 and len(vel_dbz_match) > 0 and dbz10_check[-1] == 1:
                hybrid_check_tmp, dbz_gradient, height_diff_tmp, dbz_3rdbin, dbz_bbbot = eval_hybrid_profile(mrr.Z[t,:],mrr.rangegate[t,:],min(vel_dbz_match))
                hybrid_check_ind_tmp.append(hybrid_check_tmp)
                height_diff.append(height_diff_tmp) #bottom of bb height diff
                dbz_300m.append(dbz_3rdbin)
                dbz_bottomofbb.append(dbz_bbbot)
                bb_height_tmp,bb_dbz_tmp = find_bb_height(mrr.Z[t,:],mrr.rangegate[t,:],min(vel_dbz_match))
                bb_height.append(bb_height_tmp)
                bb_dbz.append(bb_dbz_tmp)
        #Now go through 30-min period and figure out precip type
        dbz10_check = np.array(dbz10_check)
        vel_check = np.array(vel_check)
        dbz_check = np.array(dbz_check)
        hybrid_check = np.array(hybrid_check_ind_tmp)

        good_10dbz = len(np.where((dbz10_check == 1))[0])
        good_bb = np.where((dbz10_check == 1) & (vel_check == 1) & (dbz_check == 1))[0]
        good_save_cfad = goodt[np.where((dbz10_check == 1))[0]]  #good profiles to save for CFAD (meets 10dBZ criteria)

        #fill in final categories
        try:
            final_time.append(dt_time[goodt][0]) #start of time period
            final_num_good_10dbz_profiles.append(len(np.where((dbz10_check == 1))[0]))

            if good_10dbz >= 15: #don't evaluate if less than half of profiles are at least 10 dbz:
                if (len(np.where((vel_check == 1) & (dbz_check == 1))[0])/float(len(vel_check)) < 0.3) or (len(vel_dbz_match) == 0): #at least 30% of profiles have BB
                    final_category.append(0) #NON-BB
                else: #check for hybrid
                    hcit = np.array(hybrid_check_ind_tmp)
                    if len(np.where((hcit == 1))[0])/float(len(hcit)) >= 0.5:
                        final_category.append(2) #HYBRID
                    else:
                        final_category.append(1) #BB
                final_num_good_vel_profiles.append(len(np.where((vel_check == 1))[0]))
                final_num_good_dbz_profiles.append(len(np.where((dbz_check == 1))[0]))
                final_num_good_bb_profiles.append(len(good_bb))
                final_num_good_hybrid_profiles.append(len(np.where((hybrid_check == 1))[0]))
                if len(height_diff) > 0:
                    final_bb_height.append(np.mean(bb_height))
                    final_bb_dbz.append(np.mean(bb_dbz))
                else:
                    final_bb_height.append(-1)
                    final_bb_dbz.append(-1)
                
                #Add profiles to CFAD... 
                if (final_time[-1] >= ar_start[0] and final_time[-1] < ar_end[0]) or (final_time[-1] >= ar_start[1] and final_time[-1] < ar_end[1])\
                or (final_time[-1] >= ar_start[2] and final_time[-1] < ar_end[2]):
                #non-AR
                #if (final_time[-1] <= ar_start[0]) or (final_time[-1] >= ar_end[0] and final_time[-1] < ar_start[1])\
                #or (final_time[-1] >= ar_end[1] and final_time[-1] < ar_start[2]) or (final_time[-1] >= ar_end[2]):
                    bb_cfad =  np.zeros(len(dt_time[good_save_cfad])) #holder for mean bb height
                    bb_category =  np.zeros(len(dt_time[good_save_cfad])) #holder for mean bb height
                    bb_cfad[:] = final_bb_height[-1]
                    bb_category[:] = final_category[-1]
                    if len(cfad.time) == 0:
                        cfad.add_first_time(dt_time[good_save_cfad],mrr.rangegate[good_save_cfad,:],mrr.Z[good_save_cfad,:],mrr.vel[good_save_cfad,:],bb_cfad,bb_category)
                    else:
                        cfad.add_additional_times(dt_time[good_save_cfad],mrr.rangegate[good_save_cfad,:],mrr.Z[good_save_cfad,:],mrr.vel[good_save_cfad,:],bb_cfad,bb_category)
                    print 'Added to CFAD: '+final_time[-1].strftime("%Y-%b-%d %H:%M")+' '+str(len(good_save_cfad))+' elements'

            else:
                final_category.append(-1) #fail
                final_num_good_vel_profiles.append(-1) 
                final_num_good_dbz_profiles.append(-1) 
                final_num_good_bb_profiles.append(-1) 
                final_num_good_hybrid_profiles.append(-1) 
                final_bb_height.append(-1) 
                final_bb_dbz.append(-1)

            print 'Time: '+final_time[-1].strftime("%Y-%b-%d %H:%M")
            #print 'BB height: '+str(final_bb_height[-1])
            #print 'Category '+str(final_category[-1])
        except:
            print 'No profiles found for: '+starttime.strftime('%Y-%m-%d %H:%M')
 
    fg10dbz = final_num_good_10dbz_profiles
    return final_time,final_category,final_num_good_vel_profiles,final_num_good_dbz_profiles,final_num_good_bb_profiles,final_num_good_hybrid_profiles,final_bb_height,dt_time,fg10dbz
    

#---------------plot----------------------------------
def plot_mrr(mrr,final_time,final_category,final_num_good_vel_profiles,final_num_good_dbz_profiles,final_num_good_bb_profiles,final_num_good_hybrid_profiles,final_bb_height,dt_time,fg10dbz):
    graphtimestamp = get_title_timestamp(mrr.time[0])
    date = unix2timestamp(mrr.time[0])
    starttime = dt_time[0]
    fig=plt.figure(figsize=(10, 13))
  
    ax1 = fig.add_subplot(511)
    ax1.set_title(label+' MRR '+graphtimestamp)
    levels = np.arange(0,50,0.1)
    plotCF = ax1.contourf(dt_time,mrr.rangegate[-1,:], np.transpose(mrr.Z), levels,cmap=carbone, extend="both")#plt.get_cmap("spectral")
    ax1.plot(final_time,final_bb_height,'kx',label='BB height')
    cbZ=plt.colorbar(plotCF)
    cbZ.set_label('Reflectivity [dBZ]')

    ax2 = fig.add_subplot(512)
    levels = np.arange(-2,10,0.1)
    plotCF = ax2.contourf(dt_time,mrr.rangegate[-1,:], np.transpose(mrr.vel), levels,cmap=velmap, extend="both")#plt.get_cmap("spectral")
    cbZ=plt.colorbar(plotCF)
    cbZ.set_label('Fall Velocity [m s$^{-1}$]')

    ax3 = fig.add_subplot(513)
    ax3.plot(final_time,fg10dbz,'kx',label='10 dBZ')
    ax3.plot(final_time,final_num_good_vel_profiles,'rx',label='Vel')
    ax3.plot(final_time,final_num_good_dbz_profiles,'gx',label='dBZ')
    ax3.plot(final_time,final_num_good_hybrid_profiles,'bx',label='Hyb')

    ax4 = fig.add_subplot(514)
    ax4.plot(final_time,final_category,'kx',label='Final designation')

    ax5 = fig.add_subplot(515)
    ax5.plot(final_time,final_bb_height,'kx',label='BB height')

    axes = [ax1,ax2,ax3,ax4,ax5]
    axes[-1].set_xlabel('Time (UTC)')
    for ax in axes:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
        ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,2)))
        ax.set_xlim(starttime, starttime+datetime.timedelta(days=1))

        if ax == ax1 or ax == ax2 or ax == ax5:
            ax.set_ylim(0,3000)
            ax.set_ylabel('Height above radar (m)')
        else:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
  
            if ax == ax3:
                ax.set_ylim(-1,31)
                ax.set_ylabel('Number of 1-min profiles')
                #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            else:
                ax.set_ylim(-1.5,2.5)
                ax.set_ylabel('Category')

    plt.tight_layout()
    fig.savefig(outdir+loc+'_MRR_timeseries_'+date+'.png',bbox_inches='tight')
    plt.close()
    print 'Plotted: '+graphtimestamp


startdate = datetime.datetime(2015,11,7)
enddate = datetime.datetime(2015,12,10)
dt_date = datetime.timedelta(days=1)
loc = 'bishop'
label = 'Bishop'

#Full AR
ar_start = [datetime.datetime(2015,11,12,12),datetime.datetime(2015,11,16,12),datetime.datetime(2015,12,8,2)]
ar_end = [datetime.datetime(2015,11,14,0),datetime.datetime(2015,11,18,1),datetime.datetime(2015,12,9,12)]

#warm sectors
#ar_start = [datetime.datetime(2015,11,13,0),datetime.datetime(2015,11,17,1),datetime.datetime(2015,12,8,11)]
#ar_end = [datetime.datetime(2015,11,14,0),datetime.datetime(2015,11,18,11),datetime.datetime(2015,12,9,9)]

#Prefrontal
#ar_start = [datetime.datetime(2015,11,12,12),datetime.datetime(2015,11,16,12),datetime.datetime(2015,12,8,2)]
#ar_end = [datetime.datetime(2015,11,13,0),datetime.datetime(2015,11,17,1),datetime.datetime(2015,12,8,11)]

cfad = mrr_cfad()

#final variables
final_time = []
final_category = [] #0: non-bb, 1: bb, 2: hybrid, -1 = not enough data (15 profiles)
final_num_good_10dbz_profiles = []
final_num_good_vel_profiles = []
final_num_good_dbz_profiles = []
final_num_good_bb_profiles = []
final_num_good_hybrid_profiles = []
final_bb_height = [] #height of BB
final_bb_dbz = [] #dbz at final_bb_height



while startdate < enddate:
    
    yyyymm = startdate.strftime("%Y%m")
    mmdd = startdate.strftime("%m%d")
    mrrdir = '/home/disk/bob/olympex/ground/MRR/'+loc+'/NetCDF/Ave/'+yyyymm+'/'
    #Timeseries plot test (to make sure data is ok)
    outdir = '/home/disk/meso-home/jzagrod/Olympex/MRR/CFAD/Test/Plot/Timeseries/bishop/'
    #outdir = '/home/disk/user_www/jzagrod/olympex/comparison/'

    infile = mrrdir+mmdd+'.nc'

    mrr = mrr_raw()
    mrr.load_mrr(infile)

    final_time,final_category,final_num_good_vel_profiles,final_num_good_dbz_profiles,final_num_good_bb_profiles,final_num_good_hybrid_profiles,final_bb_height,dt_time,fg10dbz = brightband_algorithm(mrr)
    print len(final_time)
    plot_mrr(mrr,final_time,final_category,final_num_good_vel_profiles,final_num_good_dbz_profiles,final_num_good_bb_profiles,final_num_good_hybrid_profiles,final_bb_height,dt_time,fg10dbz)

    startdate+= dt_date

outfile = '/home/disk/meso-home/jzagrod/Olympex/MRR/CFAD/Test/Data/fishery_atmos_rivers.p'
cfad.pickle_cfad(outfile)
pdb.set_trace()
