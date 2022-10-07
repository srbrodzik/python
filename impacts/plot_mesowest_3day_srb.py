#!/usr/bin/python

'''
plot_wxdata_3day.py
Make 3-day rainfall plots (i.e. -2,-1,today)
Pull the data that you already collected in the other code
Run once per hour to update
No archiving, only real-time
Also needs to handle the case where the files don't exist
'''

from bs4 import BeautifulSoup
import os
import csv
import pdb
import time, datetime, glob
from time import gmtime, localtime, strftime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pytz
from dateutil import parser
tz = pytz.timezone('US/Pacific')
import math

#import Data
#utc = datetime.datetime.utcnow()
utc = datetime.datetime(2019,04,05,hour=0,minute=0,second=0,microsecond=1,tzinfo=None)
dt = datetime.timedelta(hours=1)
newtime = utc-dt
filetimestamp=newtime.strftime("%Y%m%d") #20141018
monthtimestamp=newtime.strftime("%Y%m") #201410
graphtimestamp=newtime.strftime("%m/%d/%y") #10/18/14
yeartimestamp = newtime.strftime("%Y") #2014
daytimestamp = newtime.strftime("%d") #18
hourtimestamp = newtime.strftime("%H") #hr
monthonlytimestamp = newtime.strftime("%m") #10
yyyymm=newtime.strftime("%Y%m") #201410

#Manually enter apu hostname: (added by Joe 9/26/14)
indir = '/home/disk/funnel/olympex/archive2/wxdata/Data/'+monthtimestamp+'/'+filetimestamp+'/'
sitelist = ['D1622','D7098','E0868','E4248','D0707','D3348','E5649','D4784','E4723','E4667','QTCW1',\
            'CWPF','CWSP','CWLM','CWYJ','CWKH','C5398','AT282','C2882','C5670','E2796','AT945','E4664',\
            'C5456','C8405','D5957','D7202','E6244','E6215','AS753','C6975','AS913','E4160',\
            'KUIL','KHQM','KSHN','KCLS','KOLM','KPWT','KCLM','KSEA','KFHR','KTIW','KNUW','KORS','KBLI',\
            'KBVS','KAWO','KBFI','KTCM','KPAE','KGRF',\
            'BKBW1','HUFW1','MIPW1','FRKW1','EMTW1','TOFW1','BKKW1','CGFW1','HWRW1','MFDW1','WFSW1','SKOW1','JEFW1','HSTW1',\
            'CDCW1','BCDW1','SKMW1','RSFW1','BRKW1','QCNW1','TT032','OLFW1','HKFW1','DOTW1','CYYJ','D4249','CMOW1','WKBW1','TILW1','ABNW1',\
            'RKCW1','CTKW1','CLSW1','D1695','EKMW1','CSHw1','HKOW1']
sitetitles = ['Copalis Beach','Aberdeen (Central)','Aberdeen (East)','Elma','Ocean Shores','Humptulips (US-101)','Brinnon','Hoodsport',\
              'Port Townsend','Sequim','Quinault CRN','Esquimalt Harbour','Sheringham','Victoria (city)',\
              'Victoria University','Malahat','Sequim (Blyn Springs Rd)','Sequim (Discovery Bay)','Anacortes',\
              'Coupeville','Camano Island','Nordland','Port Hadlock','Hansville','Kingston','Belfair','Grapeview',\
              'Poulsbo','Sequim (Dungeness)','Sequim (Bear Creek)','Rochester','Tenino','Shelton (Arcadia)',\
              'KUIL Quillayute State Airport', 'KHQM Hoquiam Bowerman Airport','Shelton Sanderson Field',\
              'Chehalis-Centralia Airport','Olympia Airport','Bremerton National Airport',\
              'Port Angeles (Fairchild Airport)','Seatac Airport','Friday Harbor Airport',\
              'Tacoma Narrows Airport','Whidbey Island Naval Air Station','Eastsound, Orcas Island Airport',\
              'Bellingham International Airport','Burlington/Mt. Vernon (Skagit Arpt)','Arlington Municipal Airport',\
              'Seattle Boeing Field','Tacoma McChord AFB','Everett (Paine Field)','Tacoma Ft Lewis',\
              'Black Knob','Humptulips (2400 ft)','Minot Peak (1768 ft)','Forks','Ellis Mountain (2305 ft)','Toms Creek (2400 ft)',\
              'Buck Knoll (1630 ft)','Cougar Mountain (2400 ft)','Haywire Ridge','Matlock','West Fork Satsop River',\
              'North Fork Skokomish River','Jefferson Creek','Hoh Rain Forest','Cedar Creek','Skookumchuck River (lower)',\
              'Skookumchuck River (upper)','Riverside Fire Authority','Brooklyn','Quilcene','Paradise Fire (5400 ft)','Owl Mountain',\
              'Huckleberry Ridge','Doty','Victoria (airport)','Sequim (Dungeness Bay)','Clam Mountain','Windy Knob','Tilton River',\
              'Abernathy Mountain','Rock Creek','Thrash Creek','Chehalis','Ocean Park','Elk Meadows','Cushman Dam','Hoko']

def read_csv_file(infile):
    global prec_accum
    with open(infile, 'r') as f:
        prec_tmp = []
        for line in f:
            data=line.split(',')
            year = int(data[0])
            month = int(data[1])
            day = int(data[2])
            hour = int(data[3])
            minute = int(data[4])
            pytime.append(datetime.datetime(year,month,day,hour,minute))
            temp.append(float(data[5]))
            wind_spd.append(float(data[7]))
            wind_gust.append(float(data[8]))
            wind_dir.append(float(data[9]))
            prec_tmp.append(float(data[10]))
            altimeter.append(float(data[11]))
            mslp.append(float(data[12]))
    if len(prec_accum) == 0: #have to add prev total if not the first file
        prec_accum += prec_tmp
    else:
        add_data = [x+prec_accum[-1] for x in prec_tmp]
        prec_accum = prec_accum + add_data
    return pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp

def add_dummy_value():
    #Adds -999.99 value if gap in obs of 2 hours or more
    #Makes plots show a gap in missing data rather than a line
    #Does not take numpy arrays, only lists
    gap = datetime.timedelta(0, 7200)
    for i,atime in enumerate(pytime):
        try:
            dt = pytime[i+1]-pytime[i]
        except:
            break
        if dt > gap:
            pytime.insert(i+1,pytime[i+1]-gap)
            temp.insert(i+1,-999.99)
            wind_spd.insert(i+1,-999.99)
            wind_gust.insert(i+1,-999.99)
            wind_dir.insert(i+1,-999.99)
            prec_accum.insert(i+1,-999.99)
            altimeter.insert(i+1,-999.99)
            mslp.insert(i+1,-999.99)


def load_data_3days(site):
    factarr = [49,25,1]
    global pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp
    pytime = []
    temp = []
    wind_spd = []
    wind_gust = []
    wind_dir = []
    prec_accum = []
    altimeter = []
    mslp = []
    for fact in factarr:
        dt = datetime.timedelta(hours=fact)
        newtime2 = utc-dt
        monthtimestamp =newtime2.strftime("%Y%m")
        filetimestamp=newtime2.strftime("%Y%m%d")
        indir = '/home/disk/funnel/olympex/archive2/wxdata/Data/'+monthtimestamp+'/'+filetimestamp+'/'
        if site[0] == 'K':
            try:
                infile = indir+'wxdata_NWS_'+site+'_'+filetimestamp+'.csv' 
                pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp = read_csv_file(infile)
            except:
                print 'File not found: '+infile
        else:
            try:
                infile = indir+'wxdata_'+site+'_'+filetimestamp+'.csv' 
                pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp = read_csv_file(infile)
            except:
                print 'File not found: '+infile

    return pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp
     
def make_3day_plot(site,sitetitle,pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp):
    #Get data using other functions
    
    #Figure out which panels to plot...
    flagarr = [1,1,1,1,1]
    #if precip doesn't exist (because it hasn't rained), make an array of all 0's:
    if max(prec_accum) < 0:
        flagarr[2] = 0
        #precip_1hr = ['0']*len(pytime_utc)
    if max(wind_spd) < 0:
        flagarr[0] = 0
        flagarr[1] = 0
    if max(mslp) < 0 and max(altimeter) < 0:
        flagarr[3] = 0
    elif max(mslp) < 0 and max(altimeter) > 0:
        flagarr[3] = 2 #use altimeter if MSLP not available
    if max(temp) < -100 and max(mslp) < 0 and max(wind_spd) < 0:
        flagarr[4] = 0 #river gauges that have rain and nothing else

    #---------Plot----------
    #5-panel plot (temp, wind, wind dir, precip accum, snow depth)
    #Precip accum will eventually be moved to another plot w/ Parsivel

    #markersize
    if len(prec_accum) > 250: 
        markersize = 0.8
    elif len(prec_accum) <= 250 and len(prec_accum) > 150: 
        markersize = 2.0
    elif len(prec_accum) <= 150 and len(prec_accum) > 75: 
        markersize = 3.5
    else:
        markersize = 4

    import matplotlib.pyplot as plt
    from matplotlib.dates import DayLocator, HourLocator, DateFormatter

    pytime = np.array(pytime)
    temp = np.array(temp)
    wind_spd = np.array(wind_spd)
    wind_gust = np.array(wind_gust)
    wind_dir = np.array(wind_dir)
    prec_accum = np.array(prec_accum)
    altimeter = np.array(altimeter)
    mslp = np.array(mslp)

    apu_str = sitetitle

    if flagarr[4] == 1:
        fig = plt.figure()
        ax1 = fig.add_subplot(5,1,1)
        fig.set_size_inches(8,8)
        graphtimestamp_start=pytime[0].strftime("%m/%d/%y") #10/18/14
        ax1.set_title(apu_str+' '+graphtimestamp_start+' - '+graphtimestamp)

        goodtemp = np.where((temp > -100))[0]
        badtemp = np.where((temp < -100))[0]
        if len(badtemp) > 0:
            temp[badtemp] = float('nan')

        ax1.plot_date(pytime,temp, 'o-',label="Temp ($^\circ$C)",linewidth=1.2, color="black",markersize=markersize)

        try:
            ymin = int(math.floor(min(temp[goodtemp]) / 10.0)) * 10 - 0.1
            ymax = int(math.ceil(max(temp[goodtemp]) / 10.0)) * 10 + 0.1
        except:
            ymin = 0
            ymax = 10
        
        if ymin == ymax:
            ymin -= 10 #in case both round to same thing
        ax1.set_ylim(ymin,ymax) #round sets range to nearest 10
        ax1.set_ylabel('Temp ($^\circ$C)')
    
        #darker line at T= 0C (if applicable)
        ax1.plot([pytime[0],pytime[-1]], [0,0], '-', lw=0.5, color="black", alpha=1)
        axes = [ax1]

    #----------Ax2 wind---------------
    if flagarr[0] == 1: #only run if wind data exists...combine speed and direction
	ax2 = plt.subplot(5,1,2)

        
	badwind = np.where((wind_dir < 0))[0]
	if len(badwind) > 0:
	    wind_dir[badwind] = float('nan')
	badspd = np.where((wind_spd < 0))[0]
	if len(badspd) > 0:
	    wind_spd[badspd] = float('nan')
	badgust = np.where((wind_gust < 0))[0]
	if len(badgust) > 0:
	    wind_gust[badgust] = float('nan')
	badgust2 = np.where((wind_gust == 0))[0]
	if len(badgust2) > 0: #no direction if 0 wind
	    wind_gust[badgust2] = float('nan')

	#because of some annoying requirement, the first element of the array can't be nan, so have to remove up to first NAN
	goodmax = np.where((wind_gust > 0))[0]
	if len(goodmax)>0:
	    pytime_wnd = pytime[goodmax[0]:]
	    wind_max_good = wind_gust[goodmax[0]:]
	else:
	    pytime_wnd = pytime
	    wind_max_good = wind_gust
	goodavg = np.where((wind_spd >= 0))[0]
	if len(goodavg)>0:
	    pytime_avg = pytime[goodavg[0]:]
	    wind_spd_good = wind_spd[goodavg[0]:]
	else:
	    pytime_avg = pytime
	    wind_spd_good = wind_spd

	ax2.plot_date(pytime_avg,wind_spd_good,'o-',label='Speed', color="blue",markersize=markersize)
	ax2.plot_date(pytime_wnd,wind_max_good,'o-',label='Gust', color="red",markersize=markersize)

	# add some text for labels, title and axes ticks
	ax2.set_ylabel('Wind (mph)')
	if len(goodmax) > 0:
	    ax2.set_ylim(-0.1*max(wind_max_good),int(max(wind_max_good)+max(wind_max_good)*0.2))
	else:
	    ax2.set_ylim(-0.1*max(wind_spd_good),int(max(wind_spd_good)+max(wind_spd_good)*0.2))
	ax2.legend(prop={'size':9},loc="upper left")


	#---------Ax3 wind direction---------
	ax3 = plt.subplot(5,1,3)

	ax3.plot_date(pytime,wind_dir,'o-',markersize=markersize,color="orange",linewidth=1.2, label="Direct")

	# add some text for labels, title and axes ticks
	ax3.set_ylabel('Wind Direction')
	ax3.set_ylim(-10,370)

	#Set major tick locs
	ax3.set_yticks([0,90,180,270,360])
        ax3.get_yaxis().set_label_coords(-0.075,0.5) #Force wind direction label to align with others...not sure why it acts weird
	#ax3.set_yticklabeks(['N','E','S','W','N'])

	axes.append(ax2)
        axes.append(ax3)
    #---------Ax4 precip---------
    if flagarr[2] == 1 and flagarr[0] == 1:
	ax4 = plt.subplot(5,1,4)
    else: #if no wind
	ax4 = plt.subplot(5,1,2)
    if flagarr[4] == 0: #occurs if this is the first and only plot
        fig = plt.figure()
        ax4 = fig.add_subplot(5,1,1)
        fig.set_size_inches(8,8)
        graphtimestamp_start=pytime[0].strftime("%m/%d/%y") #10/18/14
        ax4.set_title(apu_str+' '+graphtimestamp_start+' - '+graphtimestamp)
    if flagarr[2] == 1:
	labelname = 'Precip ('+"{:.1f}".format(prec_accum[np.where((prec_accum >= 0))[0]][-1])+' mm)'
	ax4.plot_date(pytime,prec_accum,'o-',label=labelname,markersize=markersize,color="green",linewidth=1.2)

	# add some text for labels, title and axes ticks
	ax4.set_ylabel('Precip (mm)')

	if max(prec_accum) > 0:
	    ax4.set_ylim(-0.1*max(prec_accum),max(prec_accum)+max(prec_accum)*0.2)
	else:
	    ax4.set_ylim(-0.5,5)
	ax4.legend(prop={'size':9},loc="upper left")
        if flagarr[4] == 1:
	    axes.append(ax4)
        if flagarr[4] == 0:
            axes = [ax4]
    #---------Ax5 MSLP-------- (only if it exists)
    if min(flagarr) >= 1: #all plots exist
	ax5 = plt.subplot(5,1,5) 
    elif flagarr[3] >= 1 and flagarr[0] == 0 and flagarr[2] == 0: #no wind or precip
	ax5 = plt.subplot(5,1,2) 
    elif flagarr[3] >= 1 and flagarr[0] == 1 and flagarr[2] == 0: #wind, no precip
	ax5 = plt.subplot(5,1,4)         
    elif flagarr[3] >= 1 and flagarr[0] == 0 and flagarr[2] == 1: #no wind, precip
	ax5 = plt.subplot(5,1,3)          

    if flagarr[3] == 1:
        axes.append(ax5)

	#remove bad data (usually = 234 when bad)
	goodmslp = np.where((mslp > 0))[0]
	if len(goodmslp) > 0:
            badmslp = np.where((mslp < -100))[0]
            if len(badmslp) > 0:
                mslp[badmslp] = float('nan')
	    ax5.plot_date(pytime,mslp,'o-',label='MSLP',color="red",linewidth=1.2)
	else:
	    mslp = [0]

	# add some text for labels, title and axes ticks
	ax5.set_ylabel('MSLP (hPa)')

	if max(mslp[goodmslp]) > 0:
	    ax5.set_ylim(int(min(mslp[goodmslp])-3),int(max(mslp[goodmslp])+3))
	else:
	    ax5.set_ylim(1000,1020)
        ax5.get_yaxis().get_major_formatter().set_useOffset(False) #kill exponents

    if flagarr[3] == 2: #use altimeter instead of MSLP
        axes.append(ax5)

	#remove bad data (usually = 234 when bad)
	goodalt = np.where((altimeter > 0))[0]
	if len(goodalt) > 0:
            badalt = np.where((altimeter < -100))[0]
            if len(badalt) > 0:
                altimeter[badalt] = float('nan')
	    ax5.plot_date(pytime,altimeter,'o-',label='Altimeter',markersize=markersize,color="red",linewidth=1.2)
	else:
	    altimeter = [0]

	# add some text for labels, title and axes ticks
	ax5.set_ylabel('Altimeter (hPa)')

	if max(altimeter[goodalt]) > 0:
	    ax5.set_ylim(int(min(altimeter[goodalt])-3),int(max(altimeter[goodalt])+3))
	else:
	    ax5.set_ylim(1000,1020)
        ax5.get_yaxis().get_major_formatter().set_useOffset(False) #kill exponents

    #last axis gets hour label
    axes[-1].set_xlabel('Time (UTC)')

    #Make sure x-lims match (in case of missing data)
    #delta = datetime.timedelta(minutes=10)
    min_pytime = min([pytime[0]])
    max_pytime = max([pytime[-1]])
    #adjust beginning and end times
    min_pytime = min_pytime.replace(hour=0,minute=0)
    max_pytime = max_pytime + datetime.timedelta(minutes=15)
    min_pytime = min_pytime - datetime.timedelta(minutes=15)
   

    #Tick and label stuff that repeats for every plot
    for ax in axes:
	ax.spines["top"].set_visible(False)  
	ax.spines["right"].set_visible(False)  
	ax.spines["left"].set_visible(False)
	ax.spines["bottom"].set_visible(False)

	ax.tick_params(axis='x',which='both',bottom='off',top='off')   
	ax.tick_params(axis='y',which='both',left='off',right='off')     

        #major dates (bottom axis only)
	ax.xaxis.set_major_locator( DayLocator() )
	ax.xaxis.set_major_formatter( DateFormatter('%b-%d') )
	ax.xaxis.set_minor_locator( HourLocator(np.linspace(6,18,3)) )
	ax.xaxis.set_minor_formatter( DateFormatter('%H') )

	#ax.xaxis.set_major_locator( HourLocator(np.linspace(0,24,7)))
	#ax.xaxis.set_major_formatter( DateFormatter('%H') )	

        ax.fmt_xdata = DateFormatter('Y%m%d%H%M%S')

	#Tick size
	ax.tick_params(axis = 'both', which = 'major', labelsize = 9)
	ax.tick_params(axis = 'x', which = 'minor', labelsize = 9)

        #give proper coordinates to dashed lines
        pytime_tics = np.append(min_pytime,pytime)
        pytime_tics = np.append(pytime_tics,max_pytime)

	ticlocs = ax.get_yaxis().get_majorticklocs()
	for y in ticlocs:
	    ax.plot(pytime_tics, [y] * len(range(0,len(pytime_tics))), '--', lw=0.5, color="black", alpha=0.3)

	ax.get_yaxis().set_label_coords(-0.07,0.5)
	ax.set_xlim(min_pytime,max_pytime)
        #pdb.set_trace()

    #fig.set_title(apu_str+' '+graphtimestamp)

    #saving files
    outdir_realtime = '/home/disk/user_www/brodzik/olympex/wxdata/real_time/'
    outdir_archive = '/home/disk/user_www/brodzik/olympex/wxdata/archive/'+site+'/'
    if not os.path.exists(outdir_archive):
        os.mkdir(outdir_archive)
    os.chdir(outdir_archive)
    os.system('mkdir -p ' + yyyymm)
    plt.savefig(outdir_archive+yyyymm+'/'+site+'_'+filetimestamp+'_weather_3day.png', bbox_inches='tight')
    plt.savefig(outdir_realtime+site+'_weather_3day_realtime.png', bbox_inches='tight')

    #Stacy's files
    #if site == 'KPWT' or site == 'KHQM' or site == 'KCLM' or site == 'KUIL' or site == 'KSHN' or site == 'CYYJ' or site == 'CWLM':
    #    os.chdir("/home/disk/funnel/olympex/archive/ops/asos/")
    #    os.system('mkdir -p ' + filetimestamp)
    #    stacydir = '/home/disk/funnel/olympex/archive/ops/asos/'+filetimestamp+'/'
    #    site2 = site.lower()
    #    plt.savefig(stacydir+'ops.asos.'+filetimestamp+'0000.'+site2+'.png', bbox_inches='tight') 

    plt.close('all')
    print 'Plotted:' +site


for i,site in enumerate(sitelist):
    sitetitle = sitetitles[i]
    pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp = load_data_3days(site)
    if len(pytime) > 0:
        add_dummy_value()
        make_3day_plot(site,sitetitle,pytime,temp,wind_spd,wind_gust,wind_dir,prec_accum,altimeter,mslp)
   




