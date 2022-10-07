# -*- coding: utf-8 -*-
"""
Dec 2016
Joe Zagrodnik
Load NPOL RHI
Load DOW RHI nearest to NPOL
Plot 52 degree NPOL RHI
Plot 56 degree DOW RHI on top of NPOL
Draw lines separating the two
Make appropriate adjustments to Title

"""

import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyart
import os
import glob
import pdb
import datetime


yyyymmdd = ['20151112']

#hourlist = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']
hourlist = ['21','22']

def get_time_begin(tstr):
    time_begin = datetime.datetime(int(tstr[0:4]),int(tstr[5:7]),int(tstr[8:10]),
                                   int(tstr[11:13]),int(tstr[14:16]))
    return time_begin

def load_dow(time_begin):
    yyyymmdd = time_begin.strftime('%Y%m%d')
    dd = yyyymmdd[-2:]
    hhmm = time_begin.strftime('%H%M')
    fileloc = '/home/disk/bob/olympex/cfradial/moments/dow_lo_qc/rhi/'+yyyymmdd+'/'
    fileind = fileloc+'cfrad.'+currentdate+'_'+hhmm[0:2]+'*_RHI.nc' #
    filelist = glob.glob(fileind)
    filemin = []
    for ff in filelist:
        filemin.append(int(ff.split('_')[6][2:4]))
    filemin = np.array(filemin)
    filediff = abs(int(hhmm[2:4])-filemin)
    fileloc = np.where((filediff == np.min(filediff)))[0][0]
    dow_file = filelist[fileloc]
    return dow_file

def plot_rhi(nc_file,currentdate):
    
    fields_to_plot = ['DBZ','VEL_DEALIAS'] #VEL_DEALIAS
    axis_titles = ['Reflectivity (dBZ)','Radial Velocity (m s$^{-1}$)'] #Spectrum Width (m s$^{-1}$)
    titles = [' ',' ']

    dow_fields_to_plot = ['DBZHC_F','VEL_F']#'RHOHV']
    dow_axis_titles = ['','']
    dow_titles = ['','']

    ranges = [(0,50),(-40,40)]
    subplots = [211,212]
    nplots = len(fields_to_plot)
    fig = plt.figure()
    fig.set_size_inches(6,8)
    cmaps = ['pyart_Carbone42','pyart_Carbone42']
    
    radar = pyart.io.read(nc_file)
    #Manual NPOL dBZ adjustment     
    #radar.fields['DBZ']['data']-=4.5
    
    try:
        print 'Started dealias'
        dealias_data = pyart.correct.unwrap.dealias_unwrap_phase(radar,vel_field='VEL')
        print 'Finished dealias'
        radar.add_field('VEL_DEALIAS', dealias_data)
        print 'Added field'
    except:
        print "Timed out!"
     
    time_begin = get_time_begin(radar.metadata['start_datetime'])
    dow_file = load_dow(time_begin)
    dow_radar = pyart.io.read(dow_file)
    dow_time_begin = get_time_begin(radar.metadata['start_datetime'])

    display = pyart.graph.RadarMapDisplay(radar)
    print 'Created display'
    azimuths = radar.fixed_angle['data']
    ga = np.where((azimuths >= 51) & (azimuths <= 53))[0]

    dow_display = pyart.graph.RadarDisplay(dow_radar,shift=(35000,0))
    dow_azimuths = dow_radar.fixed_angle['data']

    for rr,ray in enumerate(radar.sweep_number['data'][ga]):
    #for ray in sweep_manual[ga]:
    # plot each field
        dow_ga = np.where((dow_azimuths >= 58.0) & (dow_azimuths < 59.0))[0]
        dow_sweep = dow_radar.sweep_number['data'][dow_ga][0]

        for ind,plot_num in enumerate(xrange(nplots)):
            ax = fig.add_subplot(subplots[ind])
            field = fields_to_plot[plot_num]  
            vmin, vmax = ranges[plot_num]

            npol_plot = display.plot_rhi(field, ray,title = titles[ind], cmap=cmaps[ind], colorbar_label = axis_titles[ind],
                             vmin=vmin, vmax=vmax, axislabels=('Distance from NPOL (km)','Height above NPOL (km)'),zorder=1)
            #Plot DOW RHI on top of NPOL:
            '''
            self, field, sweep=0, mask_tuple=None,
            vmin=None, vmax=None, norm=None, cmap=None,
            mask_outside=False, title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True,
            reverse_xaxis=None, colorbar_flag=True, colorbar_label=None,
            colorbar_orient='vertical', edges=True, gatefilter=None,
            filter_transitions=True, ax=None, fig=None,
            ticks=None, ticklabs=None, **kwargs):
            ''' 
            dow_field = dow_fields_to_plot[plot_num]
            x, y, z = dow_display._get_x_y_z(dow_sweep, edges=True, filter_transitions=True)
            dowdata = dow_display._get_data(dow_field, dow_sweep, mask_tuple = ('DBZHC_F',0), filter_transitions=True, gatefilter=None) 
            dowdata[-130:,:] = ma.masked            

            R = np.sqrt(x ** 2 + y ** 2) * np.sign(y)

            #dowdata[np.where((R[0:-1,0:-1] > 65))]  = ma.masked #60 is 35 (NPOL-DOW distance) + 30 km
            #dowdata[0:10,:] = ma.masked
            
            #Adjust dow data--don't adjust for now, probably correct
            #if ind == 0:
            #    dowdata += 5
            pm = ax.pcolormesh(R, z+0.1, dowdata, vmin=vmin, vmax=vmax, cmap=cmaps[ind], norm=None,zorder=2)
            time_text = time_begin.strftime('%d-%b %H:%M UTC')
            ax.text(0.99, 0.95,time_text, horizontalalignment='right',
            verticalalignment='center',
            transform=ax.transAxes)

            display.set_limits(ylim=[0, 10],xlim=[0,80])
            ax.grid(which='major')

             
        
        #plt.subplots_adjust(wspace = 0.10, hspace = 0.12) #clean out some of the whitespace
   

        # set the figure title and show
        time_text = time_begin.strftime('%Y-%b-%d %H:%M UTC')
        azimuth = radar.fixed_angle['data'][ray]
        #radar_name = display.radar_name
        radar_name = 'NPOL '
        #title = 'RHI ' + radar_name + time_text + ' Azimuth %.1f' % (azimuth)
         
        #title2 = title.replace('T',' ').replace('NPOL1',' NPOL ').replace('Z',' UTC ') #fix title string
        #plt.suptitle(title+'$^{\circ}$', fontsize = 20)   

        #plt.tight_layout()
        os.chdir("/home/disk/meso-home/jzagrod/Misc/For_Stacy/NPOL_DOW/No_Correction")
        os.system('mkdir -p ' + currentdate)
        outdir = '/home/disk/meso-home/jzagrod/Misc/For_Stacy/NPOL_DOW/No_Correction/'+currentdate+'/'
        os.chdir(outdir)
        os.system('mkdir -p ' + str(round(azimuth,1))) #arrange dir by azimuth
        outdir_azi = outdir+str(round(azimuth,1))+'/'
        plottime = time_begin.strftime(format='%H%M%S')
        plt.savefig(outdir_azi+'olympex_npol_qc2_RHI_'+currentdate+'_'+plottime+'_'+str(round(azimuth,1))+'deg.png', bbox_inches='tight')
        print 'Plotted Date:'+currentdate+' Time:'+plottime+' UTC Azimuth:'+str(round(azimuth,1))
        plt.clf()


 

for currentdate in yyyymmdd:
    fileloc = '/home/disk/bob/olympex/cfradial/moments/npol_qc2/rhi/'+currentdate+'/' #QC2
    for hour in hourlist:
        fileind = fileloc+'cfrad.'+currentdate+'_'+hour+'*RHI_east.nc' #A is over ocean
        filelist = glob.glob(fileind) #reverse order
        for nc_file in filelist:
            if 1 == 1:
                plot_rhi(nc_file,currentdate)
                plt.close('all')
            if 1 == 2:
                print 'failure'



