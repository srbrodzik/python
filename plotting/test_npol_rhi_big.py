# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:11:08 2015

@author: jzagrod

test_npol_rhi.py
Test plot of NPOL RHI scans
Seeing if I can get it to read the formats that we have...
Update 6/16/2015: Try with different color table
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyart
import os
import glob
import pdb

yyyymmdd = ['20151113','20151114','20151203']
hourlist = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']
#hourlist = ['20','21','22','23']
hourlist = ['10','11','12','13','14','15']

def plot_rhi(nc_file,currentdate):
    
    fields_to_plot = ['DBZ','VEL','ZDR','RHOHV']
    axis_titles = ['Reflectivity (dBZ)','Radial Velocity (m s$^{-1}$)','Differential Reflectivity (dB)',
                  'Correlation coefficient'] #Spectrum Width (m s$^{-1}$)
    titles = ['Reflectivity','Velocity','Z$_{\mathregular{DR}}$',
              r'$\mathregular{\rho}_{\mathregular{HV}}$']
    ranges = [(-10,50),(-30,30),(-1,2),(0.6,1)]#(0,5)]
    subplots = [221,222,223,224]
    nplots = len(fields_to_plot)
    fig = plt.figure()
    fig.set_size_inches(30,15)
    cmaps = ['pyart_Carbone42','pyart_Carbone42',None,'pyart_Carbone42']
    
    radar = pyart.io.read(nc_file)

    #Dealias velocity...move on if code hangs up for > 10 seconds
    '''
    try:
        print 'Startind dealias'
        dealias_data = pyart.correct.unwrap.dealias_unwrap_phase(radar,vel_field='VEL')
        print 'Finished dealias'
        radar.add_field('VEL_DEALIAS', dealias_data)
        print 'Added field'
    except:
        print "Timed out!"
    '''  



    display = pyart.graph.RadarMapDisplay(radar)
    print 'Created display'
    azimuths = radar.fixed_angle['data']
    #sweep_manual = [4,5,9,10,11,12,13,14,15] #to save time, only plot 2nd ray
    #plot each azimuth
    for ray in radar.sweep_number['data']:
    #for ray in sweep_manual:
    # plot each field
        for ind,plot_num in enumerate(xrange(nplots)):
            ax = fig.add_subplot(subplots[ind])
            field = fields_to_plot[plot_num]  
            vmin, vmax = ranges[plot_num]
            display.plot_rhi(field, ray,title = titles[ind], cmap=cmaps[ind], colorbar_label = axis_titles[ind],
                             vmin=vmin, vmax=vmax, axislabels=('Distance from radar (km)','Height above radar (km)'))
            display.set_limits(ylim=[0, 10],xlim=[0,100])
            ax.grid(which='major')
             
        
        #plt.subplots_adjust(wspace = 0.10, hspace = 0.12) #clean out some of the whitespace
   

        # set the figure title and show
        time_text = display.time_begin.strftime('%Y-%b-%d %H:%M UTC')
        azimuth = radar.fixed_angle['data'][ray]
        #radar_name = display.radar_name
        radar_name = 'NPOL '
        title = 'RHI ' + radar_name + time_text + ' Azimuth %.1f' % (azimuth)
        #title2 = title.replace('T',' ').replace('NPOL1',' NPOL ').replace('Z',' UTC ') #fix title string
        plt.suptitle(title+'$^{\circ}$', fontsize = 20)   

        #plt.tight_layout()
        os.chdir("/home/disk/meso-home/jzagrod/Olympex/NPOL/RHI/RHI_4panel_100km")
        os.system('mkdir -p ' + currentdate)
        outdir = '/home/disk/meso-home/jzagrod/Olympex/NPOL/RHI/RHI_4panel_100km/'+currentdate+'/'
        os.chdir(outdir)
        os.system('mkdir -p ' + str(round(azimuth,2))) #arrange dir by azimuth
        outdir_azi = outdir+str(round(azimuth,2))+'/'
        plottime = display.time_begin.strftime(format='%H%M%S')
        plt.savefig(outdir_azi+'olympex_npol_RHI_'+currentdate+'_'+plottime+'_'+str(round(azimuth,2))+'deg.png', bbox_inches='tight')
        print 'Plotted Date:'+currentdate+' Time:'+plottime+' UTC Azimuth:'+str(round(azimuth,2))
        plt.clf()


 

for currentdate in yyyymmdd:
    fileloc = '/home/disk/bob/olympex/cfradial/moments/npol/rhi/'+currentdate+'/'
    for hour in hourlist:
        fileind = fileloc+'cfrad.'+currentdate+'_'+hour+'*B_RHI.nc' #A is over ocean
        filelist = glob.glob(fileind) #reverse order
        for nc_file in filelist:
            if 1 == 1:
                plot_rhi(nc_file,currentdate)
                plt.close('all')
            if 1 == 2:
                print 'failure'



