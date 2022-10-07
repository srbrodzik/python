'''
Read an APR-3 HDF file in Python
Plot 5-panel with 3-frequency and 2 maps
Plot 5-panel with frequency differences and 2 maps

Original template by Randy Chase, Univ of Illinois
Modified by Joe Zagrodnik, UW, Oct 2017
'''
from pyhdf.SD import SD, SDC
from matplotlib import dates
import pyart
from pyart.graph import cm
from mpl_toolkits.basemap import Basemap, maskoceans
import numpy as np
import datetime
import pdb
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.dates import MinuteLocator, DateFormatter
import os

#custom colormap
from matplotlib.colors import LinearSegmentedColormap
cMap = []
for value, colour in zip([0,75,275,425,575,1000],['#809980','#80B280','#EAD480','#E0C08C','#BAA38C','#FFFFFF']):
    cMap.append((value/1000.0, colour))
JoeTerrainFade = LinearSegmentedColormap.from_list("custom", cMap)

#topographic map
indir = '/home/disk/funnel/olympex/archive2/wxdata/Plotcode/topo/'
infile = indir+'topo_arrays_final.npz'
indata = np.load(infile)

lons = indata['lons']
lats = indata['lats']
topoin = indata['topoin']

def apr3read(filename):
    
    """
    ===========
    
    This is for reading in apr3 hdf files from OLYMPEX and return them all in one dictionary
    
    ===========
    
    filename = filename of the apr3 file
    """
    
    apr = {}
    flag = 0
    ##Radar varibles in hdf file found by hdf.datasets
    radar_freq = 'zhh14' #Ku
    radar_freq2 = 'zhh35' #Ka
    radar_freq3 = 'zvv95' #W
    radar_freq4 = 'ldr14' #LDR
    vel_str = 'vel14' #Doppler
    ##

    hdf = SD(filename, SDC.READ)
    
    listofkeys = hdf.datasets().keys()

    alt = hdf.select('alt3D')
    lat = hdf.select('lat')
    lon = hdf.select('lon')
    time = hdf.select('scantime').get()
    surf = hdf.select('surface_index').get()
    isurf = hdf.select('isurf').get()
    plane = hdf.select('alt_nav').get()
    radar = hdf.select(radar_freq)
    radar2 = hdf.select(radar_freq2)
    radar4 = hdf.select(radar_freq4)
    vel = hdf.select(vel_str)
    lon3d = hdf.select('lon3D')
    lat3d = hdf.select('lat3D')
    alt3d = hdf.select('alt3D')
    lat3d_scale = hdf.select('lat3D_scale').get()[0][0]
    lon3d_scale = hdf.select('lon3D_scale').get()[0][0]
    alt3d_scale = hdf.select('alt3D_scale').get()[0][0]
    lat3d_offset = hdf.select('lat3D_offset').get()[0][0]
    lon3d_offset = hdf.select('lon3D_offset').get()[0][0]
    alt3d_offset = hdf.select('alt3D_offset').get()[0][0]
 
    alt = alt.get()
    ngates = alt.shape[0]
    lat = lat.get()
    lon = lon.get()
    
    lat3d = lat3d.get()
    lat3d = (lat3d/lat3d_scale) + lat3d_offset
    lon3d = lon3d.get()
    lon3d = (lon3d/lon3d_scale) + lon3d_offset
    alt3d = alt3d.get()
    alt3d = (alt3d/alt3d_scale) + alt3d_offset
    
    radar_n = radar.get()
    radar_n = radar_n/100.
    radar_n2 = radar2.get()
    radar_n2 = radar_n2/100.
    
    #Search for W vv
    listofkeys = hdf.datasets().keys()
    if 'zvv95' in listofkeys:
        if 'zhh95' in listofkeys:
            radar_nadir = hdf.select('zvv95')
            radar_scanning = hdf.select('zhh95')
            radar_nadir = radar_nadir.get()
            radar_scanning = radar_scanning.get()
            radar_n3 = radar_scanning
            ##uncomment if you want high sensativty as nadir scan (WARNING, CALIBRATION)
            #radar_n3[:,12,:] = radar_nadir[:,12,:]
            radar_n3 = radar_n3/100.
        else:
            radar3 = hdf.select('zvv95')
            radar_n3 = radar3.get()
            radar_n3 = radar_n3/100.
            print('No hh, using just nadir vv')
    elif 'zhh95' in listofkeys:
        radar3 = hdf.select('zhh95')
        radar_n3 = radar3.get()
        radar_n3 = radar_n3/100.
        print('No vv, using hh')
    else:
        radar_n3 = np.array([])
        flag = 1
        print('No W band')
        
    radar_n4 = radar4.get()
    radar_n4 = radar_n4/100.
    vel_n = vel.get()
    vel_n = vel_n/100.
    

    #Quality control (masked where invalid
    radar_n = np.ma.masked_where(radar_n <= -99.99,radar_n)
    radar_n2 = np.ma.masked_where(radar_n2 <= -99.99,radar_n2)
    radar_n3 = np.ma.masked_where(radar_n3 <= -99.99,radar_n3)
    radar_n4 = np.ma.masked_where(radar_n4 <= -99.99,radar_n4)
    
    ##convert time to datetimes
    time_dates = np.empty(time.shape,dtype=object)
    for i in np.arange(0,time.shape[0]):
        for j in np.arange(0,time.shape[1]):
            tmp = datetime.datetime.utcfromtimestamp(time[i,j])
            time_dates[i,j] = tmp
            
    #Create a time at each gate (assuming it is the same down each ray, there is a better way to do this)      
    time_gate = np.empty(lat3d.shape,dtype=object)
    for k in np.arange(0,550):
        for i in np.arange(0,time_dates.shape[0]):
            for j in np.arange(0,time_dates.shape[1]):
                time_gate[k,i,j] = time_dates[i,j]        

    apr['Ku'] = radar_n
    apr['Ka'] = radar_n2
    apr['W'] = radar_n3
    apr['DFR_1'] = radar_n - radar_n2 #Ku - Ka
    
    if flag == 0:
        apr['DFR_3'] = radar_n2 - radar_n3 #Ka - W
        apr['DFR_2'] = radar_n - radar_n3 #Ku - W
        apr['info'] = 'The shape of these arrays are: Radar[Vertical gates,Time/DistanceForward]'
    else:
        apr['DFR_3'] = np.array([]) #Ka - W
        apr['DFR_2'] = np.array([]) #Ku - W
        apr['info'] = 'The shape of these arrays are: Radar[Vertical gates,Time/DistanceForward], Note No W band avail'
        
    apr['ldr'] = radar_n4
    apr['vel'] = vel_n
    apr['lon'] = lon
    apr['lat'] = lat
    apr['alt_gate'] = alt3d
    apr['alt_plane'] = plane
    apr['surface'] = isurf 
    apr['time']= time
    apr['timedates']= time_dates
    apr['time_gate'] = time_gate
    apr['lon_gate'] = lon3d
    apr['lat_gate'] = lat3d
    
    fileheader = hdf.select('fileheader')
    roll = hdf.select('roll').get()
    pitch = hdf.select('pitch').get()
    drift = hdf.select('drift').get()
    
    apr['fileheader'] = fileheader
    apr['ngates'] = ngates
    apr['roll'] = roll
    apr['pitch'] = pitch
    apr['drift'] = drift
    
    _range = np.arange(15,550*30,30)
    _range = np.asarray(_range,float)
    ind = np.where(_range >= plane.mean())
    _range[ind] = np.nan
    apr['range'] = _range
    
    return apr

def get_nexrad_ppis(aprtime):
    #get two display objects for Langley and Camano PPIs
    #nearest time to the flight track
    yyyymmddhhmm = aprtime.strftime('%Y%m%d%H%M')
    radardir = '/home/disk/bob/olympex/cfradial/moments/'
    lgx_files = glob.glob('{}klgx/{}/*'.format(radardir,yyyymmddhhmm[0:8]))
    atx_files = glob.glob('{}katx/{}/*'.format(radardir,yyyymmddhhmm[0:8]))
    lgx_datetimes = []
    atx_datetimes = []
    for fi in lgx_files:
        tstr = fi.split('/')[9].split('.')[1].split('_')[1]
        lgx_datetimes.append(datetime.datetime(int(yyyymmddhhmm[0:4]),int(yyyymmddhhmm[4:6]),int(yyyymmddhhmm[6:8]),
                             int(tstr[0:2]),int(tstr[2:4])))
    for fi in atx_files:
        tstr = fi.split('/')[9].split('.')[1].split('_')[1]
        atx_datetimes.append(datetime.datetime(int(yyyymmddhhmm[0:4]),int(yyyymmddhhmm[4:6]),int(yyyymmddhhmm[6:8]),
                             int(tstr[0:2]),int(tstr[2:4])))
    lgx_datetimes = np.array(lgx_datetimes)
    atx_datetimes = np.array(atx_datetimes)
    lgx_diff = abs(lgx_datetimes - aprtime)
    lgx_loc = np.where((lgx_diff == np.min(lgx_diff)))[0][0]
    atx_diff = abs(atx_datetimes - aprtime)
    atx_loc = np.where((atx_diff == np.min(atx_diff)))[0][0]
    #load pyart
    lgxradar = pyart.io.read(lgx_files[lgx_loc])
    atxradar = pyart.io.read(atx_files[atx_loc])
    lgxdisplay = pyart.graph.RadarMapDisplay(lgxradar)
    atxdisplay = pyart.graph.RadarMapDisplay(atxradar)
    return lgxdisplay,atxdisplay

def apr3plot(filename):
    apr = apr3read(filename)
    starttime = apr['timedates'][12,0]
    endtime = apr['timedates'][12,-1]

    #make 2d array of times
    time1d = np.repeat(apr['timedates'][12,:],550,axis=0)
    xlen = np.shape(apr['timedates'])[1]
    ylen = np.shape(apr['Ku'])[0]
    time2d = time1d.reshape((xlen,ylen)).T

    fig = plt.figure()
    fig.set_size_inches(15,8.5)

    #ax0 = plt.subplot2grid((3, 4), (0, 0), colspan=3)
    #ax1 = plt.subplot2grid((3, 4), (1, 0), colspan=3)
    #ax2 = plt.subplot2grid((3, 4), (2, 0), colspan=3)
    #ax3 = plt.subplot2grid((3, 4), (0, 3), rowspan=3)

    ax0 = plt.subplot2grid((6, 4), (0, 0), colspan=3,rowspan=2)
    ax1 = plt.subplot2grid((6, 4), (2, 0), colspan=3,rowspan=2)
    ax2 = plt.subplot2grid((6, 4), (4, 0), colspan=3,rowspan=2)
    ax3 = plt.subplot2grid((6, 4), (0, 3), rowspan=3)
    ax4 = plt.subplot2grid((6, 4), (3, 3), rowspan=3)

    try:
        W_band_data = apr['W'][:,12,:]
    except:
        W_band_data = np.zeros(np.shape(apr['Ku'][:,12,:]))
        W_band_data[:,:] = float('nan')

    #ax0 = fig.add_subplot(311)
    ax0.set_title('OLYMPEX APR-3 Reflectivity (dBZ)  {} - {} UTC'.format(starttime.strftime('%Y-%b-%d %H:%M'),endtime.strftime('%H:%M')),size=15)
    plt0 = ax0.pcolor(time2d,apr['alt_gate'][:,12,:]/1000.,apr['Ku'][:,12,:],cmap='pyart_Carbone42',vmin=-20,vmax=50)
    #ax1 = fig.add_subplot(312)
    plt1 = ax1.pcolor(time2d,apr['alt_gate'][:,12,:]/1000.,apr['Ka'][:,12,:],cmap='pyart_Carbone42',vmin=-20,vmax=50)
    #ax2 = fig.add_subplot(313)
    plt2 = ax2.pcolor(time2d,apr['alt_gate'][:,12,:]/1000.,W_band_data,cmap='pyart_Carbone42',vmin=-20,vmax=50)
    axes = [ax0,ax1,ax2]
    pltnames = [plt0,plt1,plt2]
    pltlabels = ['(a) Ku','(b) Ka','(c) W']
    for i,ax in enumerate(axes):
        ax.set_ylim(0,12)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="1%", pad=0.15)
        cbar = plt.colorbar(pltnames[i], cax=cax)
        cbar.ax.tick_params(labelsize=10) 
        #cbar.set_label('Reflectivity (dBZ)', fontsize=10)
        ax.annotate(pltlabels[i], xy=(0.02, 0.96), xycoords='axes fraction', fontsize=16,
                horizontalalignment='left', verticalalignment='top')   
        ax.set_ylabel('Height (km)')
	ax.xaxis.set_major_locator( MinuteLocator() )
	ax.xaxis.set_major_formatter( DateFormatter('%H:%M') )

    fontsize=12
    lat = apr['lat'][12,:]
    lon = apr['lon'][12,:]
    corners = [np.around(np.min(lon)-0.2,decimals=1),np.around(np.min(lat)-0.2,decimals=1),
               np.around(np.max(lon)+0.2,decimals=1),np.around(np.max(lat)+0.2,decimals=1)]
    #merdians = np.linspace(corners[0]+1,corners[2]-1,3)
    #parallels = np.linspace(corners[1],corners[3],5)

    

    #m = Basemap(llcrnrlon=corners[0],llcrnrlat=corners[1],urcrnrlon=corners[2],urcrnrlat=corners[3],ax=ax3,resolution='h',projection='cass')
    m = Basemap(
		      lon_0 = -124.829133,
		      lat_0 = 47.653729,
		      width = 380000,
		      height = 320000,ax=ax3,resolution='h',projection='cass')

    #basemap for zoomed flight track
    m2 = Basemap(lon_0 = np.mean(lon), lat_0 = np.mean(lat), llcrnrlon=corners[0],llcrnrlat=corners[1],urcrnrlon=corners[2],urcrnrlat=corners[3],ax=ax4,resolution='h',projection='cass')
    
    nx = int((m2.xmax-m2.xmin)/100.)+1; ny = int((m2.ymax-m2.ymin)/100.)+1
    topodat = m2.transform_scalar(topoin,lons,lats,nx,ny)
    im = m2.imshow(topodat,'Greys',vmin = 0, vmax = 2500,zorder=-1) #cm.GMT_haxby vmax = 2500
  
    #coastline...
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch
    patches = []
    shapedir2 = '/home/disk/meso-home/jzagrod/Olympex/RiverBasins/Olympex_basins/Coast/'
    shapefile2 = shapedir2+'pacnw_water'
    #m2.readshapefile(shapefile2,'pacnw_water', drawbounds=False)
    #for info, shape in zip(m2.pacnw_water, m2.pacnw_water):
    #         patches.append(Polygon(np.array(shape), True) )
    #ax4.add_collection(PatchCollection(patches, facecolor= None, edgecolor='k', linewidths=0.8,zorder=1))

   #outline of lake quinault
    shapedir2 = '/home/disk/meso-home/jzagrod/Olympex/RiverBasins/Olympex_basins/LakeQuinault/'
    shapefile2 = shapedir2+'lake_quinault'
    m2.readshapefile(shapefile2,'lake_quinault', drawbounds = False)
    patches   = []
    for info, shape in zip(m2.lake_quinault_info, m2.lake_quinault):
            patches.append( Polygon(np.array(shape), True) )
    ax4.add_collection(PatchCollection(patches, facecolor= 'white', edgecolor='k', linewidths=0.8,zorder=1))

    xnp,ynp = m2(-124.21,47.27)
    m2.scatter(xnp, ynp, marker='D',color='orange',s=30,edgecolor = 'k',zorder=10)

    xnp,ynp = m2(-123.498,47.97)
    m2.scatter(xnp, ynp, marker='D',color='red',s=30,edgecolor = 'k',zorder=10)

    #Get plot dot at lat/lon of each minute in the track
    alltimes = apr['timedates'][12,:]
    allmins = np.array([y.minute for y in alltimes])
    allsecs = np.array([y.second for y in alltimes])
    #delete beginning elements if second > 3
    if allsecs[0] > 3:
        alltimes = np.delete(alltimes,np.where((allmins == allmins[0]))[0])
        allsecs = np.delete(allsecs,np.where((allmins == allmins[0]))[0])
        lonind = np.delete(lon,np.where((allmins == allmins[0]))[0])
        latind = np.delete(lat,np.where((allmins == allmins[0]))[0])
        allmins = np.delete(allmins,np.where((allmins == allmins[0]))[0])
    else:
        lonind = lon
        latind = lat
    uniq_mins, min_indices = np.unique(allmins,return_index=True)

    #plot radars
    try:
        lgx_disp,atx_disp = get_nexrad_ppis(alltimes[len(alltimes)/2])
        atxplt = atx_disp.plot_ppi_map('DBZ',0,vmin=-20, vmax=50, cmap="pyart_Carbone42", 
                         title=None, mask_outside = True, colorbar_flag = False, ax=ax3, projection='cass',
                         basemap = m,zorder=2)
        lgxplt = lgx_disp.plot_ppi_map('DBZ',1,vmin=-20, vmax=50, cmap="pyart_Carbone42", 
                         title=None, mask_outside = True, colorbar_flag = False, ax=ax3, projection='cass',
                         basemap = m,zorder=2)
    except:
        print 'No ground radars available'


    #m.fillcontinents(color='#e6e6e6',zorder=1)
    #ax3.text('KLGX+KATX 0.5$^\circ$ PPI (dBZ)')

    mapaxes = [ax3,ax4]
    maps = [m,m2]
    linewidths = [1,2]
    markersizes = [10,30]
    ax3.set_title('Wide Flight Track',size=15)
    ax4.set_title('Zoomed Flight Track',size=15)
    for ii,mm in enumerate(maps):

        mm.drawmapboundary()
        parallels = np.arange(0.,90,0.5)
        mm.drawparallels(parallels,labels=[0,1,0,0],fontsize=10, linewidth=1) #0 hides them
        meridians = np.arange(180.,360.,1)
        mm.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10, linewidth=1)
        mm.drawcoastlines(linewidth=1,color='k')
        x,y = mm(lon,lat)
        mm.plot(x,y,'-',color='navy',lw=linewidths[ii],zorder=5)

        #loop through unique indices and plot an 'X' at the first one and a dot at others
        for ind,mn in enumerate(uniq_mins):
            x1,y1 = mm(lonind[min_indices[ind]],latind[min_indices[ind]])
            if min_indices[ind] == 0:
                label = alltimes[min_indices[ind]].strftime('%H:%M UTC')
                mm.scatter(x1,y1,marker='X',color='navy',label=label,s=markersizes[ii]*2.2,zorder=5)
            else:
                mm.scatter(x1,y1,marker='o',color='navy',s=markersizes[ii],zorder=5)

        mapaxes[ii].legend(prop={'size':10},loc="upper left",scatterpoints=1)
    outdate = starttime.strftime('%Y%m%d')
    outdatehhmm = starttime.strftime('%Y%m%d_%H%M%S')
    outdatestacy = starttime.strftime('%Y%m%d%H%M%S')
    os.chdir(outdir)
    os.system('mkdir -p ' + outdate)
    os.chdir(outdir+outdate)
    plt.tight_layout(h_pad=0.3)
    plt.savefig('APR3_{}.png'.format(outdatehhmm), bbox_inches='tight',dpi=200)
    stacyfile = '/home/disk/funnel/olympex/archive/research/dc8/{}/research.dc8.{}.apr3.png'.format(outdate,outdatestacy)
    plt.savefig(stacyfile, bbox_inches='tight')

    print 'Plotted {}'.format(outdatehhmm)
    plt.close('all')

def apr3plot_diff(filename):
    apr = apr3read(filename)
    starttime = apr['timedates'][12,0]
    endtime = apr['timedates'][12,-1]

    #make 2d array of times
    time1d = np.repeat(apr['timedates'][12,:],550,axis=0)
    xlen = np.shape(apr['timedates'])[1]
    ylen = np.shape(apr['Ku'])[0]
    time2d = time1d.reshape((xlen,ylen)).T

    fig = plt.figure()
    fig.set_size_inches(15,8.5)

    #ax0 = plt.subplot2grid((3, 4), (0, 0), colspan=3)
    #ax1 = plt.subplot2grid((3, 4), (1, 0), colspan=3)
    #ax2 = plt.subplot2grid((3, 4), (2, 0), colspan=3)
    #ax3 = plt.subplot2grid((3, 4), (0, 3), rowspan=3)

    ax0 = plt.subplot2grid((6, 4), (0, 0), colspan=3,rowspan=2)
    ax1 = plt.subplot2grid((6, 4), (2, 0), colspan=3,rowspan=2)
    ax2 = plt.subplot2grid((6, 4), (4, 0), colspan=3,rowspan=2)
    ax3 = plt.subplot2grid((6, 4), (0, 3), rowspan=3)
    ax4 = plt.subplot2grid((6, 4), (3, 3), rowspan=3)

    try:
        W_band_data = apr['W'][:,12,:]
    except:
        W_band_data = np.zeros(np.shape(apr['Ku'][:,12,:]))
        W_band_data[:,:] = float('nan')

    #ax0 = fig.add_subplot(311)
    ax0.set_title('OLYMPEX APR-3 Reflectivity (dBZ)  {} - {} UTC'.format(starttime.strftime('%Y-%b-%d %H:%M'),endtime.strftime('%H:%M')),size=15)
    plt0 = ax0.pcolor(time2d,apr['alt_gate'][:,12,:]/1000.,apr['Ku'][:,12,:]-apr['Ka'][:,12,:],cmap='seismic',vmin=-10,vmax=10)
    #ax1 = fig.add_subplot(312)
    plt1 = ax1.pcolor(time2d,apr['alt_gate'][:,12,:]/1000.,apr['Ka'][:,12,:]-W_band_data,cmap='seismic',vmin=-20,vmax=20)
    #ax2 = fig.add_subplot(313)
    plt2 = ax2.pcolor(time2d,apr['alt_gate'][:,12,:]/1000.,apr['Ku'][:,12,:]-W_band_data,cmap='seismic',vmin=-30,vmax=30)
    axes = [ax0,ax1,ax2]
    pltnames = [plt0,plt1,plt2]
    pltlabels = ['(a) Ku-Ka','(b) Ka-W','(c) Ku-W']
    for i,ax in enumerate(axes):
        ax.set_ylim(0,12)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="1%", pad=0.15)
        cbar = plt.colorbar(pltnames[i], cax=cax)
        cbar.ax.tick_params(labelsize=10) 
        #cbar.set_label('Reflectivity (dBZ)', fontsize=10)
        ax.annotate(pltlabels[i], xy=(0.02, 0.96), xycoords='axes fraction', fontsize=16,
                horizontalalignment='left', verticalalignment='top')   
        ax.set_ylabel('Height (km)')
	ax.xaxis.set_major_locator( MinuteLocator() )
	ax.xaxis.set_major_formatter( DateFormatter('%H:%M') )

    fontsize=12
    lat = apr['lat'][12,:]
    lon = apr['lon'][12,:]
    corners = [np.around(np.min(lon)-0.2,decimals=1),np.around(np.min(lat)-0.2,decimals=1),
               np.around(np.max(lon)+0.2,decimals=1),np.around(np.max(lat)+0.2,decimals=1)]
    #merdians = np.linspace(corners[0]+1,corners[2]-1,3)
    #parallels = np.linspace(corners[1],corners[3],5)

    

    #m = Basemap(llcrnrlon=corners[0],llcrnrlat=corners[1],urcrnrlon=corners[2],urcrnrlat=corners[3],ax=ax3,resolution='h',projection='cass')
    m = Basemap(
		      lon_0 = -124.829133,
		      lat_0 = 47.653729,
		      width = 380000,
		      height = 320000,ax=ax3,resolution='h',projection='cass')

    #basemap for zoomed flight track
    m2 = Basemap(lon_0 = np.mean(lon), lat_0 = np.mean(lat), llcrnrlon=corners[0],llcrnrlat=corners[1],urcrnrlon=corners[2],urcrnrlat=corners[3],ax=ax4,resolution='h',projection='cass')
    
    nx = int((m2.xmax-m2.xmin)/100.)+1; ny = int((m2.ymax-m2.ymin)/100.)+1
    topodat = m2.transform_scalar(topoin,lons,lats,nx,ny)
    im = m2.imshow(topodat,'Greys',vmin = 0, vmax = 2500,zorder=-1) #cm.GMT_haxby vmax = 2500
  
    #coastline...
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch
    patches = []
    shapedir2 = '/home/disk/meso-home/jzagrod/Olympex/RiverBasins/Olympex_basins/Coast/'
    shapefile2 = shapedir2+'pacnw_water'
    #m2.readshapefile(shapefile2,'pacnw_water', drawbounds=False)
    #for info, shape in zip(m2.pacnw_water, m2.pacnw_water):
    #         patches.append(Polygon(np.array(shape), True) )
    #ax4.add_collection(PatchCollection(patches, facecolor= None, edgecolor='k', linewidths=0.8,zorder=1))

   #outline of lake quinault
    shapedir2 = '/home/disk/meso-home/jzagrod/Olympex/RiverBasins/Olympex_basins/LakeQuinault/'
    shapefile2 = shapedir2+'lake_quinault'
    m2.readshapefile(shapefile2,'lake_quinault', drawbounds = False)
    patches   = []
    for info, shape in zip(m2.lake_quinault_info, m2.lake_quinault):
            patches.append( Polygon(np.array(shape), True) )
    ax4.add_collection(PatchCollection(patches, facecolor= 'white', edgecolor='k', linewidths=0.8,zorder=1))

    xnp,ynp = m2(-124.21,47.27)
    m2.scatter(xnp, ynp, marker='D',color='orange',s=30,edgecolor = 'k',zorder=10)

    xnp,ynp = m2(-123.498,47.97)
    m2.scatter(xnp, ynp, marker='D',color='red',s=30,edgecolor = 'k',zorder=10)

    #Get plot dot at lat/lon of each minute in the track
    alltimes = apr['timedates'][12,:]
    allmins = np.array([y.minute for y in alltimes])
    allsecs = np.array([y.second for y in alltimes])
    #delete beginning elements if second > 3
    if allsecs[0] > 3:
        alltimes = np.delete(alltimes,np.where((allmins == allmins[0]))[0])
        allsecs = np.delete(allsecs,np.where((allmins == allmins[0]))[0])
        lonind = np.delete(lon,np.where((allmins == allmins[0]))[0])
        latind = np.delete(lat,np.where((allmins == allmins[0]))[0])
        allmins = np.delete(allmins,np.where((allmins == allmins[0]))[0])
    else:
        lonind = lon
        latind = lat
    uniq_mins, min_indices = np.unique(allmins,return_index=True)

    #plot radars
    try:
        lgx_disp,atx_disp = get_nexrad_ppis(alltimes[len(alltimes)/2])
        atxplt = atx_disp.plot_ppi_map('DBZ',0,vmin=-20, vmax=50, cmap="pyart_Carbone42", 
                         title=None, mask_outside = True, colorbar_flag = False, ax=ax3, projection='cass',
                         basemap = m,zorder=2)
        lgxplt = lgx_disp.plot_ppi_map('DBZ',1,vmin=-20, vmax=50, cmap="pyart_Carbone42", 
                         title=None, mask_outside = True, colorbar_flag = False, ax=ax3, projection='cass',
                         basemap = m,zorder=2)
    except:
        print 'No ground radars available'


    #m.fillcontinents(color='#e6e6e6',zorder=1)
    #ax3.text('KLGX+KATX 0.5$^\circ$ PPI (dBZ)')

    mapaxes = [ax3,ax4]
    maps = [m,m2]
    linewidths = [1,2]
    markersizes = [10,30]
    ax3.set_title('Wide Flight Track',size=15)
    ax4.set_title('Zoomed Flight Track',size=15)
    for ii,mm in enumerate(maps):

        mm.drawmapboundary()
        parallels = np.arange(0.,90,0.5)
        mm.drawparallels(parallels,labels=[0,1,0,0],fontsize=10, linewidth=1) #0 hides them
        meridians = np.arange(180.,360.,1)
        mm.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10, linewidth=1)
        mm.drawcoastlines(linewidth=1,color='k')
        x,y = mm(lon,lat)
        mm.plot(x,y,'-',color='navy',lw=linewidths[ii],zorder=5)

        #loop through unique indices and plot an 'X' at the first one and a dot at others
        for ind,mn in enumerate(uniq_mins):
            x1,y1 = mm(lonind[min_indices[ind]],latind[min_indices[ind]])
            if min_indices[ind] == 0:
                label = alltimes[min_indices[ind]].strftime('%H:%M UTC')
                mm.scatter(x1,y1,marker='X',color='navy',label=label,s=markersizes[ii]*2.2,zorder=5)
            else:
                mm.scatter(x1,y1,marker='o',color='navy',s=markersizes[ii],zorder=5)

        mapaxes[ii].legend(prop={'size':10},loc="upper left",scatterpoints=1)
    outdate = starttime.strftime('%Y%m%d')
    outdatehhmm = starttime.strftime('%Y%m%d_%H%M%S')
    outdatestacy = starttime.strftime('%Y%m%d%H%M%S')
    os.chdir(outdir)
    os.system('mkdir -p ' + outdate)
    os.chdir(outdir+outdate)
    plt.tight_layout(h_pad=0.3)
    plt.savefig('APR3_{}_difference.png'.format(outdatehhmm), bbox_inches='tight',dpi=200)
    stacyfile = '/home/disk/funnel/olympex/archive/research/dc8/{}/research.dc8.{}.apr3_difference.png'.format(outdate,outdatestacy)
    #plt.savefig(stacyfile, bbox_inches='tight')

    print 'Plotted {}'.format(outdatehhmm)
    plt.close('all')
 

apr3_dates = ['20151202','20151203','20151204','20151205','20151208','20151210','20151212',
              '20151213','20151218','20151219','20151112']
for dd in apr3_dates:
    indir = '/home/disk/bob/olympex/apr3/{}/'.format(dd)
    stacydir = '/home/disk/funnel/olympex/archive/research/dc8/'
    outdir = '/home/disk/meso-home/jzagrod/Olympex/APR3/Plot/Diff/'
    infiles = sorted(glob.glob(indir+'*.HDF'))
    #aprdata = apr3read(infile)
    for infile in infiles:
        if 1 == 1:
            apr3plot_diff(infile)
        if 1 == 2:
            print 'Failure'
    #pdb.set_trace()
    #apr3flighttrack(infile)



pdb.set_trace()
