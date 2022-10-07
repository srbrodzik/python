from pyhdf.SD import SD, SDC
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib import dates
import pyresample as pr
from scipy.spatial import cKDTree
from pyproj import Proj

def apr3tocit(apr3filename,fl,sphere_size,query_k = 1,plotson=False):

    """
    =================
    
    This function finds either the closest gate or averages over a number of gates (query_k) nearest to 
    the citation aircraft in the radar volume of the APR3. It will return a dict of the original full length
    arrays and the matched arrays. 
    
    =================
    """
    
    cit_time = fl['time']['data']
    apr = apr3read(apr3filename)
    time_dates = apr['timedates'][:,:]
    
    if time_dates[12,:].shape[0] < 50:
        print('Limited radar gates in time')
        return
    
    fontsize=14
    #Varibles needed for the kdtree
    leafsize = 16
    query_eps = 0
    query_p=2
    query_distance_upper_bound = sphere_size
    query_n_jobs =1
    Barnes = True
    K_d = 1e3
    #


    #Pre-Determine arrays
    Ku_gate = np.array([])
    Ka_gate = np.array([])
    DFR_gate = np.array([])
    lon_c = np.array([])
    lat_c = np.array([])
    alt_c = np.array([])
    t_c = np.array([])
    lon_r = np.array([])
    lat_r = np.array([])
    alt_r = np.array([])
    t_r = np.array([])
    dis_r = np.array([])
    ind_r = np.array([])
    conc_hvps3 = np.array([])
    #

    
    #Set reference point (Currently Mount Olympus, Washington)
    lat_0 = 47.7998
    lon_0 = -123.7066
    #Set up map projection to calculate distances for this radar
    p = Proj(proj='laea', zone=10, ellps='WGS84',
             lat_0=lat_0,
             lon_0=lon_0)


    td = np.ravel(time_dates)
    datestart = td[0]
    dateend = td[td.shape[0]-1] 

    #Constrain Citation data to radar time
    ind = np.where(cit_time > datestart)
    ind2 = np.where(cit_time < dateend)
    ind3 = np.intersect1d(ind,ind2)
    #

    cit_time2 = fl['time']['data'][ind3]
    cit_lon = fl['longitude']['data'][ind3]
    cit_lat = fl['latitude']['data'][ind3]
    cit_alt = fl['altitude']['data'][ind3]
    bigins = fl['HVPS3_hor']['data'][ind3]
    
    #Print out number of potential points
    print(cit_time2.shape)
    #
    

    #Make 1-D arrays of radar spatial files
    apr_x = np.ravel(apr['lon_gate'][:,:,:])
    apr_y = np.ravel(apr['lat_gate'][:,:,:])
    apr_alt = np.ravel(apr['alt_gate'][:,:,:])
    apr_t = np.ravel(apr['time_gate'][:,:,:])
    #
    
    ##Make 1-D arrays of radar 
    apr_ku = np.ravel(apr['Ku'][:,:,:])
    apr_ka = np.ravel(apr['Ka'][:,:,:])
    apr_DFR = apr_ku - apr_ka
    ##
    
    #Use projection to get cartiesian distances
    apr_x2,apr_y2 = p(apr_x,apr_y)
    cit_x2,cit_y2 = p(cit_lon,cit_lat)
    #
    
    #Kdtree things
    kdt = cKDTree(zip(apr_x2, apr_y2, apr_alt), leafsize=leafsize)

    prdistance, prind1d = kdt.query(zip(cit_x2,cit_y2,cit_alt),k=query_k, eps=query_eps, p=query_p,
                            distance_upper_bound=query_distance_upper_bound,n_jobs=query_n_jobs)
    
    #
    
    
    #if query_k >1 means you are considering more than one gate and an average is needed
    
    if query_k > 1:
        
        #Issue with prind1d being the size of apr_ku...
        
        ind = np.where(prind1d == apr_ku.shape[0])
        if len(ind[0]) > 0 or len(ind[1]) > 0:
            print('gate was outside distance upper bound, eliminating those instances')
            #mask values outside search area
            prind1d[ind] = np.ma.masked
            prdistance[ind] = np.ma.masked
            
        #Barnes weighting     
        W_d_k = np.ma.array(np.exp(-1*prdistance**2./K_d**2.))
        W_d_k2 = np.ma.masked_where(np.ma.getmask(apr_ku[prind1d]), W_d_k.copy())
        w1 = np.ma.sum(W_d_k2 * 10. **(apr_ku[prind1d] / 10.),axis=1)
        w2 = np.ma.sum(W_d_k2, axis=1)
        ku_temp = 10. * np.ma.log10(w1/w2)
        
        W_d_k = np.ma.array(np.exp(-1*prdistance**2./K_d**2.))
        W_d_k2 = np.ma.masked_where(np.ma.getmask(apr_ka[prind1d]), W_d_k.copy())
        w1 = np.ma.sum(W_d_k2 * 10. **(apr_ka[prind1d] / 10.),axis=1)
        w2 = np.ma.sum(W_d_k2, axis=1)
        ka_temp = 10. * np.ma.log10(w1/w2)

        W_d_k2 = np.ma.masked_where(np.ma.getmask(apr_DFR[prind1d]), W_d_k.copy())
        w1 = np.ma.sum(W_d_k2 * 10. **((apr_DFR[prind1d]) / 10.),axis=1)
        w2 = np.ma.sum(W_d_k2, axis=1)
        dfr_temp = 10. * np.ma.log10(w1/w2)

        W_d_k = np.ma.array(np.exp(-1*prdistance**2./K_d**2.))
        W_d_k2 = np.ma.masked_where(np.ma.getmask(prdistance), W_d_k.copy())
        w1 = np.ma.sum(W_d_k2 * 10. **(prdistance/ 10.),axis=1)
        w2 = np.ma.sum(W_d_k2, axis=1)
        dis_temp = 10. * np.ma.log10(w1/w2)
    
        Ku_gate = ku_temp
        Ka_gate = ka_temp
        DFR_gate = dfr_temp

        #append current lat,lon and alt of the citation plane
        lat_c = np.append(lat_c,cit_lat)
        lon_c = np.append(lon_c,cit_lon)
        alt_c = np.append(alt_c,cit_alt)
        t_c = np.append(t_c,cit_time2)
        #

        #Use plane location for barnes averaged radar value
        lat_r = cit_lat
        lon_r = cit_lon
        alt_r = cit_alt
        t_r = cit_time2
        #
        dis_r = dis_temp
        ind_r = np.nan
        
        
        t_tiled = np.empty([t_c.shape[0],query_k],dtype=object)
        for i in np.arange(0,t_c.shape[0]):
            t_tiled[i,:] = t_c[i]
        diftime = apr_t[prind1d] - t_tiled
        diftime2 = np.empty(diftime.shape)
        for i in np.arange(0,diftime.shape[0]-1):
            for j in np.arange(0,diftime.shape[1]-1):
                diftime2[i,j] = diftime[i,j].total_seconds()
        
        W_d_k = np.ma.array(np.exp(-1*prdistance**2./K_d**2.))
        W_d_k2 = np.ma.masked_where(np.ma.getmask(diftime2), W_d_k.copy())
        w1 = np.ma.sum(W_d_k2 * 10. **(diftime2/ 10.),axis=1)
        w2 = np.ma.sum(W_d_k2, axis=1)
        dif_temp = 10. * np.ma.log10(w1/w2)
        
        dif_t = dif_temp
        
   
    else:
            
        #If gate outside sphere will need to remove flaged data == apr_ku.shape[0]
        ind = np.where(prind1d == apr_ku.shape[0])
        if len(ind[0]) > 0:
            print('gate was outside distance upper bound, eliminating those instances')
            #mask ind and distances that are outside the search area
            prind1d[ind] = np.ma.masked
            prdistance[ind] = np.ma.masked
               
        ku_temp = apr_ku[prind1d]
        ka_temp = apr_ka[prind1d]
        dfr_temp = ku_temp - ka_temp
        Ku_gate = np.append(Ku_gate,ku_temp)
        Ka_gate = np.append(Ka_gate,ka_temp)
        DFR_gate = np.append(DFR_gate,dfr_temp)
        #

        #append current lat,lon and alt of the citation plane
        lat_c = np.append(lat_c,cit_lat)
        lon_c = np.append(lon_c,cit_lon)
        alt_c = np.append(alt_c,cit_alt)
        t_c = np.append(t_c,cit_time2)
        conc_hvps3 = np.append(conc_hvps3,bigins)
        #

        #Get radar gate info and append it
        lat_r = np.append(lat_r,apr_y[prind1d])
        lon_r = np.append(lon_r,apr_x[prind1d])
        alt_r = np.append(alt_r,apr_alt[prind1d])
        t_r = np.append(t_r,apr_t[prind1d])
        dis_r = np.append(dis_r,prdistance)
        ind_r = np.append(ind_r,prind1d)

        dif_t = np.nan

    
    #Make lists
    matcher = {}
    Cit = {}
    APR = {}
    matched = {}
    kdtree = {}
    
    #Pack values in lists for export
    kdtree['prind1d'] = prind1d
    kdtree['prdistance'] = prdistance
    kdtree['query_k'] = query_k
    
    Cit['lat'] = cit_lat
    Cit['lon'] = cit_lon
    Cit['alt'] = cit_alt
    Cit['time'] = cit_time2
    APR['lat'] = apr_y
    APR['lon'] = apr_x
    APR['alt'] = apr_alt
    APR['Ku'] = apr_ku
    APR['Ka'] = apr_ka
    APR['DFR'] = apr_ku - apr_ka
    APR['time'] = apr_t

    matched['Ku'] = Ku_gate
    matched['Ka'] = Ka_gate
    matched['DFR'] = DFR_gate
    matched['lat_r'] = lat_r
    matched['lon_r'] = lon_r
    matched['alt_r'] = alt_r
    matched['lat_c'] = lat_c
    matched['lon_c'] = lon_c
    matched['alt_c'] = alt_c
    matched['time_r'] = t_r
    matched['time_c'] = t_c
    matched['dist'] = dis_r
    matched['dif_t'] = dif_t
    matched['array index'] = ind_r
    matched['conc_hvps3'] = conc_hvps3
    
    matcher['Cit'] = Cit
    matcher['APR'] = APR
    matcher['matched'] = matched
    matcher['kdtree'] = kdtree
    
    #Several plots to visualize data
    if plotson:
        fontsize=fontsize
        matched = matcher
        
        if query_k <= 1:
            diftime = matched['matched']['time_r'] - matched['matched']['time_c']
            diftime2 = np.array([])
            for i in np.arange(0,diftime.shape[0]):
                diftime2 = np.append(diftime2,diftime[i].total_seconds())
        else:
            diftime2= matched['matched']['dif_t']
            

        fig1,axes = plt.subplots(1,2,)
        
        #ax1 is the histogram of times
        ax1 = axes[0]
        ax1.hist(diftime2/60.,facecolor='b',edgecolor='k')
        ax1.set_xlabel('$t_{gate} - t_{Cit}, [min]$')
        ax1.set_ylabel('Number of gates')
        ax1.set_title(matched['matched']['time_r'][0])
        #ax2 is the histogram of distances
        ax2 = axes[1]
        distances = matched['matched']['dist']
        ax2.hist(distances,facecolor='r',edgecolor='k')
        ax2.set_xlabel('Distance, $[m]$')
        ax2.set_ylabel('Number of gates')
        ax2.set_title(matched['matched']['time_r'][0])

        plt.tight_layout()

        #Print some quick stats
        print(distances.shape[0],np.nanmean(diftime2)/60.,np.nanmean(distances))
        #
        
        fig = plt.figure()
        #ax3 is the swath plot to show radar and plane location
        ax3 = plt.gca()
        apr = apr3read(apr3filename)
        lat3d = apr['lat_gate']
        lon3d = apr['lon_gate']
        alt3d = apr['alt_gate']
        radar_n = apr['Ku']

        lon_s = np.empty(alt3d.shape[1:])
        lat_s = np.empty(alt3d.shape[1:])
        swath = np.empty(alt3d.shape[1:])
        for i in np.arange(0,alt3d.shape[2]):
            for j in np.arange(0,alt3d.shape[1]):
                ind = np.where(alt3d[:,j,i]/1000. > 3.5)
                ind2 = np.where(alt3d[:,j,i]/1000. < 3.6)
                ind3 = np.intersect1d(ind,ind2)
                ind3= ind3[0]
                l1 = lat3d[ind3,j,i]
                l2 = lon3d[ind3,j,i]
                k1 = radar_n[ind3,j,i]
                lon_s[j,i] = l2
                lat_s[j,i] = l1
                swath[j,i] = k1

        area_def = pr.geometry.AreaDefinition('areaD', 'IPHEx', 'areaD',
                                          {'a': '6378144.0', 'b': '6356759.0',
                                           'lat_0': '47.7998', 'lat_ts': '47.7998','lon_0': '-123.7066', 'proj': 'stere'},
                                          400, 400,
                                          [-70000., -70000.,
                                           70000., 70000.])
        bmap = pr.plot.area_def2basemap(area_def,resolution='l',ax=ax3)
        bmap.drawcoastlines(linewidth=2)
        bmap.drawstates(linewidth=2)
        bmap.drawcountries(linewidth=2)
        parallels = np.arange(-90.,90,4)
        bmap.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
        meridians = np.arange(180.,360.,4)
        bmap.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)
        bmap.drawmapboundary(fill_color='aqua')
        bmap.fillcontinents(lake_color='aqua')

        x,y = bmap(lon_s,lat_s)
        swath[np.where(swath < 0)] = np.nan
        pm1 = bmap.pcolormesh(x,y,swath,vmin=0,vmax=40,zorder=11,cmap='seismic')
        cbar1 = plt.colorbar(pm1,label='$Z_m, [dBZ]$')

        x2,y2 = bmap(matched['matched']['lon_c'],matched['matched']['lat_c'])
        pm2 = bmap.scatter(x2,y2,c=diftime2/60.,marker='o',zorder=12,cmap='PuOr',edgecolor=[],vmin=-10,vmax=10)
        cbar2 = plt.colorbar(pm2,label = '$\Delta{t}, [min]$')

        ax3.set_ylabel('Latitude',fontsize=fontsize,labelpad=20)
        ax3.set_xlabel('Longitude',fontsize=fontsize,labelpad=20)

        plt.tight_layout()
        plt.show()
        
        #Plot timeseries of barnes averaged or closest gate.
        plt.figure()
        plt.plot(matched['matched']['time_c'],matched['matched']['Ku'],'b',label='Ku')
        plt.plot(matched['matched']['time_c'],matched['matched']['Ka'],'r',label='Ka')
        plt.xlabel('Time')
        plt.ylabel('Z, [dBZ]')
        plt.legend()
       
    print('done')
    return matcher





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
    if 'zvv95' in listofkeys:
        radar3 = hdf.select(radar_freq3)
        radar_n3 = radar3.get()
        radar_n3 = radar_n3/100.
    else:
        radar_n3 = np.array([])
        flag = 1
        print('No W band')
    
    
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
    #alt = alt[:,scan,:]
    lat = lat.get()
    #lat = lat[scan,:]
    lon = lon.get()
    #lon = lon[scan,:]
    
    lat3d = lat3d.get()
    lat3d = (lat3d/lat3d_scale) + lat3d_offset
    lon3d = lon3d.get()
    lon3d = (lon3d/lon3d_scale) + lon3d_offset
    alt3d = alt3d.get()
    alt3d = (alt3d/alt3d_scale) + alt3d_offset
    
    #time = time[scan,:]
    #surf = surf[scan,:]
    #isurf = isurf[scan,:]
    #plane = plane[scan,:]
    radar_n = radar.get()
    radar_n = radar_n/100.
    radar_n2 = radar2.get()
    radar_n2 = radar_n2/100.

    radar_n4 = radar4.get()
    radar_n4 = radar_n4/100.
    vel_n = vel.get()
    vel_n = vel_n/100.
    
    #Quality control
    radar_n[radar_n <= -99] = np.nan
    radar_n2[radar_n2 <= -99] = np.nan
    radar_n4[radar_n4 <= -99] = np.nan
    
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
        apr['DFR_2'] = radar_n2 - radar_n3 #Ka - W
        apr['DFR_3'] = radar_n - radar_n3 #Ku - W
        apr['info'] = 'The shape of these arrays are: Radar[Vertical gates,Time/DistanceForward]'
    else:
        apr['DFR_2'] = np.array([]) #Ka - W
        apr['DFR_3'] = np.array([]) #Ku - W
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
