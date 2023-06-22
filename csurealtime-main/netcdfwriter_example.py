import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
import glob
import os
from skewt import SkewT
from csu_blended_rain_olympex import calc_blended_rain_olympex
from csu_blended_rain_julie import calc_blended_rain

from csu_radartools import csu_blended_rain
from csu_radartools import csu_fhc
#from uw_raintype import raintype
from datetime import datetime
import cPickle as pickle

from sound_interp import check_sounding_for_montonic
from sound_interp import interpolate_sounding_to_grid


def write_output(infile_nc, outfile_nc, dza, dra, kda, rha,fhta,rra,metha):#,nwa,d0a):
    infile = Dataset(infile_nc, 'r')
    lat0 = infile.variables['lat0']
    lon0 = infile.variables['lon0']
    z0 = infile.variables['z0']
    x0 = infile.variables['x0']
    y0 = infile.variables['y0']
    time0 = infile.variables['time']
    MISSING_VALUE = -9999.
    nx = len(x0)
    ny = len(y0)
    nz = len(z0)
    mygroup = Dataset(outfile_nc, 'w', format='NETCDF3_64BIT')

    mygroup.createDimension('time', 0)
    mygroup.createDimension('x0', nx)
    mygroup.createDimension('y0', ny)
    mygroup.createDimension('z0', nz)

    time = mygroup.createVariable('time', 'f8', ('time',))
    y = mygroup.createVariable('y0', 'f4', ('y0',))
    x = mygroup.createVariable('x0', 'f4', ('x0',))
    z = mygroup.createVariable('z0', 'f4', ('z0',))

    lon = mygroup.createVariable('lon0', 'f4', ('y0','x0'), fill_value=MISSING_VALUE)
    lat = mygroup.createVariable('lat0', 'f4', ('y0','x0'), fill_value=MISSING_VALUE)
    dz = mygroup.createVariable('CZ', 'f4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
#     d0 = mygroup.createVariable('D0', 'f4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
#     nw = mygroup.createVariable('NW', 'f4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
    dr = mygroup.createVariable('DR', 'f4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
    kd = mygroup.createVariable('KD', 'f4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
    rh = mygroup.createVariable('RH', 'f4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
    rr = mygroup.createVariable('RR', 'f4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
    fht = mygroup.createVariable('FHT', 'i4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)
    meth = mygroup.createVariable('METH', 'i4', ('time','z0', 'y0', 'x0'), fill_value=MISSING_VALUE)

    time[:] = time0
    lon[:] = lon0
    lat[:] = lat0
    z[:] = z0
    x[:] = x0
    y[:] = y0
    
    dza=dza[np.newaxis,...]
#     d0a=d0a[np.newaxis,...]
#     nwa=nwa[np.newaxis,...]
    dra=dra[np.newaxis,...]
    kda=kda[np.newaxis,...]
    rra=rra[np.newaxis,...]
    rha=rha[np.newaxis,...]
    fhta=fhta[np.newaxis,...]
    metha=metha[np.newaxis,...]
    
    dz[...] = dza[...]
#     d0[...] = d0a[...]
#     nw[...] = nwa[...]
    dr[...] = dra[...]
    kd[...] = kda[...]
    rr[...] = rra[...]
    rh[...] = rha[...]
    fht[...] = fhta[...]
    meth[...] = metha[...]

    write_standard_attributes(time, lat, lon, x,y,z, dz, dr, kd, rh,fht,rr,meth,MISSING_VALUE)
    mygroup.close()
    infile.close()


def write_standard_attributes(time, lat, lon, x,y,z, dz, dr, kd, rh,fht,rr,meth,MISSING_VALUE):
    # attributes
    time.standard_name = "time"
    time.units = "time_elapsed"
    time.calendar = "standard"

    lat.standard_name = "latitude"
    lat.long_name = "latitude"
    lat.units = "degrees_north"
    lat.axis = "Y"

    lon.standard_name = "longitude"
    lon.long_name = "longitude"
    lon.units = "degrees_east"
    lon.axis = "X"

    z.standard_name = "altitude"
    z.long_name = "constant altitude levels"
    z.units = "km"
    z.positive="up"
    
    dz.long_name = "Reflectivity"
    dz.units = "dBZ"
    dz.missing_value = MISSING_VALUE

    dr.long_name = "Differential Reflectivity"
    dr.units = "dB"
    dr.missing_value = MISSING_VALUE

    kd.long_name = "Specific differential reflectivity"
    kd.units = "deg/km"
    kd.missing_value = MISSING_VALUE

    rh.long_name = "correlation H-to_v"
    rh.units = " "
    rh.missing_value = MISSING_VALUE


    fht.long_name = "Hydrometer Identification"
    fht.units = " "
    fht.missing_value = MISSING_VALUE

    rr.long_name = "Rain Rate CSU-OLYMPEX"
    rr.units = "mm /hr "
    rr.missing_value = MISSING_VALUE

    meth.long_name = "Rain Rate CSU-OLYMPEX Method"
    meth.units = "unitless "
    meth.missing_value = MISSING_VALUE

#     nw.long_name = "Normalized Intercept Parameter NASA"
#     nw.units = "log(Nw) unitless "
#     nw.missing_value = MISSING_VALUE
# 
#     d0.long_name = "Median Drop Diameter NASA"
#     d0.units = "mm "
#     d0.missing_value = MISSING_VALUE


def check_sounding_for_montonic(sounding):
    """
    So the sounding interpolation doesn't fail, force the sounding to behave
    monotonically so that z always increases. This eliminates data from
    descending balloons.
    """
    snd_T = sounding.soundingdata['temp']  # In old SkewT, was sounding.data
    snd_z = sounding.soundingdata['hght']  # In old SkewT, was sounding.data
    dummy_z = []
    dummy_T = []
    if not snd_T.mask[0]: #May cause issue for specific soundings
        dummy_z.append(snd_z[0])
        dummy_T.append(snd_T[0])
        for i, height in enumerate(snd_z):
            if i > 0:
                if snd_z[i] > snd_z[i-1] and not snd_T.mask[i]:
                    dummy_z.append(snd_z[i])
                    dummy_T.append(snd_T[i])
        snd_z = np.array(dummy_z)
        snd_T = np.array(dummy_T)
    return snd_T, snd_z


def get_freezing_altitude(T, z):
    """
    T = temperature (deg C), z = altitude (any units)
    Assumes T & z are aligned, unmasked ndarrays of same shape.
    Returns highest possible freezing altitude.
    Partial soundings (i.e., no pass thru freezing altitude) get fitted to
    a regression line and Hfrz is then estimated via the intercept.
    """
    if '0.0' in T:
    # Maybe there's already one or more zero points in sounding
        return np.max(z[T == 0])
    else:
        sin = np.argsort(z)
        zs = z[sin]
        Tz = T[sin]
        def find_nearest(a, a0):
            "Element in nd array `a` closest to the scalar value `a0`"
            idx = np.abs(a - a0).argmin()
            return idx
        T_nearest= find_nearest(Tz,0)
        print Tz[T_nearest]
        if Tz[T_nearest] < 0:
            T2 = Tz[T_nearest-1:T_nearest+1]
            z2 = zs[T_nearest-1:T_nearest+1]
#            print T2,z2
            interpolated = z2[0]+((0.0-T2[1])/(T2[0]-T2[1]))*(z2[1]-z2[0])
            return interpolated
        else:
            T2 = Tz[T_nearest:T_nearest+2]
            z2 = zs[T_nearest:T_nearest+2]
#            print T2,z2
            interpolated = z2[0]+((0.0-T2[1])/(T2[0]-T2[1]))*(z2[1]-z2[0])
            return interpolated            


###START the program


#dates=['20151114','20151117','20151203','20151211','20151214','20151219','20151222','20160123','20160130']
dates=['20151211','20151219','20160111']
#dates=['20151219']

for date1 in dates:
    #Read in the sounding data
    sndfile = (glob.glob('{d}_*.txt'.format(d=date1)))[0]

    radar ='KLGX'
#    fdir= '/Volumes/andes/data/olympex/OLYMPEX_NPOL/OVERPASS/1km/{d}/'.format(d=date1)
    fdir= '/Volumes/andes/data/olympex/OLYMPEX_NPOL/KLGX/'.format(d=date1)

    odir= '/Volumes/andes/data/olympex/OLYMPEX_NPOL/KLGX/'
#    codir= '/Users/bdolan/scratch/OLYMPEX_NPOL/'
    codir= '/Users/bdolan/scratch/OLYMPEX_NPOL/'

    cdates=[]
    bflist=[]
    cflist = sorted(glob.glob('{f}KLGX_ppi_cart2_{d:}*.nc'.format(f=fdir,d=date1)))
    for v,cname in enumerate(cflist):
        base = os.path.basename(cname)
        temp=base[15:30]
        print temp
        cdates.append(datetime.strptime(temp,'%Y%m%d_%H%M%S'))



        #sndfile = '{d}_*.txt'.format(d=date1)
        print sndfile
        sounding = SkewT.Sounding(sndfile)

        cfile = Dataset(cname)

        dz = np.squeeze(cfile.variables['CZ'][:])
        dr = np.squeeze(cfile.variables['DR'][:])
        rh = np.squeeze(cfile.variables['RHO'][:])
        kd = np.squeeze(cfile.variables['KD'][:])
     #    rr = np.squeeze(cfile.variables['RR'][:])
# #         nw = np.squeeze(cfile.variables['NW'][:])
#         d0 = np.squeeze(cfile.variables['D0'][:])

        radx = cfile.variables['x0'][:]
        rady = cfile.variables['y0'][:]
        radz = cfile.variables['z0'][:]

        TT,ZZ = check_sounding_for_montonic(sounding)
        freezing_altitude = get_freezing_altitude(TT,ZZ)
        freezing_altitude = freezing_altitude/np.float(1000)
        #print freezing_altitude


        # Get everything in the same array size with interpolate_sounding_to_grid and make sure it's a 4-D array
        dum_dbz=dz
        z = np.array(radz)
    #    print np.shape(z)
    #    print np.shape(dum_dbz)
        radar_T, radar_z = interpolate_sounding_to_grid(sounding,z,dum_dbz)
        new_T = np.zeros(dz.shape)
        new_z = np.zeros(dz.shape)
        dbzzz =  dz.shape[0]
    #    print np.shape(new_T)
    
        CSV_DIR='/Users/bdolan/scratch/csu_radartools/beta_function_parameters/'

        scorecs = csu_fhc.csu_fhc_summer(dz=dz, zdr=dr, rho=rh, kdp=kd, use_temp=True, band='S',
                                        T=radar_T)
        fhc = np.argmax(scorecs, axis=0) + 1
        fhc = np.ma.masked_where(dz.mask,fhc)

        rr,meth=calc_blended_rain_olympex(dz=dz,zdr=dr,kdp=kd,band='S')


        outfile_nc = '{o}KLGX_OLY_{d}.cdf'.format(o=odir,d=temp)
        write_output(cname, outfile_nc, dz, dr, kd, rh,fhc,rr,meth)#,nw,d0)

        print 'done writing Netcdf! Moving on...'

        fig, ax = plt.subplots(3,2,figsize=(16,16))
        axf = ax.flatten()
        level =0
        label = ['Drizzle', 'Rain', 'Ice Crystals', 'Aggregates',
                                   'Wet Snow', 'Vertical Ice', 'LD Graupel',
                                   'HD Graupel', 'Hail', 'Big Drops']
        hid_colors = ['White', 'LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                      'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
        rr_colors = ['White', 'LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
                      'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']

        cmaphid = colors.ListedColormap(hid_colors)

        boundshid=np.arange(0,np.max(np.shape(hid_colors))+1,1)
        normhid = colors.BoundaryNorm(boundshid, cmaphid.N)

        #print np.shape(hiddiff)
        hidcmap = colors.ListedColormap(hid_colors)
        dbz_hand=axf[0].contourf(radx,rady,dz[level,:,:],cmap='jet', levels=np.arange(-10,60,1))
        cb = fig.colorbar(dbz_hand, ax = axf[0],cmap='jet')
        axf[0].set_title('Reflectivity')

        hid_hand =axf[1].imshow(fhc[level,:,:],cmap=cmaphid,origin='lower',aspect='auto',interpolation='none',norm=normhid)
        cb = fig.colorbar(hid_hand, ax = axf[1],cmap=cmaphid)
        cb.set_ticks(boundshid)
        labs = np.concatenate((np.array(['']), np.array(label)))
        cb.set_ticklabels(labs)
        axf[1].set_title('HID')


        zdr_hand=axf[2].contourf(radx,rady,dr[level,:,:],cmap='jet', levels=np.arange(-2,6,0.1))
        cb = fig.colorbar(zdr_hand, ax = axf[2],cmap='jet')
        axf[2].set_title('ZDR')


        kdp_hand=axf[3].contourf(radx,rady,kd[level,:,:],cmap='jet', levels=np.arange(-1,3,0.1))
        cb = fig.colorbar(kdp_hand, ax = axf[3],cmap='jet')
        axf[3].set_title('Kdp')


        rho_hand=axf[4].contourf(radx,rady,rh[level,:,:],cmap='jet', levels=np.arange(0.8,1.01,0.01))
        cb = fig.colorbar(rho_hand, ax = axf[4],cmap='jet')
        axf[4].set_title('RHO')

        rr_hand=axf[5].pcolormesh(radx,rady,rr[level,:,:],cmap='jet', norm = colors.LogNorm(vmin=0.1,vmax=40.))
        #rr_hand = axf[5].pcolormesh(radx, rady, rr[level,:,:], vmin=0.1, vmax=40., cmap='jet')
        cb = fig.colorbar(rr_hand, ax = axf[5])
        cb.set_ticks([0.1,0.5,1.0,5.0,10.,25.0,50.])
        cb.set_ticklabels([0.1,0.5,1.0,5.0,10.,25.0,50.])
        axf[5].set_title('RR')

        plt.tight_layout()

        plt.savefig('OLYMPEX_KLGX_{c:%Y%m%d_%H%M%S}_6panel.png'.format(c=cdates[v],r=radar),dpi=200)

        radar_rain_all = rr
        radar_method_all = meth
        vol_time = 15.
        radar = 'KLGX'

        fig, ax = plt.subplots(1,2,figsize=(16,16))
        axf = ax.flatten()


        kdsub = kd[kd > 0.]
        print type(kdsub)
        dat=axf[0].hist(np.ma.compressed(kdsub),bins=np.arange(0,1.5,0.02),align='left',normed=1)

        zdsub = dr[dr > 0.]
        print type(zdsub)
        dat=axf[1].hist(np.ma.compressed(zdsub),bins=np.arange(0,2.0,0.02),align='left',normed=1)

        plt.savefig('KdpZdrhist_olympex_KLGX_{c:%Y%m%d_%H%M%S}.png'.format(c=cdates[v]),dpi=200)

        #Code for plotting rain volume and methods: (where raing and methg came from the call to csu_blended_rain)
        # rainvg=np.zeros(np.shape(dis_rain))
        rainvg=np.zeros(np.shape(radar_rain_all))
        # rainvg[dis_rain > 0]=dis_rain[dis_rain>0]*10./60. #Sum up accumulation over 1 minute for the Disdrometer, change to volume time for radar accumulations
        rainvg[radar_rain_all > 0]=radar_rain_all[radar_rain_all>0]*vol_time

        #Calculate the contributions to the total rain volume from each method.
        tot_rainv=np.ma.sum(rainvg)
        rainvm1=np.ma.sum(rainvg[radar_method_all==1])/tot_rainv*100.
        rainvm2=np.ma.sum(rainvg[radar_method_all==2])/tot_rainv*100.
        rainvm3=np.ma.sum(rainvg[radar_method_all==3])/tot_rainv*100.
        rainvm4=np.ma.sum(rainvg[radar_method_all==4])/tot_rainv*100.
        #rainvm5=np.sum(rainvg[radar_method==5])/tot_rainv*100.
        #rainvm6=np.sum(rainvg[radar_method==6])/tot_rainv*100.
        #ind=np.arange(1,7,1)
        ind = np.arange(1,5,1)
        rainvolmethod=[rainvm1,rainvm2,rainvm3,rainvm4]

        plt.bar(ind,rainvolmethod)
        x = np.arange(1.5,7.5,1)
        plt.title('Rain volume by method',y=1.04)
        # plt.title('Disdrometer rain volume by method {d}'.format(d=date),y=1.04)
        plt.xticks(x)
        plt.ylabel('Contribution to rain volume ')
        plt.gca().set_xticklabels(['R(Kdp,Zdr)', 'R(Kdp)', 'R(Z,Zdr)', 'R(Z)'])
        # plt.show()
        plt.savefig('{r}_rainvol_bymeth_KLGX_{c:%Y%m%d_%H%M%S}.png'.format(c=cdates[v],r=radar),dpi=200,bbox_inches='tight')
        # plt.savefig('disdrometer_rainvol_bymeth_{d}_end.png'.format(d=date),dpi=200,bbox_inches='tight')


        dat=plt.hist(np.ma.compressed(radar_method_all[:,:,:]),bins=[1,2,3,4,5],normed=1,align='left')
        # dat=plt.hist(np.ma.compressed(dis_method[:]),bins=[1,2,3,4,5],normed=1,align='left')
        plt.xticks()
        x = np.arange(1,8,1)
        plt.xticks(x)
        plt.gca().set_xticklabels(['R(Kdp,Zdr)', 'R(Kdp)', 'R(Z,Zdr)', 'R(Z)'])
        plt.ylabel('Frequency of occurrence')
        plt.title("{r} Blended Rain Estimator frequency".format(r=radar),y=1.04)
        # plt.title("Disdrometer Blended Rain Estimator frequency {d}".format(d=date),y=1.04)
        # plt.show()
        plt.savefig('{r}_rainrelationfreq_KLGX_{c:%Y%m%d_%H%M%S}.png'.format(c=cdates[v],r=radar),dpi=200,bbox_inches='tight')
        # plt.savefig('disdrometer_rainrelationfreq_{d}_end.png'.format(d=date),dpi=200,bbox_inches='tight')
        print 'dz:',np.shape(dz)

        dz1=dz[...,np.newaxis]
        dr1=dr[...,np.newaxis]
        rr1=rr[...,np.newaxis]
        meth1=meth[...,np.newaxis]
        kdp1=kd[...,np.newaxis]
        hid1=fhc[...,np.newaxis]
        rh1=rh[...,np.newaxis]
#         d01=d0[...,np.newaxis]
#         nw1=nw[...,np.newaxis]


        if v == 0:
            dz4d=dz[...,np.newaxis]
            dr4d=dr[...,np.newaxis]
            rr4d=rr[...,np.newaxis]
            meth4d=meth[...,np.newaxis]
            kdp4d=kd[...,np.newaxis]
            hid4d=fhc[...,np.newaxis]
            rh4d=rh[...,np.newaxis]
#             d04d=d0[...,np.newaxis]
#             nw4d=nw[...,np.newaxis]



    
        else:
            dz4d=np.ma.concatenate([dz4d,dz1],axis=3)
            dr4d=np.ma.concatenate([dr4d,dr1],axis=3)
            rr4d=np.ma.concatenate([rr4d,rr1],axis=3)
            kdp4d=np.ma.concatenate([kdp4d,kdp1],axis=3)
            hid4d=np.ma.concatenate([hid4d,hid1],axis=3)
            rh4d=np.ma.concatenate([rh4d,rh1],axis=3)
#             d04d=np.ma.concatenate([d04d,d01],axis=3)
#             nw4d=np.ma.concatenate([nw4d,nw1],axis=3)
            meth4d=np.ma.concatenate([meth4d,meth1],axis=3)
        
        print 'dz4d:',np.shape(dz4d)
        
    mydata={'dbz':dz4d,'hid':hid4d,'kdp':kdp4d,'rr':rr4d,'meth':meth4d,'zdr':dr4d,
        'Time':cdates,'rho':rh4d}	

    pickle.dump(mydata,open('{o}{d}_data_f_overpass_KLGX.p'.format(o=odir,d=date1),'wb'))
