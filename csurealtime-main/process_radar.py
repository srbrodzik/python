# Brody Fuchs, Sept 2017
# brfuchs@atmos.colostate.edu

# Just some code to get started reading in radar data and doing some plotting and rainrate stuff
# This version will use the csuram stuff to make it easier

import numpy as np 
import matplotlib.pyplot as plt 
from netCDF4 import Dataset
from glob import glob
import matplotlib as mpl
from matplotlib.colors import LogNorm

from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, 
                            csu_kdp, csu_misc, csu_blended_rain_tropical)

from mpl_toolkits.basemap import Basemap 

import argparse
from copy import deepcopy
import datetime
import time
import raintype
from csuram import RadarData, Cell, general_tools, Case, RadarConfig
from skewt import SkewT
import argparse


#parser = argparse.ArgumentParser(
#    description='Unfold radar velocities for a day and write the results in CF/Radial format.')
#parser.add_argument('fname', type=str, help='The name of the file')
#parser.add_argument('startfile', type=str, help='The file to get the start sounding volume from')
#args = parser.parse_args()


start = time.time()

parser = argparse.ArgumentParser(description='Put in a file to be processed')

#parser.add_argument('--noarg', action="store_true", default=False)
parser.add_argument('--realtime', action="store", dest="realtime", type=bool, default=False)



cs_colors =['#FFFFFF', 'DodgerBlue', 'Red', 'Khaki']
cs_cmap = mpl.colors.ListedColormap(cs_colors)
cs_bounds = np.arange(0,5)
cs_norm = mpl.colors.BoundaryNorm(cs_bounds, cs_cmap.N)


map_limits = {'lat': [39.0, 41.5], 'lon': [-106.0, -103]}
map_res = 'c'


radar_path = '/Users/brodyfuchs-macbook/code/rainfallmap/realtime_gridded'

# find all radar files in the path
radar_files = glob('%s/*.nc'%(radar_path))


sounding_file = '/Users/brodyfuchs-macbook/code/rainfallmap/sample_tropical_sounding.txt'
# ******** Need to read this file and get it into the HID calculation **********

snd = SkewT.Sounding(sounding_file)


# defining the accumulated rainfall array, or list rather
accum_rainfall = []
accum_times = []

accum_minute_thresh1 = 30.
accum_minute_thresh2 = 60.
rainfall_age_thresh = 5.0 # 10 mm/hr
dtime = 4.0 # minute

cb_pad = 0.05
cb_frac = 0.046



xl = [-105.0, -103.5]
yl = [39.5, 40.5]


#fraction=0.046, pad=0.04

# defining the C/S algorithm parameters
minZdiff = 20 
deepcoszero = 40
shallowconvmin = 28
truncZconvthres = 43;
dBZformaxconvradius = 46;
# truncZconvthres = self.zconv
# dBZformaxconvradius = self.zconv+3;
weakechothres = 7
backgrndradius = 5       #(in km)
maxConvRadius = 10       #(in km)
minsize = 8              #(in km^2)
startslope = 50          #(in km^2)
maxsize = 2000           #(in km^2)


for irf, rf in enumerate(radar_files):
    # read in the file
#    print 'starting: {}'.format(time.time() - start)

    #data = Dataset(rf, 'r')
    radar = RadarData.RadarData(radar_file=rf, x='x0', y='y0', z='z0', band='C', squeeze=True, 
                    lat='lat0', lon='lon0', dz='DZQC', kdp='KDP', zdr='DR', rho='RH')


    radar.add_sounding_object(snd)
    radar.interp_sounding()
    radar.set_hid(band='C', use_temp=True) # this calls the appropriate functions and uses the sounding that we've
        # attributed to the radar object and uses the proper radar band


    last_dot = rf.rfind('.')

    dt_string = rf[last_dot-15:last_dot]
    print dt_string
    file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
    radar.date = file_time

    accum_times.append(file_time)

    # need to run clutter filter on the data using CSU radartools

    # rainfall algorithms

    # Now do some plotting
    # first plot the reflecitivity?
    lat = radar.data['lat0'][:]
    lon = radar.data['lon0'][:]
    z = radar.data['z0'][:]



    z_val = 2.0
    z_ind = np.argmin(np.abs(z-z_val))

    # Back out location of the radar using the x/y/lat/lon info
    x_0 = np.where(radar.data[radar.x_name] == 0)[0][0]
    y_0 = np.where(radar.data[radar.y_name] == 0)[0][0]

    radar_lat = radar.data[radar.lat_name][y_0, x_0]
    radar_lon = radar.data[radar.lon_name][y_0, x_0]





    dbz = radar.data['DZQC'][:].squeeze()[z_ind, :, :]
    zdr = radar.data['DR'][:].squeeze()[z_ind, :, :]
    kdp = radar.data['KDP'][:].squeeze()[z_ind, :, :]
    rho = radar.data['RH'][:].squeeze()[z_ind, :, :]


    #radar.set_hid(use_temp=False)


    # scores = csu_fhc.csu_fhc_summer(dz=dbz, zdr=zdr, rho=rho, kdp=kdp, use_temp=True, band='C')
    # fh = np.argmax(scores, axis=0) + 1


    ##### Okay, now try to do the C/S partitioning


    cs_raintype, cs_rtypes = raintype.raintype(dbz, refl_missing_val=dbz.min(), 
                                   refl_dx=1.0, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                   shallowconvmin=shallowconvmin, truncZconvthres=truncZconvthres,
                                   dBZformaxconvradius=dBZformaxconvradius,
                                   weakechothres=weakechothres, backgrndradius=backgrndradius,
                                   maxConvRadius=maxConvRadius, minsize=minsize,
                                   startslope=startslope, maxsize=maxsize)



    # ***** Probably need to convert his categories into ours, maybe here before we plot??
    # ISO_CONV_CORE: 4, 'CONVECTIVE': 2, ISO_CS_CORE: 9 ARE CONVECTIVE WHICH WILL BE ASSINGED A 1
    # WEAK_ECHO: 6, STRATIFORM: 1, ISO_CONV_FRINGE: 5 ARE STRATIFORM AND GET A -1
    # UNCERTAIN: 3 IS mixed and will get a 0

    # first doing the convection
    conv = (cs_raintype == 4) | (cs_raintype == 2) | (cs_raintype == 9)
    strat = (cs_raintype == 6) | (cs_raintype == 1) | (cs_raintype == 5)
    mixed = cs_raintype == 3


    cs_arr = deepcopy(cs_raintype)
    cs_arr[conv] = 2
    cs_arr[strat] = 1
    cs_arr[mixed] = 3

    cs_arr[dbz<0] = 0

    radar.add_field(cs_arr, 'CS')


    trop_rr, trop_rr_meth = csu_blended_rain_tropical.calc_blended_rain_tropical(dz=dbz, zdr=zdr, kdp=kdp, 
                        cs=cs_arr, fhc=radar.data['HID'][z_ind])


    radar.add_field(trop_rr, 'RR')
    radar.add_field(trop_rr_meth, 'Rmethod')
    # Need to add in the rain rate as a field to the RadarData object, also need to add the convective/stratiform

    accum_rainfall.append(trop_rr.filled(0.0)*dtime/60.0)


    # okay now do accumulated rainfall for the last 30 minutes
    tdiffs = np.array([np.abs((_ - file_time).total_seconds()) for _ in accum_times])
    #print 'tdiffs: {}'.format(tdiffs)

    good_files1 = np.where(tdiffs <= accum_minute_thresh1*60.0)
    good_files2 = np.where(tdiffs <= accum_minute_thresh2*60.0)

    #print 'good files: {}'.format(good_files)

    accum_rainfall_good_times1 = np.sum(np.array(accum_rainfall)[good_files1], axis=0)
    accum_rainfall_good_times2 = np.sum(np.array(accum_rainfall)[good_files2], axis=0)

    radar.add_field(accum_rainfall_good_times1, 'ACCRR30')
    radar.add_field(accum_rainfall_good_times2, 'ACCRR60')



    # ***** Alright, now how do I do the rainfall age thing??

    rainfall_age = np.zeros_like(trop_rr)
    rainfall_age = np.ma.masked_where(rainfall_age == 0.0, rainfall_age)

    # set to 0's here
    accrain_array = np.array(accum_rainfall)
    enough_rain = np.where(accrain_array >= rainfall_age_thresh)

    # now need to loop thru each entry in enough_rain
    xy_pairs = zip(enough_rain[1], enough_rain[2])

    uniq_xy_pairs = list(set(xy_pairs))



    for pt in uniq_xy_pairs:
        # pt is an x,y pair
        all_pts = np.where( (enough_rain[1] == pt[0]) & (enough_rain[2] == pt[1]) )
        t_inds_thresh = enough_rain[0][all_pts[0]]
        last_t = t_inds_thresh[-1]
        this_age = (len(accum_rainfall) - last_t)*dtime - (dtime*0.9)
        rainfall_age[pt] = this_age





    radar.add_field(rainfall_age, 'RRAGE')



    # ****** NOW PLOTTING STUFF **********

    fig, ax = radar.cappi_multiplot(z=2.0, varlist=['DZQC', 'DR', 'KDP', 'RH', 'CS', 'HID'], coords='ll', xlim=xl, ylim=yl)


    for a in ax.flatten():
        a.scatter(radar_lon, radar_lat, s=50, facecolor='MediumPurple')
        a.grid(True)

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)

    plt.savefig('figures/%s_variables.png'%(dt_string))
#    print 'figure saved: {}'.format(time.time() - start)

    plt.close(fig)


    #quit()

#     #rr_hidro, rr_hidro_method = csu_blended_rain.csu_hidro_rain(dz=dz, zdr=dr, kdp=kdp, fhc=fh)




    fig, ax = plt.subplots(2,2, figsize=(10,8))
    axf = ax.flatten()

    rr_pc = radar.cappi('RR', z=z_val, ax=axf[0], cmap=plt.cm.BuPu, norm=LogNorm(vmin=1e-1, vmax=100), coords='ll', xlim=xl, ylim=yl)
    axf[0].set_title('Instantaneous rain rate (mm/hr)')

    arr_pc1 = radar.cappi('ACCRR30', z=z_val, ax=axf[1], cmap=plt.cm.afmhot_r, norm=LogNorm(vmin=1e-1, vmax=100), coords='ll', xlim=xl, ylim=yl)

#    arr_cb1 = plt.colorbar(arr_pc1, orientation='vertical', ax=axf[1], fraction=cb_frac, pad=cb_pad)
    axf[1].set_title('30 min accum rainfall (mm)')

    arr_pc2 = radar.cappi('ACCRR60', z=z_val, ax=axf[2], cmap=plt.cm.afmhot_r, norm=LogNorm(vmin=1e-1, vmax=100), coords='ll', xlim=xl, ylim=yl)

    # arr_pc2 = axf[2].pcolormesh(radar.data['lon0'], radar.data['lat0'], accum_rainfall_good_times2, 
    #                 cmap=plt.cm.afmhot_r, norm=LogNorm(vmin=1e-1, vmax=100))
#    arr_cb2 = plt.colorbar(arr_pc2, orientation='vertical', ax=axf[2], fraction=cb_frac, pad=cb_pad)
    axf[2].set_title('60 min accum rainfall (mm)')

    rrage_pc = radar.cappi('RRAGE', z=z_val, ax=axf[3], cmap=plt.cm.rainbow, vmin=0, vmax=120, coords='ll', xlim=xl, ylim=yl)

#    rrage_cb3 = plt.colorbar(rrage_pc, orientation='vertical', ax=axf[3], fraction=cb_frac, pad=cb_pad)
    axf[3].set_title('Rainfall age using %d mm threshold'%(rainfall_age_thresh))

    for a in axf:
        a.grid(True)
        a.scatter(radar_lon, radar_lat, s=50, facecolor='MediumPurple')


    plt.tight_layout()



    plt.savefig('rain_figs/%s_rraccum.png'%(dt_string))

    plt.close(fig)

















