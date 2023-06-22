# Brody Fuchs, Sept 2017
# brfuchs@atmos.colostate.edu

# Just some code to get started reading in radar data and doing some plotting and rainrate stuff

import numpy as np 
import matplotlib.pyplot as plt 
from netCDF4 import Dataset
from glob import glob
import matplotlib as mpl

from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, 
                            csu_kdp, csu_misc, csu_blended_rain_tropical)

from mpl_toolkits.basemap import Basemap 

from copy import deepcopy
import datetime
import time
import raintype
from csuram import RadarData, Cell, general_tools, Case, RadarConfig


start = time.time()


cs_colors =['#FFFFFF', 'DodgerBlue', 'Red', 'Khaki']
cs_cmap = mpl.colors.ListedColormap(cs_colors)
cs_bounds = np.arange(0,5)
cs_norm = mpl.colors.BoundaryNorm(cs_bounds, cs_cmap.N)


map_limits = {'lat': [39.0, 41.5], 'lon': [-106.0, -103]}
map_res = 'c'


radar_path = '/Users/brodyfuchs-macbook/code/rainfallmap/output_chil/20170901'

# find all radar files in the path
radar_files = glob('%s/*.nc'%(radar_path))

# defining the accumulated rainfall array, or list rather
accum_rainfall = []
accum_times = []

accum_minute_thresh = 30.
dtime = 4.0 # minute

cb_pad = 0.1
cb_frac = 0.046

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


for rf in radar_files:
    # read in the file
#    print 'starting: {}'.format(time.time() - start)

    data = Dataset(rf, 'r')
    last_dot = rf.rfind('.')

    dt_string = rf[last_dot-15:last_dot]
    print dt_string
    file_time = datetime.datetime.strptime(dt_string, '%Y%m%d_%H%M%S')
    accum_times.append(file_time)

    # need to run clutter filter on the data using CSU radartools

    # rainfall algorithms

    # Now do some plotting

    # first plot the reflecitivity?
    lat = data.variables['lat0'][:]
    lon = data.variables['lon0'][:]
    z = data.variables['z0'][:]

    z_val = 2.0
    z_ind = np.argmin(np.abs(z-z_val))


    dbz = data.variables['DZ'][:].squeeze()[z_ind, :, :]
    zdr = data.variables['DR'][:].squeeze()[z_ind, :, :]
    kdp = data.variables['KD'][:].squeeze()[z_ind, :, :]
    rho = data.variables['RH'][:].squeeze()[z_ind, :, :]


    scores = csu_fhc.csu_fhc_summer(dz=dbz, zdr=zdr, rho=rho, kdp=kdp, use_temp=True, band='C')
    fh = np.argmax(scores, axis=0) + 1



    fig, ax = plt.subplots(3, 2, figsize=(9,11))
    axf = ax.flatten()


    dbz_pc = axf[0].pcolormesh(lon, lat, dbz, vmin=0, vmax=70)
    dbz_cb = plt.colorbar(dbz_pc, orientation='horizontal', ax=axf[0], fraction=cb_frac, pad=cb_pad)
    axf[0].set_title('Reflectivity (dBZ)')

    # now try the despeckle thing?
    mask_ds = csu_misc.despeckle(dbz, ngates=4)

    # Not sure how to do the radar QC yet, the CSU radartools tutorial on Tim Lang's git page
    # use standard deviation of phase, which I can't do right now. Just gonna skip it for now.

#    print 'Rainfall algorithm: {}'.format(time.time() - start)

    rr_blend, rr_blend_method = csu_blended_rain.calc_blended_rain(dz=dbz, zdr=zdr, kdp=kdp)

    accum_rainfall.append(rr_blend.filled(0.0)*dtime/60.0)

    #rr_hidro, rr_hidro_method = csu_blended_rain.csu_hidro_rain(dz=dz, zdr=dr, kdp=kdp, fhc=fh)



    zdr_pc = axf[1].pcolormesh(lon, lat, zdr, cmap=plt.cm.PRGn, vmin=-4.0, vmax=4.0)
    zdr_cb = plt.colorbar(zdr_pc, orientation='horizontal', ax=axf[1], fraction=cb_frac, pad=cb_pad)
    axf[1].set_title('ZDR')

    kdp_pc = axf[2].pcolormesh(lon, lat, kdp, cmap=plt.cm.PuRd, vmin=-1.0, vmax=5.0)
    kdp_cb = plt.colorbar(kdp_pc, orientation='horizontal', ax=axf[2], fraction=cb_frac, pad=cb_pad)
    axf[2].set_title('Kdp')

    rr_pc = axf[3].pcolormesh(lon, lat, rr_blend, cmap=plt.cm.afmhot_r, vmin=1.0, vmax=50)
    rr_cb = plt.colorbar(rr_pc, orientation='horizontal', ax=axf[3], fraction=cb_frac, pad=cb_pad)
    axf[3].set_title('Rain rate (mm/hr)')


#    print 'Plotting fields: {}'.format(time.time() - start)




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





    # m4 = Basemap(projection='cyl', resolution=map_res, llcrnrlon=map_limits['lon'][0],  urcrnrlon=map_limits['lon'][1],
    #                   llcrnrlat=map_limits['lat'][0], urcrnrlat=map_limits['lat'][1], ax=axf[4])

    # m4.drawstates()
    # m4.drawcounties()
    # m4.drawrivers(color='blue')

    # x4, y4 = m4(lon, lat)
    cs_pc = axf[4].pcolormesh(lon, lat, cs_arr, vmin=0, vmax=3, cmap=cs_cmap, norm=cs_norm)
    cb_pc = plt.colorbar(cs_pc, orientation='horizontal', ax=axf[4], fraction=cb_frac, pad=cb_pad)

    axf[4].set_title('Convective/Stratiform')


    # okay now do accumulated rainfall for the last 30 minutes
    tdiffs = np.array([np.abs((_ - file_time).total_seconds()) for _ in accum_times])
    #print 'tdiffs: {}'.format(tdiffs)

    good_files = np.where(tdiffs <= accum_minute_thresh*60.0)
    #print 'good files: {}'.format(good_files)

    accum_rainfall_good_times = np.sum(np.array(accum_rainfall)[good_files], axis=0)

    m5 = Basemap(projection='cyl', resolution=map_res, llcrnrlon=map_limits['lon'][0],  urcrnrlon=map_limits['lon'][1],
                      llcrnrlat=map_limits['lat'][0], urcrnrlat=map_limits['lat'][1], ax=axf[5])

    # m5.drawstates()
    # m5.drawcounties()
    # m5.drawrivers(color='blue')

    # x5, y5 = m5(lon, lat)
    arr_pc = axf[5].pcolormesh(lon, lat, accum_rainfall_good_times, cmap=plt.cm.afmhot_r, vmin=1.0, vmax=25.0)
    arr_pc = plt.colorbar(arr_pc, orientation='horizontal', ax=axf[5], fraction=cb_frac, pad=cb_pad)

    axf[5].set_title('%d min accumulated rainfall (mm)'%(accum_minute_thresh))




    plt.tight_layout()
    plt.savefig('figures/%s_variables.png'%(dt_string))
#    print 'figure saved: {}'.format(time.time() - start)

    plt.close(fig)

    fig, ax = plt.subplots(1,2, figsize=(8,5))

    arr_pc1 = ax[1].pcolormesh(lon, lat, accum_rainfall_good_times, cmap=plt.cm.afmhot_r, vmin=1.0, vmax=25.0)
    arr_cb1 = plt.colorbar(arr_pc1, orientation='horizontal', ax=ax[1], fraction=cb_frac, pad=cb_pad)

    plt.tight_layout()
    plt.savefig('rain_accum_figs/%s_rraccum.png'%(dt_string))

    plt.close(fig)


    # now try the tropical rainfall stuff
    trop_rr, trop_rr_meth = csu_blended_rain_tropical.calc_blended_rain_tropical(dz=dbz, zdr=zdr, kdp=kdp, cs=cs_arr)






    quit()


    del data

    # plt.tight_layout()
    # plt.savefig('figures/%s_cs.png'%(dt_string))














