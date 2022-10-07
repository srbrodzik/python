# Imports and functions
import os
import netCDF4 as nc4
import numpy as np
#import numpy.ma as ma
#import sys

# For plotting
%matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

# User inputs
indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05_testing'
#indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v05'
years = ['2016']
months = ['12']
VERT_RES = 0.125  #km
STRA = 1
CONV = 2

# ORIG
#THRESHOLD = 38.0  #dBZ
#DELTA_HT_FOR_SLOPE = 1.5  #km
#SLOPE = 8.0
#MAX_BB_HT = 5.0   #km
#REFL_MAX_OFFSET = 2.0 #dBZ

# TRY1
#THRESHOLD = 0.0  #dBZ
#DELTA_HT_FOR_SLOPE = 1.0  #km
#SLOPE = 5.0
#MAX_BB_HT = 5.0   #km
#REFL_MAX_OFFSET = 2.0 #dBZ

# TRY2
THRESHOLD = 0.0  #dBZ
DELTA_HT_FOR_SLOPE = 1.5  #km
SLOPE = 7.5
MAX_BB_HT = 5.0   #km
REFL_MAX_OFFSET = 2.0 #dBZ

for iyear in years:
    for imonth in months:
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
            if fname.endswith('nc'):

                iyear = years[0]
                imonth = months[0]

                fname1 = 'GPM2Ku5_uw3_20140630.214655_to_20140630.215308_001917_NAM.nc'
                fname2 = 'GPM2Ku5_uw3_20161212.214204_015857_NAM.nc'
                fname3 = 'GPM2Ku5_uw3_20161213.172250_015870_NAM.nc'
                fname4 = 'GPM2Ku5_uw3_20161215.021118_015891_NAM.nc'
                fname5 = 'GPM2Ku5_uw3_20161225.195121_016058_NAM.nc'
                fname6 = 'GPM2Ku5_uw3_20170801.201702_to_20170801.201950_019465_NAM.nc'
                
                fname = fname6
                print fname

                #ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'a')
                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'r')

                refl = np.array(ncid.variables['refl'])
                #refl = refl[:,1:,:,:]  # get rid of lowest level due to suspicious data
                (ntime,nalt,nlat,nlon) = refl.shape
                #refl_missing = ncid.variables['refl'].missing_value
                refl_missing = ncid.variables['refl']._FillValue
                refl[:,0,:,:] = refl_missing  # get rid of lowest level due to suspicious data
                bin_missing = int(refl_missing)

                rt = np.array(ncid.variables['rain_type'])
                #rt_missing = ncid.variables['rain_type'].missing_value
                rt_missing = ncid.variables['rain_type']._FillValue
                #rt[(rt==rt_missing)] = np.nan

                rt_uw = np.array(rt, copy=True)

                srt = np.array(ncid.variables['shallow_rain_type'])
                #srt_missing = ncid.variables['shallow_rain_type'].missing_value
                srt_missing = ncid.variables['shallow_rain_type']._FillValue
                #srt[(srt==srt_missing)] = np.nan

                # find max_refl and max_refl_ht at each lat/lon
                max_refl = np.zeros(rt.shape,dtype=float)+refl_missing
                max_refl = np.amax(refl,axis=1)
                max_refl_bin = np.zeros(rt.shape,dtype=int)+bin_missing
                max_refl_bin = np.argmax(refl,axis=1)
                max_refl_ht = max_refl_bin * VERT_RES

                max_refl_bin[(max_refl==refl_missing)] = bin_missing
                max_refl_ht[(max_refl==refl_missing)] = refl_missing

                # initialize reference refl/bin/ht
                # these reflect highest refl,bin,ht above max_refl,bin.ht and within +/- REFL_MAX_OFFSET
                reference_refl = np.array(max_refl, copy=True)
                reference_refl_bin = np.array(max_refl_bin, copy=True)
                reference_refl_ht = np.array(max_refl_ht, copy=True)
                
                # find bin & height of highest good refl at each lat/lon
                max_ht_bin = np.zeros(rt.shape,dtype=int) + bin_missing
                for lev in range(nalt-1,-1,-1):
                    max_ht_bin = np.where(np.logical_and((refl[:,lev,:,:] != refl_missing),
                                                         (max_ht_bin == bin_missing) ), 
                                          lev, max_ht_bin )
                max_ht = max_ht_bin * (max_ht_bin != bin_missing) * VERT_RES
                max_ht[(max_ht_bin==bin_missing)] = refl_missing

                # For rt_uw, change rain_type to stratiform if:
                # 1. rain_type = convective and
                # 2. shallow_rain_type = 0 and
                # 3. max_refl in column >= THRESHOLD and   ## DON'T USE THIS CRITERIA ANYMORE
                # 4. [(dbz @ reference point) - (dbz @ point 1.5km higher)]/
                #    [(ht  @ reference point) - (ht  @ point 1.5km higher)] >= abs(SLOPE)
                #  as long as ht @ reference point <= MAX_BB_HT
                
                # Determine indices of interest 
                #change = np.logical_and(np.logical_and(rt==2,srt==0),max_refl>=THRESHOLD)
                change = np.logical_and(rt==2,srt==0)
                ind_time,ind_lat,ind_lon = np.where(change==True)
                npixels = ind_time.size

                # for each index
                #    start at max_refl_ht & go up to find reference refl within +/- REFL_MAX_OFFSET
                #    if height of new reference refl is less than MAX_BB_HT
                #        determine slope between that point and point DELTA_HT_FOR_SLOPE higher or highest good point in column
                #        if slope exceeds SLOPE, change raintype from convective to stratiform                    
                for ind in range(0,npixels):
                    #itime = ind_time[ind]
                    #ilat = ind_lat[ind]
                    #ilon = ind_lon[ind]
                    # CONV - fname1
                    #itime = 0
                    #ilat = 36
                    #ilon = 308
                    # STRA - fname6
                    #itime = 0
                    #ilat = 27
                    #ilon = 81
                    # TWEENER - fname2
                    #itime = 0
                    #ilat = 190
                    #ilon = 128
                    # TWEENER - fname3
                    #itime = 0
                    #ilat = 197
                    #ilon = 138
                    # TWEENER - fname4
                    #itime = 0
                    #ilat = 489
                    #ilon = 89
                    # TWEENER - fname5
                    #itime = 0
                    #ilat = 137
                    #ilon = 81
                    # TWEENER1 - fname6
                    #itime = 0
                    #ilat = 38
                    #ilon = 141
                    # TWEENER2 - fname6
                    itime = 0
                    ilat = 39
                    ilon = 166
                    
                    if max_refl[itime,ilat,ilon] != refl_missing:
                        #print('Max refl point: ',max_refl[itime,ilat,ilon],max_refl_bin[itime,ilat,ilon],max_refl_ht[itime,ilat,ilon])
                        
                        next_bin = reference_refl_bin[itime,ilat,ilon]+1
                        next_ht = reference_refl_ht[itime,ilat,ilon]+VERT_RES
                        next_refl = refl[itime,next_bin,ilat,ilon]
                        #print('Reference point:',reference_refl[itime,ilat,ilon],reference_refl_bin[itime,ilat,ilon],reference_refl_ht[itime,ilat,ilon])
                        #print(next_refl,next_bin,next_ht)
                        while abs(max_refl[itime,ilat,ilon]-next_refl) <= REFL_MAX_OFFSET:
                            reference_refl[itime,ilat,ilon]=next_refl
                            reference_refl_bin[itime,ilat,ilon] = next_bin
                            reference_refl_ht[itime,ilat,ilon]=next_ht  
                            next_bin = reference_refl_bin[itime,ilat,ilon]+1
                            next_ht = reference_refl_ht[itime,ilat,ilon]+VERT_RES
                            next_refl = refl[itime,next_bin,ilat,ilon]
                            #print('Reference point:',reference_refl[itime,ilat,ilon],reference_refl_bin[itime,ilat,ilon],reference_refl_ht[itime,ilat,ilon])
                            #print(next_refl,next_bin,next_ht)

                        #print('Reference point: ',reference_refl[itime,ilat,ilon],reference_refl_bin[itime,ilat,ilon],reference_refl_ht[itime,ilat,ilon])
                        
                        if reference_refl_ht[itime,ilat,ilon] < MAX_BB_HT:
                            if (reference_refl_ht[itime,ilat,ilon] + DELTA_HT_FOR_SLOPE) <= max_ht[itime,ilat,ilon]:
                                compare_bin = reference_refl_bin[itime,ilat,ilon] + int(DELTA_HT_FOR_SLOPE/VERT_RES)
                                compare_ht = reference_refl_ht[itime,ilat,ilon] + DELTA_HT_FOR_SLOPE
                                compare_refl = refl[itime,compare_bin,ilat,ilon]
                            else:
                                compare_bin = max_ht_bin[itime,ilat,ilon]
                                compare_ht = max_ht[itime,ilat,ilon]
                                compare_refl = refl[itime,compare_bin,ilat,ilon]
                            #print('Compare point: ',compare_refl,compare_bin,compare_ht)
                            if (abs(reference_refl_ht[itime,ilat,ilon] - compare_ht)) != 0.0:
                                if ( abs(reference_refl[itime,ilat,ilon] - compare_refl)/
                                     abs(reference_refl_ht[itime,ilat,ilon] - compare_ht) ) >= SLOPE:
                                    rt_uw[itime,ilat,ilon] = STRA

                # PLOT COLUMN

                #some computations
                delta_refl = reference_refl[itime,ilat,ilon] - compare_refl
                delta_ht = compare_ht - reference_refl_ht[itime,ilat,ilon]
                slope = delta_refl/delta_ht
                
                #y = range(1, 30)
                x = refl[0,0:54,ilat,ilon]
                y = np.linspace(0.0,6.625,num=54)
                fig = plt.figure()
                plt.plot(x, y, '-',color='black')
                #plt.plot(x, y, 'o')
                #plt.xlim(15, 50);
                plt.xlim(15, 60);
                plt.ylim(0.0, 7.0);
                plt.xlabel('refl (dBZ)')
                plt.ylabel('ht (km)')

                x1 = range(15,61)
                y1 = np.array([MAX_BB_HT for i in xrange(len(x1))])
                #plt.plot(x1, y1, 'r-',label='MAX_BB_HT') 
                #plt.plot()

                x2 = range(15,61)
                y2 = np.array([max_refl_ht[itime,ilat,ilon] for i in xrange(len(x2))])
                plt.plot(x2, y2, 'm-', label='max refl pt') 
                plt.plot()

                x3 = range(15,61)
                y3 = np.array([reference_refl_ht[itime,ilat,ilon] for i in xrange(len(x3))])
                plt.plot(x3, y3, 'b--', label='ref refl pt') 
                plt.plot()

                x4 = range(15,61)
                y4 = np.array([compare_ht for i in xrange(len(x4))])
                plt.plot(x4, y4, 'g--', label='compare refl pt') 
                plt.plot()

                y5 = range(0,8)
                x5 = np.array([reference_refl[itime,ilat,ilon] for i in xrange(len(y5))])
                plt.plot(x5, y5, 'b--') 
                plt.plot()

                y6 = range(0,8)
                x6 = np.array([compare_refl for i in xrange(len(y6))])
                plt.plot(x6, y6, 'g--') 
                plt.plot()

                plt.legend()

                plt.fill_between(x1,y1,color='gray',alpha=0.3)

                # CONV-file1
                plt.title('Convective Reflectivity Profile - 20140630/2147')
                # STRA-file1
                plt.title('Stratiform Reflectivity Profile - 20170801/2018')
                # TWEENER-file2/3/4/5/6
                plt.title('Reclassified Stratiform Reflectivity Profile - 20140630/2147')
                plt.title('Reclassified Stratiform Reflectivity Profile - 20161212/2143')
                plt.title('Reclassified Stratiform Reflectivity Profile - 20161213/1723')
                plt.title('Reclassified Stratiform Reflectivity Profile - 20161215/0212')
                plt.title('Reclassified Stratiform Reflectivity Profile - 20161225/1952')
                plt.title('Reclassified Stratiform Reflectivity Profile - 20170801/2018')
                plt.text(18,reference_refl_ht[itime,ilat,ilon]+delta_ht/2,'Delta_ht='+str(round(delta_ht,2)))
                plt.text(compare_refl+delta_refl/2,5,'Delta_refl='+str(round(delta_refl,2)),rotation=90)
                plt.text(18,0.5,'SLOPE='+str(round(delta_refl/delta_ht,2)))
                
                #TESTING - reset some conv pixels to stra
                rt_uw[itime,127,93] = CONV
                #rt_uw[itime,126,93] = CONV
                #rt_uw[itime,126,91] = CONV

                # Replace old UW raintype field with new one in input netcdf file
                ncid.variables['rain_type_uw'][:] = rt_uw
    
                ncid.close()

# FOR TESTING DATA
#TESTING
refl[refl < 0.0] = np.nan
fig = plt.figure
plt.matshow(refl[0,18,:,:])
plt.colorbar(orientation='vertical')
