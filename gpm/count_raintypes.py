import os
import netCDF4 as nc4
import numpy as np
import datetime as dt
from matplotlib.dates import DateFormatter, MonthLocator

# For plotting
%matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt

STRA = 1
CONV = 2
OTHER = 3

NHW = 0
NHC = 1
NHE = 2

#User inputs
indir_v05 = '/home/disk/bob/gpm/nht_ku/classify/ex_data_v05'
indir_ite132 = '/home/disk/bob/gpm/nht_ku/classify/ex_data_ite132'
indir_ite601 = '/home/disk/bob/gpm/nht_ku/classify/ex_data_ite601'

months = np.array(['201504','201505','201506','201507','201508','201509',
                   '201510','201511','201512','201601','201602','201603'])

#Declare output arrays
stra_counts_v05 = np.zeros((3,12),dtype=int)
stra_counts_v05_uw = np.zeros((3,12),dtype=int)
stra_counts_ite132 = np.zeros((3,12),dtype=int)
stra_counts_ite601 = np.zeros((3,12),dtype=int)

conv_counts_v05 = np.zeros((3,12),dtype=int)
conv_counts_v05_uw = np.zeros((3,12),dtype=int)
conv_counts_ite132 = np.zeros((3,12),dtype=int)
conv_counts_ite601 = np.zeros((3,12),dtype=int)

other_counts_v05 = np.zeros((3,12),dtype=int)
other_counts_v05_uw = np.zeros((3,12),dtype=int)
other_counts_ite132 = np.zeros((3,12),dtype=int)
other_counts_ite601 = np.zeros((3,12),dtype=int)

# array for datetime objects used for x-axis labels
x_values = []
#len(x_values)

num_months = months.size
for i in range(0,num_months):
    year = months[i][0:4]
    month = months[i][4:6]
    print 'year=',year,' and month=',month
 
    # turn month into datetime object
    curr_month = dt.datetime(int(year),int(month),1)
    x_values.append(curr_month)

    print 'Going through v05 data . . .'
    for fname in os.listdir(indir_v05+'/'+year+'/'+month):
        #print fname
        #open fname and read rain_type and rain_type_uw
        ncid = nc4.Dataset(indir_v05+'/'+year+'/'+month+'/'+fname,'r')
        rt = np.array(ncid.variables['rain_type'])
        rt_uw = np.array(ncid.variables['rain_type_uw'])
        ncid.close()
        if fname.endswith('NHW.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_v05[NHW,i] += np.count_nonzero(rt == STRA)
            conv_counts_v05[NHW,i] += np.count_nonzero(rt == CONV)
            other_counts_v05[NHW,i] += np.count_nonzero(rt == OTHER)
            stra_counts_v05_uw[NHW,i] += np.count_nonzero(rt_uw == STRA)
            conv_counts_v05_uw[NHW,i] += np.count_nonzero(rt_uw == CONV)
            other_counts_v05_uw[NHW,i] += np.count_nonzero(rt_uw == OTHER)
        elif fname.endswith('NHC.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_v05[NHC,i] += np.count_nonzero(rt == STRA)
            conv_counts_v05[NHC,i] += np.count_nonzero(rt == CONV)
            other_counts_v05[NHC,i] += np.count_nonzero(rt == OTHER)
            stra_counts_v05_uw[NHC,i] += np.count_nonzero(rt_uw == STRA)
            conv_counts_v05_uw[NHC,i] += np.count_nonzero(rt_uw == CONV)
            other_counts_v05_uw[NHC,i] += np.count_nonzero(rt_uw == OTHER)
        elif fname.endswith('NHE.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_v05[NHE,i] += np.count_nonzero(rt == STRA)
            conv_counts_v05[NHE,i] += np.count_nonzero(rt == CONV)
            other_counts_v05[NHE,i] += np.count_nonzero(rt == OTHER)
            stra_counts_v05_uw[NHE,i] += np.count_nonzero(rt_uw == STRA)
            conv_counts_v05_uw[NHE,i] += np.count_nonzero(rt_uw == CONV)
            other_counts_v05_uw[NHE,i] += np.count_nonzero(rt_uw == OTHER)
    print 'Going through ite132 data . . .'
    for fname in os.listdir(indir_ite132+'/'+year+'/'+month):
        #print fname
        #open fname and read rain_type
        ncid = nc4.Dataset(indir_ite132+'/'+year+'/'+month+'/'+fname,'r')
        rt = np.array(ncid.variables['rain_type'])
        ncid.close()
        if fname.endswith('NHW.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_ite132[NHW,i] += np.count_nonzero(rt == STRA)
            conv_counts_ite132[NHW,i] += np.count_nonzero(rt == CONV)
            other_counts_ite132[NHW,i] += np.count_nonzero(rt == OTHER)
        elif fname.endswith('NHC.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_ite132[NHC,i] += np.count_nonzero(rt == STRA)
            conv_counts_ite132[NHC,i] += np.count_nonzero(rt == CONV)
            other_counts_ite132[NHC,i] += np.count_nonzero(rt == OTHER)
        elif fname.endswith('NHE.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_ite132[NHE,i] += np.count_nonzero(rt == STRA)
            conv_counts_ite132[NHE,i] += np.count_nonzero(rt == CONV)
            other_counts_ite132[NHE,i] += np.count_nonzero(rt == OTHER)
    print 'Going through ite601 data . . .'
    for fname in os.listdir(indir_ite601+'/'+year+'/'+month):
        #print fname
        #open fname and read rain_type
        ncid = nc4.Dataset(indir_ite601+'/'+year+'/'+month+'/'+fname,'r')
        rt = np.array(ncid.variables['rain_type'])
        ncid.close()
        if fname.endswith('NHW.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_ite601[NHW,i] += np.count_nonzero(rt == STRA)
            conv_counts_ite601[NHW,i] += np.count_nonzero(rt == CONV)
            other_counts_ite601[NHW,i] += np.count_nonzero(rt == OTHER)
        elif fname.endswith('NHC.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_ite601[NHC,i] += np.count_nonzero(rt == STRA)
            conv_counts_ite601[NHC,i] += np.count_nonzero(rt == CONV)
            other_counts_ite601[NHC,i] += np.count_nonzero(rt == OTHER)
        elif fname.endswith('NHE.nc'):
            #unique, counts = np.unique(rt, return_counts=True)
            stra_counts_ite601[NHE,i] += np.count_nonzero(rt == STRA)
            conv_counts_ite601[NHE,i] += np.count_nonzero(rt == CONV)
            other_counts_ite601[NHE,i] += np.count_nonzero(rt == OTHER)
                
# plot convective output
fig = plt.figure()
ax = fig.gca()
#plt.xlim(1, 12)
#plt.ylim();
plt.plot(x_values,conv_counts_v05[0,:], '-', color='black', label='NHW_v05')
plt.plot()
plt.plot(x_values,conv_counts_v05_uw[0,:], '-', color='red', label='NHW_v05_uw')
plt.plot()
plt.plot(x_values,conv_counts_ite132[0,:], '-', color='blue', label='NHW_ite132')
plt.plot()
plt.plot(x_values,conv_counts_ite601[0,:], '-', color='green', label='NHW_ite601')
plt.plot()
plt.plot(x_values,conv_counts_v05[1,:], ':', color='black', label='NHC_v05')
plt.plot()
plt.plot(x_values,conv_counts_v05_uw[1,:], ':', color='red', label='NHC_v05_uw')
plt.plot()
plt.plot(x_values,conv_counts_ite132[1,:], ':', color='blue', label='NHC_ite132')
plt.plot()
plt.plot(x_values,conv_counts_ite601[1,:], ':', color='green', label='NHC_ite601')
plt.plot()
plt.plot(x_values,conv_counts_v05[2,:], '--', color='black', label='NHE_v05')
plt.plot()
plt.plot(x_values,conv_counts_v05_uw[2,:], '--', color='red', label='NHE_v05_uw')
plt.plot()
plt.plot(x_values,conv_counts_ite132[2,:], '--', color='blue', label='NHE_ite132')
plt.plot()
plt.plot(x_values,conv_counts_ite601[2,:], '--', color='green', label='NHE_ite601')
plt.plot()
plt.grid(linestyle=":")
plt.legend()
# Set which months to put tickmarks at
ax.set_xlim([x_values[0], x_values[-1]])
#ax.xaxis.set_major_locator( MonthLocator(0,1,2,3,4,5,6,7,8,9,10,11))
# Set format of months (YYYYMM)
#ax.xaxis.set_major_formatter( DateFormatter('%Y%m') )
# Rotate label text to vertical
plt.setp(plt.xticks()[1], rotation=45)
plt.xlabel('months')
plt.ylabel('pixel counts')
plt.title('Convective Counts - Near Equatorial, Northern Hemisphere')
plt.tight_layout()

# plot stratiform output
fig = plt.figure()
ax = fig.gca()
#plt.xlim(1, 12)
#plt.ylim();
plt.plot(x_values,stra_counts_v05[0,:], '-', color='black', label='NHW_v05')
plt.plot()
plt.plot(x_values,stra_counts_v05_uw[0,:], '-', color='red', label='NHW_v05_uw')
plt.plot()
plt.plot(x_values,stra_counts_ite132[0,:], '-', color='blue', label='NHW_ite132')
plt.plot()
plt.plot(x_values,stra_counts_ite601[0,:], '-', color='green', label='NHW_ite601')
plt.plot()
plt.plot(x_values,stra_counts_v05[1,:], ':', color='black', label='NHC_v05')
plt.plot()
plt.plot(x_values,stra_counts_v05_uw[1,:], ':', color='red', label='NHC_v05_uw')
plt.plot()
plt.plot(x_values,stra_counts_ite132[1,:], ':', color='blue', label='NHC_ite132')
plt.plot()
plt.plot(x_values,stra_counts_ite601[1,:], ':', color='green', label='NHC_ite601')
plt.plot()
plt.plot(x_values,stra_counts_v05[2,:], '--', color='black', label='NHE_v05')
plt.plot()
plt.plot(x_values,stra_counts_v05_uw[2,:], '--', color='red', label='NHE_v05_uw')
plt.plot()
plt.plot(x_values,stra_counts_ite132[2,:], '--', color='blue', label='NHE_ite132')
plt.plot()
plt.plot(x_values,stra_counts_ite601[2,:], '--', color='green', label='NHE_ite601')
plt.plot()
plt.grid(linestyle=":")
plt.legend()
# Set which months to put tickmarks at
ax.set_xlim([x_values[0], x_values[-1]])
#ax.xaxis.set_major_locator( MonthLocator(0,1,2,3,4,5,6,7,8,9,10,11))
# Set format of months (YYYYMM)
#ax.xaxis.set_major_formatter( DateFormatter('%Y%m') )
# Rotate label text to vertical
plt.setp(plt.xticks()[1], rotation=45)
plt.xlabel('months')
plt.ylabel('pixel counts')
plt.title('Stratiform Counts - Near Equatorial, Northern Hemisphere')
plt.tight_layout()

# plot other output (sanity check -- should all be the same)
fig = plt.figure()
ax = fig.gca()
#plt.xlim(1, 12)
#plt.ylim();
plt.plot(x_values,other_counts_v05[0,:], '-', color='black', label='NHW_v05')
plt.plot()
plt.plot(x_values,other_counts_v05_uw[0,:], '-', color='red', label='NHW_v05_uw')
plt.plot()
plt.plot(x_values,other_counts_ite132[0,:], '-', color='blue', label='NHW_ite132')
plt.plot()
plt.plot(x_values,other_counts_v05[1,:], ':', color='black', label='NHC_v05')
plt.plot()
plt.plot(x_values,other_counts_v05_uw[1,:], ':', color='red', label='NHC_v05_uw')
plt.plot()
plt.plot(x_values,other_counts_ite132[1,:], ':', color='blue', label='NHC_ite132')
plt.plot()
plt.plot(x_values,other_counts_v05[2,:], '--', color='black', label='NHE_v05')
plt.plot()
plt.plot(x_values,other_counts_v05_uw[2,:], '--', color='red', label='NHE_v05_uw')
plt.plot()
plt.plot(x_values,other_counts_ite132[2,:], '--', color='blue', label='NHE_ite132')
plt.plot()
plt.grid(linestyle=":")
plt.legend()
# Set which months to put tickmarks at
#ax.xaxis.set_major_locator( MonthLocator(0,1,2,3,4,5,6,7,8,9,10,11))
# Set format of months (YYYYMM)
#ax.xaxis.set_major_formatter( DateFormatter('%Y%m') )
# Rotate label text to vertical
#plt.setp(plt.xticks()[1], rotation=45)
plt.xlabel('months')
plt.ylabel('pixel counts')
plt.title('Other Counts - Near Equatorial, Northern Hemisphere')
plt.tight_layout()

