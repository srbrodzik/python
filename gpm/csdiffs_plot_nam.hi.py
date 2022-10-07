import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator
%matplotlib

months = np.array(['201403','201404','201405','201406','201407','201408','201409','201410','201411','201412',
                   '201501','201502','201503','201504','201505','201506','201507','201508','201509','201510','201511','201512',
                   '201601','201602','201603','201604','201605','201606','201607','201608','201609','201610','201611','201612',
                   '201701','201702','201703','201704','201705','201706','201707','201708'])
num_months = months.size

# Turn months into datetime objects
x_values = []
for i in range(0,num_months):
    year = months[i][0:4]
    month = months[i][4:6]
    curr_month = dt.datetime(int(year),int(month),1)
    x_values.append(curr_month)

#-------------------------------------------------------------------
# CHOOSE ONE OF THESE SETS OF DATA - BS, DWC, DCC, WCC
#BS data
data_nasa = np.array([12,28,35,43,34,45,49,77,48,39,
                      47,44,50,46,31,35,46,27,48,60,58,59,
                      48,44,43,46,41,30,30,32,58,53,51,46,
                      56,34,52,47,51,48,39,40])
#data_uw   = np.array([12,28,35,45,34,46,49,78,51,42,
#                      47,46,51,46,32,36,46,27,51,62,59,60,
#                      49,44,43,46,42,30,33,33,59,56,52,48,
#                      56,34,52,48,51,51,39,41])
data_uw2  = np.array([12,30,38,48,36,48,50,79,53,42,
                      47,49,53,47,35,39,51,30,56,64,61,63,
                      50,46,45,49,44,31,35,36,63,57,56,48,
                      56,36,53,49,56,53,40,43])
#title = 'NAM BS (High Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM BS (High Thresholds)- NASA/UW - v05'

#DWC
data_nasa = np.array([2,5,7,31,19,16,14,3,1,0,
                      1,0,0,11,9,14,31,27,6,2,0,3,
                      2,0,5,3,8,26,30,18,10,4,1,1,
                      0,0,1,14,17,9,32,12])
#data_uw   = np.array([2,5,7,31,20,14,14,3,1,0,
#                      1,0,0,11,9,14,31,27,6,2,0,3,
#                      2,0,5,3,7,26,29,16,10,4,1,1,
#                      0,0,1,14,16,9,31,12])
data_uw2  = np.array([2,5,7,30,20,15,14,3,1,0,
                      1,0,0,11,9,14,28,27,6,2,0,3,
                      2,0,5,3,6,26,30,16,10,4,1,1,
                      0,0,1,14,16,9,32,13])
#title = 'NAM DWC (High Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM DWC (High Thresholds)- NASA/UW - v05'

#DCC
data_nasa = np.array([0,4,14,25,29,47,21,1,1,0,
                      2,0,0,5,22,8,62,39,7,2,4,2,
                      1,0,0,8,14,37,50,29,17,2,2,0,
                      0,1,5,5,13,11,54,33])
#data_uw   = np.array([0,4,13,24,29,48,21,1,1,0,
#                      2,0,0,5,22,8,61,39,7,1,4,2,
#                      1,0,0,8,15,37,50,30,17,2,2,0,
#                      0,1,5,5,14,11,53,33])
data_uw2  = np.array([0,4,14,26,29,48,21,1,1,0,
                      2,0,0,5,22,8,65,39,7,2,4,2,
                      2,0,0,8,15,37,50,31,17,2,2,0,
                      0,1,5,5,14,11,54,33])
#title = 'NAM DCC (High Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM DCC (High Thresholds)- NASA/UW - v05'

#WCC
data_nasa = np.array([16,42,65,82,55,85,59,83,66,38,
                      60,41,23,52,47,62,73,68,95,55,31,60,
                      63,37,44,72,69,84,88,86,75,67,47,24,
                      29,40,54,62,82,85,83,74])
#data_uw   = np.array([13,29,43,67,44,69,52,64,54,15,
#                      47,30,17,39,41,53,63,58,71,46,24,38,
#                      44,25,33,57,49,69,72,77,65,52,31,18,
#                      19,26,42,45,68,64,72,65])
data_uw2  = np.array([9,23,32,53,38,64,45,51,45,9,
                      31,13,10,33,32,43,49,49,59,35,18,27,
                      25,18,17,40,45,57,61,71,55,39,24,13,
                      11,8,24,36,53,53,64,56])
#title = 'NAM WCC (High Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM WCC (High Thresholds)- NASA/UW - v05'
#-------------------------------------------------------------------

# Make plot
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.plot(x_values,data_nasa,'r')
#ax1.plot(x_values,data_uw,'g')
ax1.plot(x_values,data_uw2,'b')

# Set which months to put tickmarks at
ax1.xaxis.set_major_locator( MonthLocator([3,6,9,12,15,18,21,24,27,30,33,36,39,41]))
#set format of months (YYYYMM)
ax1.xaxis.set_major_formatter( DateFormatter('%Y%m') )
#rotate label text to vertical
plt.setp(plt.xticks()[1], rotation=90)

#-- add a legend
#plt.legend(('NASA','UW','UW2'),loc='upper left')
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title(title)
#plt.xlabel('Month')
plt.ylabel('Count')

plt.tight_layout()
plt.show()
