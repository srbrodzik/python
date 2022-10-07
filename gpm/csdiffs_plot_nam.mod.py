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
#BS
data_nasa = np.array([20,47,55,71,53,66,65,100,67,53,
                      62,64,66,70,45,53,65,51,70,84,80,75,
                      63,60,58,65,61,54,55,56,89,77,68,56,
                      66,48,67,64,79,71,61,57])
#data_uw   = np.array([20,47,55,74,53,66,66,103,68,56,
#                      63,67,69,70,48,52,69,51,73,84,81,76,
#                      63,60,60,66,63,54,59,57,89,78,68,56,
#                      67,50,67,67,81,72,62,58])
data_uw2  = np.array([21,48,58,81,55,67,70,109,74,56,
                      66,68,70,71,50,57,73,53,75,86,84,78,
                      65,63,64,71,68,56,61,60,92,80,69,56,
                      67,50,68,70,85,78,63,59])
#title = 'NAM BS (Moderate Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM BS (Moderate Thresholds)- NASA/UW - v05'

#DWC
data_nasa = np.array([11,44,126,232,211,299,234,119,78,20,
                      41,17,17,89,140,128,291,286,178,116,55,39,
                      23,12,40,75,108,265,340,272,190,118,58,31,
                      17,17,38,88,133,175,325,263])
#data_uw   = np.array([11,44,120,220,196,283,223,117,76,18,
#                      40,18,18,86,138,126,277,271,171,112,52,39,
#                      25,10,39,71,102,255,321,260,188,113,57,31,
#                      16,15,37,87,132,170,310,250])
data_uw2  = np.array([12,42,115,215,187,253,205,110,70,16,
                      33,18,17,81,124,119,270,257,169,109,48,39,
                      25,10,37,67,98,241,304,241,172,105,53,28,
                      19,14,34,86,128,156,296,235])
#title = 'NAM DWC (Moderate Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM DWC (Moderate Thresholds)- NASA/UW - v05'

#DCC
data_nasa = np.array([26,88,248,701,773,1047,644,234,102,36,
                      50,13,25,184,494,384,1102,1156,516,352,126,41,
                      12,6,64,158,162,815,1261,823,548,248,102,34,
                      10,18,54,86,253,439,1148,1071])
#data_uw   = np.array([25,87,252,711,788,1048,637,235,105,38,
#                      49,13,25,185,494,382,1104,1164,527,354,129,41,
#                      12,7,62,160,167,820,1279,837,550,252,101,34,
#                      11,20,57,89,249,442,1163,1073])
data_uw2  = np.array([26,91,262,731,812,1097,680,244,111,39,
                      57,14,27,191,513,393,1128,1193,528,360,133,41,
                      15,8,67,164,176,837,1300,868,575,259,107,36,
                      12,24,60,93,261,462,1178,1103])
#title = 'NAM DCC (Moderate Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM DCC (Moderate Thresholds)- NASA/UW - v05'

#WCC
data_nasa = np.array([46,104,155,205,135,185,189,193,142,124,
                      149,117,94,132,112,145,172,175,193,190,140,163,
                      169,111,119,137,127,213,190,191,182,241,189,105,
                      91,135,144,112,144,164,168,224])
#data_uw   = np.array([29, 73,117,158,103,139,143,133,110, 73,
#                      114, 93,58, 94, 91,103,127,125,136,141,114,120,
#                      121, 83, 86,105, 87,141,135,145,135,179,126,78,
#                      61,105,114, 75,112,126,120,174])
data_uw2  = np.array([14,33,65,85,71,86,85,80,70,33,
                      59,43,29,52,50,64,72,72,87,89,80,66,
                      66,51,34,62,46,92,83,85,90,110,70,56,
                      28,51,48,41,77,70,73,116])
#title = 'NAM WCC (Moderate Thresholds)- NASA/UW/UW2 - v05'
title = 'NAM WCC (Moderate Thresholds)- NASA/UW - v05'
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
