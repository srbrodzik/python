import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, date, time, timedelta
import pandas as pd

# use magic command to allow multiple plots
%matplotlib

months = ('201403','201404','201405','201406','201407','201408')
#months = ('201403','201404','201405','201406','201407','201408','201409','201410','201411','201412',
#          '201501','201502','201503','201504','201505','201506','201507','201508','201509','201510','201511','201512',
#          '201601','201602','201603','201604','201605','201606','201607','201608','201609','201610','201611','201612',
#          '201701','201702','201703','201704','201705','201706','201707','201708')
t = np.arange(np.size(months))
#t = np.array( np.size(months),np.dtype(datetime) )

#-------------------------------------------------------------------------
t1 = datetime.strptime(t)
dates = matplotlib.dates.date2num(datetimes)
plt.plot_date(dates,data1)

#t = [datetime.datetime(2017, 8, 17, 16, 56, 25, 817369),
# datetime.datetime(2017, 8, 17, 17, 56, 25, 817387),
# ...]
i = 0
for month in months:
    print month
    t[i] = datetime.strptime(month,'%Y%m')
    print t
    i=i+1

base = datetime(2014, 3, 1)
arr = np.array([base + timedelta(months=i) for i in xrange(10)])
t = np.arange(datetime(1985,7,1), datetime(2015,7,1), timedelta(days=1)).astype(datetime)
t = pd.to_datetime(['201403','201404','201405','201406','201407','201408','201409','201410','201411','201412'])

t = datetimes
plot(t,data1)
plt.gcf().autofmt_xdate()
#-------------------------------------------------------------------------

#BS
data_nasa = np.array([25,51,89,80,95,86])
#data_nasa = np.array([25,51,89,80,95,86,71,72,59,48,
#                      44,34,45,69,80,89,127,98,70,79,52,55,
#                      45,39,52,68,94,63,119,116,75,57,65,63,
#                      44,51,58,70,84,100,108,89])
data_uw   = np.array([25,52,90,82,98,87])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA BS (High Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

#DWC
data_nasa = np.array([1,8,19,12,13,8])
#data_nasa = np.array([1,8,19,12,13,8,4,0,0,0,
#                      0,1,2,12,17,13,21,17,2,1,1,0,
#                      1,0,6,7,16,24,17,9,4,6,0,0,
#                      0,0,7,7,17,15,32,14])
data_uw   = np.array([1,7,19,12,13,7])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA DWC (High Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

#DCC
data_nasa = np.array([7,16,34,40,28,26])
#data_nasa = np.array([7,16,34,40,28,26,16,1,1,0,
#                      0,0,4,17,28,21,52,36,14,4,1,0,
#                      1,3,11,9,23,42,46,34,6,7,0,0,
#                      0,1,8,22,35,34,72,39])
data_uw   = np.array([6,17,33,39,26,26])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA DCC (High Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

#WCC
data_nasa = np.array([15,34,71,141,109,150])
#data_nasa = np.array([15,34,71,141,109,150,91,48,57,39,
#                      30,15,23,40,93,144,190,106,72,53,44,36,
#                      26,13,24,72,88,134,131,168,94,75,52,48,
#                      37,32,48,42,60,133,151,120])
data_uw   = np.array([12,23,54,113,97,128])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA WCC (High Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

