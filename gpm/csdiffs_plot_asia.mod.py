import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, date, time, timedelta
import pandas as pd

# use magic command to allow multiple plots
%matplotlib

months = ('201403','201404','201405','201406','201407','201408','201409','201410','201411','201412',
          '201501','201502','201503','201504','201505','201506','201507','201508','201509','201510','201511','201512',
          '201601','201602','201603','201604')
#          '201601','201602','201603','201604','201605','201606','201607','201608','201609','201610','201611','201612',
#          '201701','201702','201703','201704','201705','201706','201707','201708')
t = np.arange(np.size(months))
#t = np.array( np.size(months),np.dtype(datetime) )

#-----------------------------------------------------------------------------
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
#-----------------------------------------------------------------------------

#BS
data_nasa = np.array([31,71,121,127,158,122,114,112,73,61,
                      57,51,61,86,118,135,168,156,111,112,75,67,
                      53,50,70,92])
#                      53,50,70,92,130,101,173,168,107,81,85,77,
#                      55,68,72,92,134,146,154,140])
data_uw   = np.array([32,72,122,128,160,122,115,113,74,61,
                      58,51,62,87,119,140,172,159,113,112,75,67,
                      53,51,70,93])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA BS (Moderate Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

#DWC
data_nasa = np.array([22,79,184,413,376,419,239,83,51,5,
                      1,7,30,106,212,305,501,368,192,116,46,8,
                      6,11,52,100])
#                      6,11,52,100,170,450,492,477,274,141,44,17,
#                      8,10,54,98,202,340,573,407])
data_uw   = np.array([22,75,173,396,355,404,228,81,51,5,
                      2,7,30,104,207,291,485,346,185,107,46,6,
                      5,11,48,96])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA DWC (Moderate Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

#DCC
data_nasa = np.array([81,270,536,1737,1468,1626,802,316,129,8,
                      3,32,89,283,716,764,1955,1410,478,342,96,9,
                      6,44,149,278])
#                      6,44,149,278,407,1803,2149,1384,873,321,78,20,
#                      5,12,102,243,631,989,2204,2009])
data_uw   = np.array([81,275,542,1743,1464,1608,793,313,128,8,
                      3,31,88,285,710,757,1940,1406,472,350,95,11,
                      7,40,153,281])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA DCC (Moderate Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

#WCC
data_nasa = np.array([42,102,174,328,346,368,213,204,169,142,
                      132,120,95,110,189,284,493,388,211,191,157,142,
                      97,78,95,129])
#                      97,78,95,129,174,359,377,414,257,258,172,143,
#                      137,94,146,87,185,286,373,391])
data_uw   = np.array([32,76,133,253,270,292,168,153,122,93,
                      97,81,70,79,144,232,404,291,155,151,116,104,
                      62,45,65,85])

plt.plot(t,data_nasa,'r')
plt.plot(t,data_uw,'b')

#-- add a legend
plt.legend(('NASA','UW'),loc='upper left')

#-- add a annotations
plt.title('ASIA WCC (Moderate Thresholds)- NASA vs UW - v05')

plt.xlabel('YYYYMM')
plt.ylabel('Count')

plt.show()

