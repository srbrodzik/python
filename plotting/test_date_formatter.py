'''
test_date_formatter.py
Sample plot with YYYYMM plotting
'''

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MonthLocator
%matplotlib

#first make some fake data
x_values = []
monthnow = dt.datetime(2017,10,1)

#make list of incremental months
for i in range(0,40):
    x_values.append(monthnow)
    monthnow+=dt.timedelta(days=35)
    monthnow = monthnow.replace(day=1)

#random y-values to go with our months
y_values = np.random.randn(40)

#make plot
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.plot(x_values,y_values)

#set which months to put tick marks at
ax1.xaxis.set_major_locator( MonthLocator([1,3,5,7,9,11]))

#set format of months (YYYYMM)
ax1.xaxis.set_major_formatter( DateFormatter('%Y-%m') )

#rotate label text to vertical
plt.setp(plt.xticks()[1], rotation=90)

plt.tight_layout()
plt.show()

