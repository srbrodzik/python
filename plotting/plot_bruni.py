#!/usr/bin/python3

#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.dates import DayLocator, HourLocator, MinuteLocator, DateFormatter
plt.style.use('seaborn-whitegrid')
import numpy as np
from datetime import datetime
from datetime import timedelta

# From http://lextoumbourou.github.io/blog/posts/get-the-last-day-of-every-month-in-a-range-with-python/
def months_in_range(start_date, end_date):
    """Get the last day of every month in a range between two datetime values.
    Return a generator
    """
    start_month = start_date.month
    end_months = (end_date.year-start_date.year)*12 + end_date.month

    for month in range(start_month, end_months + 1):
        # Get years in the date range, add it to the start date's year
        year = (month-1)//12 + start_date.year
        month = (month-1) % 12 + 1

        yield datetime(year, month, 1)-timedelta(days=1)        

ira = [243049,240911,247988,    # 2017 Q3
       249625,256227,259851,    # 
       266664,249925,245129,    # 2018
       240806,248100,251134,
       260087,265968,265529,
       252463,253483,230052,
       251099,256974,260046,    # 2019
       271546,252419,264204,
       265532,257779,265051,
       268637,281812,284677,
       275386,252098,181038,    # 2020
       205328,208865,213014,
       212747,229529,224171,
       228326,273685,294023,
       306334,339630,356485,    # 2021
       371200,379663,373098,
       370246,377340,371667,
       389521,372150,390330]
roth = [100707, 99817,102743,   # 2017 Q3
        103350,106099,107641,
        110453,103778,101838,   # 2018
         99953,102871,104105,
        107802,110039,109947,
        104730,105051, 95324,
        103652,106000,107219,   # 2019
        112058,104636,109362,
        109855,107072,109885,
        111502,116733,117691,
        114275,104885, 76499,   # 2020
         85928, 87424, 89074,
         88866, 95864, 93707,
         95448,114163,122656,
        127597,141229,147889,   # 2021
        153951,157030,154034,
        152782,155596,152343,
        159375,152563,160146]
        
start_date_text = '20170801'
start_date = datetime.strptime(start_date_text,'%Y%m%d')

end_date_text = '20220101'
end_date = datetime.strptime(end_date_text,'%Y%m%d')

date_list = []
for date in months_in_range(start_date,end_date):
    date_list.append(date)

# Define the date format
#date_form = DateFormatter('%Y/%m')
#ax.xaxis.set_major_formatter(date_form)

#---------
# IRA Plot
#---------
# Create figure and plot space
fig = plt.figure()
ax = plt.axes()

# Plot the data and define labels
# Add this to ax.set if you want to set first and last date plotted
# xlim=['YYYY-MM-DD','YYYY-MM-DD']
plt.plot(date_list,ira)
ax.set(xlabel='Month',
       ylabel='Value ($)',
       title='Bruni IRA - 201707-202112')

# Define the date format
date_form = DateFormatter('%Y/%m')
#date_form = DateFormatter('%Y/%m/%d')
ax.xaxis.set_major_formatter(date_form)

# Ensure a minor tick for every month (this doesn't work)
#ax.xaxis.set_minor_locator(mdates.MonthLocator())
#ax.minorticks_on()
ax.grid(which='major',linestyle=':',color='black')
#ax.grid(which='minor',linestyle=':',color='black')

plt.tight_layout()
#plt.show()
plt.savefig('/home/disk/shear2/brodzik/Bruni_IRA_201707-202112.png')


#--------------
# Roth IRA Plot
#--------------
# Create figure and plot space
fig = plt.figure()
ax = plt.axes()

# Plot the data and define labels
plt.plot(date_list,roth)
ax.set(xlabel='Month',
       ylabel='Value ($)',
       title='Bruni Roth IRA - 201707-202112')

# Define the date format
date_form = DateFormatter('%Y/%m')
#date_form = DateFormatter('%Y/%m/%d')
ax.xaxis.set_major_formatter(date_form)

# Ensure a minor tick for every month (this doesn't work)
#ax.xaxis.set_minor_locator(mdates.MonthLocator())
#ax.minorticks_on()
ax.grid(which='major',linestyle=':',color='black')
#ax.grid(which='minor',linestyle=':',color='black')

plt.tight_layout()
#plt.show()
plt.savefig('/home/disk/shear2/brodzik/Bruni_RIRA_201707-202112.png')

