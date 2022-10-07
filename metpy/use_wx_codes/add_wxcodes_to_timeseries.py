#!/usr/bin/python3

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from metpy.plots import StationPlot, current_weather
import numpy as np

# Create some dummy data, including WMO codes for current weather
cur_wx_codes = [65, 85, 95, 99, 54, 57, 59, 67, 88, 82, 79, 75]
dummy_y_vals = np.zeros(len(cur_wx_codes))  # only to place on plot
dates = [datetime.utcnow() + timedelta(hours=i) for i in range(len(cur_wx_codes))]
temps = np.random.randn(len(dates)) * 10 + 70

# Create a plot using matplotlib's StationPlot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sp = StationPlot(ax, dates, dummy_y_vals)
ax.plot(dates, temps)
sp.plot_symbol('C', cur_wx_codes, current_weather, fontsize=16, color='red')
plt.savefig('test.png')
plt.close()

# With my data
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
sp = StationPlot(ax, times, dummy_y_vals)
ax.plot(times, precip_accums)
sp.plot_symbol('C', wxcodes, current_weather, fontsize=16, color='red')
plt.savefig('test.png')
plt.close()

cur_wx_codes = [71,71,71,71,71,71,71,71,71,71,
                71,71,71,71,71,71,71,71,71,71,
                71,71,71,71,71,71,71,71,71,71,
                71,71,71,71,71,71,71,71,71,71,
                71,71,71,71,71,71,71,71,71,71,
                71,71,71,71,71,71,71,71,71,71,
                71,71,71,71,71,71,71,71,71,71,
                71,71]
dummy_y_vals = np.zeros(len(cur_wx_codes))  # only to place on plot
dates = ['2020-02-05 00:00:00', '2020-02-05 01:00:00',
         '2020-02-05 02:00:00', '2020-02-05 03:00:00',
         '2020-02-05 04:00:00', '2020-02-05 05:00:00',
         '2020-02-05 06:00:00', '2020-02-05 07:00:00',
         '2020-02-05 08:00:00', '2020-02-05 09:00:00',
         '2020-02-05 10:00:00', '2020-02-05 11:00:00',
         '2020-02-05 12:00:00', '2020-02-05 13:00:00',
         '2020-02-05 14:00:00', '2020-02-05 15:00:00',
         '2020-02-05 16:00:00', '2020-02-05 17:00:00',
         '2020-02-05 18:00:00', '2020-02-05 19:00:00',
         '2020-02-05 20:00:00', '2020-02-05 21:00:00',
         '2020-02-05 22:00:00', '2020-02-05 23:00:00',
         '2020-02-06 00:00:00', '2020-02-06 01:00:00',
         '2020-02-06 02:00:00', '2020-02-06 03:00:00',
         '2020-02-06 04:00:00', '2020-02-06 05:00:00',
         '2020-02-06 06:00:00', '2020-02-06 07:00:00',
         '2020-02-06 08:00:00', '2020-02-06 09:00:00',
         '2020-02-06 10:00:00', '2020-02-06 11:00:00',
         '2020-02-06 12:00:00', '2020-02-06 13:00:00',
         '2020-02-06 14:00:00', '2020-02-06 15:00:00',
         '2020-02-06 16:00:00', '2020-02-06 17:00:00',
         '2020-02-06 18:00:00', '2020-02-06 19:00:00',
         '2020-02-06 20:00:00', '2020-02-06 21:00:00',
         '2020-02-06 22:00:00', '2020-02-06 23:00:00',
         '2020-02-07 00:00:00', '2020-02-07 01:00:00',
         '2020-02-07 02:00:00', '2020-02-07 03:00:00',
         '2020-02-07 04:00:00', '2020-02-07 05:00:00',
         '2020-02-07 06:00:00', '2020-02-07 07:00:00',
         '2020-02-07 08:00:00', '2020-02-07 09:00:00',
         '2020-02-07 10:00:00', '2020-02-07 11:00:00',
         '2020-02-07 12:00:00', '2020-02-07 13:00:00',
         '2020-02-07 14:00:00', '2020-02-07 15:00:00',
         '2020-02-07 16:00:00', '2020-02-07 17:00:00',
         '2020-02-07 18:00:00', '2020-02-07 19:00:00',
         '2020-02-07 20:00:00', '2020-02-07 21:00:00',
         '2020-02-07 22:00:00', '2020-02-07 23:00:00']

temps = 