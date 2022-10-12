# Copyright (c) 2016 MetPy Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause

# Modified by S Brodzik, Univ of Washington, 11 Oct 2022

"""
==========
plot_skewt
==========

Layout a Skew-T plot with an optional hodograph inset into the plot.
"""

import sys
sys.path.append('/home/disk/meso-home/brodzik/python/metpy/skewt')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, Hodograph, SkewT
from metpy.units import units
from sndg_io import read_infile
from calculations import calcs
import matplotlib.patheffects as mpatheffects
import numpy as np

#degC = '$^{\circ}$C'
degC = '\N{DEGREE CELSIUS}'
hodo = False

if len(sys.argv) != 4:
    print('Usage: {} [input path] [input file] [format (CSU|cdf|MUnc|MUtxt_ws|raw|SBUnc|SBUnc_static|SCOUT|UIUC|UIUCnc|UWYO]'.format(sys.argv[0]))
    sys.exit()
else:
    inpath = sys.argv[1]
    infile = sys.argv[2]
    format = sys.argv[3]
    print('indir =',indir,'infile =',infile,'and format =',format)

# Test inputs
inpath = '/home/disk/bob/impacts/upperair/uill/20220117'
infile = 'upperair.UILL_sonde.202201170600.nc'
format = 'UIUCnc'

###########################################
# Read sounding data into dataframe with column names as indicated
# col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']

(df,out_fname,figtitle) = read_infile(inpath,infile,format)

###########################################
# We will pull the data out of the example dataset into individual variables and
# assign units.

hght = df['height'].values * units.hPa
p = df['pressure'].values * units.hPa
T = df['temperature'].values * units.degC
Td = df['dewpoint'].values * units.degC
wind_speed = df['speed'].values * units.knots
wind_dir = df['direction'].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir)

###########################################
# Calculate integer heights at specified pressure levels
# Input pressure vals must be increasing so I flip them

p_vals = np.flip(df['pressure'].values)
h_vals = np.flip(df['height'].values)
p_levels = [100,200,300,400,500,600,700,800,900,1000]
h_levels = np.interp(p_levels,p_vals,h_vals)

###########################################
# Create a new figure. The dimensions here give a good aspect ratio
fig = plt.figure(figsize=(10, 9))

# Write data values in right column of plot
fig.text(0.86,0.85,'Sounding Params', fontsize=10)
(h0,h10,h20,h30,h40,wcd,pw,lcl_height,shear6km,shear850200mb,shear700mb) = calcs(df)
fig.text(0.86,0.83,  '0%s: %.0f m'%(degC, h0))
fig.text(0.86,0.81,'-10%s: %.0f m'%(degC, h10))
fig.text(0.86,0.79,'-20%s: %.0f m'%(degC, h20))
fig.text(0.86,0.77,'-30%s: %.0f m'%(degC, h30))
fig.text(0.86,0.75,'-40%s: %.0f m'%(degC, h40))
fig.text(0.86,0.73,'sfc LCL: %.0f m'%lcl_height)
fig.text(0.86,0.71,'WCD: %d m'%wcd)
fig.text(0.86,0.69,'PW: %.1f mm'%pw)
fig.text(0.86,0.67,'6 km shear: %d kts'%shear6km)
fig.text(0.86,0.63,'850-200mb shear: \n       %d kts'%shear850200mb)
fig.text(0.86,0.59,'SFC-700mb shear: \n       %d kts'%shear700mb)
fig.text(0.88,0.20,'Surface')
fig.text(0.88,0.18,'P: %.1f hPa'%df['pressure'][0])
fig.text(0.88,0.16,'H: %.0f m'%(df['height'][0]))
fig.text(0.88,0.14,'T: %.1f %s'%(df['temperature'][0], degC))
fig.text(0.88,0.12,'T$_{D}$: %.1f %s'%(df['dewpoint'][0], degC))

#add_metpy_logo(fig, 115, 100)

# Grid for plots
skew = SkewT(fig, rotation=45)
skew.ax.set_title(figtitle,fontsize=13)
skew.ax.set_ylabel('Pressure (hPa)')
skew.ax.set_xlabel('Temperature (\N{DEGREE CELSIUS})')

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p, T, 'r',linewidth=2.5)
skew.plot(p, Td, 'g',linewidth=2.5)

# Plot barbs
# Set spacing interval--Every 50 mb from 1000 to 100 mb
my_interval = np.arange(100, 1000, 40) * units('mbar')
# Get indexes of values closest to defined interval
ix = mpcalc.resample_nn_1d(p, my_interval)
# Plot only values nearest to defined interval values
skew.plot_barbs(p[ix], u[ix], v[ix])

# Plot height values at specified pressure levels on y axis
fig.text(0.18,0.883,'%.0f m'%(h_levels[0]),fontsize=8)
fig.text(0.18,0.65,'%.0f m'%(h_levels[1]),fontsize=8)
fig.text(0.18,0.515,'%.0f m'%(h_levels[2]),fontsize=8)
fig.text(0.18,0.42,'%.0f m'%(h_levels[3]),fontsize=8)
fig.text(0.18,0.345,'%.0f m'%(h_levels[4]),fontsize=8)
fig.text(0.18,0.285,'%.0f m'%(h_levels[5]),fontsize=8)
fig.text(0.18,0.23,'%.0f m'%(h_levels[6]),fontsize=8)
fig.text(0.18,0.187,'%.0f m'%(h_levels[7]),fontsize=8)
fig.text(0.18,0.147,'%.0f m'%(h_levels[8]),fontsize=8)
fig.text(0.18,0.113,'%.0f m'%(h_levels[9]),fontsize=8)

# Add the relevant special lines
da_temps = units.Quantity(np.linspace(-40,150,20), 'degC')
skew.plot_dry_adiabats(t0=da_temps,colors='red',alpha=0.3)
ma_temps = units.Quantity(np.array([-40, -30, -20, -10, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36]), 'degC')
ma_pres_levels = units.Quantity(np.linspace(1000.,220.), 'mbar')
skew.plot_moist_adiabats(t0=ma_temps,pressure=ma_pres_levels,colors='green',alpha=0.3)
skew.plot_mixing_lines(colors='purple',alpha=0.4)

# Add labels to moist adiabats
skew.ax.text(-86,215,'4',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-81,215,'8',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-75,215,'12',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-68,215,'16',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-58,215,'20',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-47,215,'24',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-35,215,'28',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-24,215,'32',color='green',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])

# Add labels to mixing lines
skew.ax.text(-33.5,600,'0.4',color='purple',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-24,600,'1',color='purple',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-15,600,'2',color='purple',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(-6.1,600,'4',color='purple',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(1,600,'7',color='purple',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(5.75,600,'10',color='purple',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])
skew.ax.text(11.5,600,'16',color='purple',fontsize=9,
             path_effects=[mpatheffects.withStroke(foreground='white', linewidth=3)])

# Good bounds
skew.ax.set_xlim(-40, 40)
#skew.ax.set_ylim(1000, 100)
skew.ax.set_ylim(1050, 100)

# Create a hodograph
if hodo:
    ax_hod = inset_axes(skew.ax, '25%', '25%', loc=1)
    h = Hodograph(ax_hod, component_range=80.)
    h.add_grid(increment=20)
    h.plot_colormapped(u, v, hght)

# Show the plot
plt.show()
