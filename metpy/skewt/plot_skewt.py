#!/usr/bin/python3

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

degC = '\N{DEGREE CELSIUS}'
hodo = False
outpath = '/tmp'
min_pres = 100.
max_pres = 1050.
min_temp = -40.
max_temp = 40.
barb_spacing = 25

if len(sys.argv) != 6:
    print('Usage: {} [input path] [input file] [output path] [fmt (Albany|CSU|MUnc|MUtxt_ws|NCSU|Purdue|SBUnc|SBUnc_mobile|UIUCnc|UNCA|UWYO|Valpo] [hodograph (True|False)]'.format(sys.argv[0]))
    sys.exit()
else:
    inpath = sys.argv[1]
    infile = sys.argv[2]
    outpath = sys.argv[3]
    fmt = sys.argv[4]
    hodo   = sys.argv[5]
    print('indir =',indir,'infile =',infile,'and fmt =',fmt)

# Test inputs
inpath = '/home/disk/bob/impacts/upperair/ualbany/20220225'
infile = 'upperair.sounding.202202250900.Albany-ESSX.txt'
fmt = 'Albany'

inpath = '/home/disk/monsoon/relampago/raw/sounding/CSU/FIELD/20181105'
infile = 'edt_20181105_0858.txt'
inpath = '/home/disk/monsoon/precip/raw/soundings/20220811'
infile = 'edt_20220811_1200.txt'
fmt = 'CSU'

inpath = '/home/disk/bob/impacts/upperair/umill/20220219'
infile = 'upperair.UMILL_sonde.20220219.nc'
fmt = 'MUnc'

inpath = '/home/disk/bob/impacts/upperair/umill/20220225'
infile = 'upperair.UMILL_windsonde1.202202250600.txt'
fmt = 'MUtxt_ws'

inpath = '/home/disk/bob/impacts/upperair/ncsu/20220213'
infile = ''
fmt = 'NCSU'

inpath = '/home/disk/bob/impacts/upperair/purdue/20220217'
infile = '2022-02-17_1509.sounding.csv'
fmt = 'Purdue'

inpath = '/home/disk/bob/impacts/upperair/rtso/20220213'
infile = '20220213.1500.rtso'
fmt = 'RTSO'

inpath = '/home/disk/bob/impacts/upperair/sbu/20220225'
infile = 'upperair.SBU_sonde_SBUSouthP.202202250525.nc'
fmt = 'SBUnc'

inpath = '/home/disk/bob/impacts/upperair/sbu/20220129'
infile = 'upperair.SBU_sonde_RadarTruck.202201290356.nc'
fmt = 'SBUnc_mobile'

inpath = '/home/disk/bob/impacts/upperair/uill/20220117'
infile = 'upperair.UILL_sonde.202201170600.nc'
fmt = 'UIUCnc'

inpath = '/home/disk/bob/impacts/upperair/unca/20220116'
infile = 'upperair.UNCA_sonde.202201160600.txt'
fmt = 'UNCA'

inpath = '/home/disk/bob/impacts/upperair/nws/20221018'
infile = 'upperair.SkewT.202210181200.ALB.new'
fmt = 'UWYO'

inpath = '/home/disk/bob/impacts/upperair/valpo/20220217'
infile = 'upperair.VALPO_sonde.202202172100.csv'
fmt = 'Valpo'

###########################################
# Read sounding data into dataframe with column names as indicated
# col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']

(df,out_fname,figtitle) = read_infile(inpath,infile,fmt)
if df.empty:
    sys.exit('Unrecognized fmt =',fmt)

###########################################
# We will pull the data out of the example dataset into individual variables and
# assign units.

hght = df['height'].values * units.m
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

# THIS CALL DOESN'T WORK FOR RTSO DATA - can't calculate pw
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
#skew.ax.grid(axis='x', color='tan')
skew.ax.grid(axis='x', color='darkgoldenrod')

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p, T, 'r',linewidth=2.5)
skew.plot(p, Td, 'g',linewidth=2.5)

# Plot barbs
# Set spacing interval--Every 50 mb from 1000 to 100 mb
my_interval = np.arange(min_pres,max_pres,barb_spacing) * units('mbar')
# Get indexes of values closest to defined interval
ix = mpcalc.resample_nn_1d(p, my_interval)
# Plot only values nearest to defined interval values
skew.plot_barbs(p[ix], u[ix], v[ix], linewidth=0.5)

# Plot height values at specified pressure levels on y axis
fig.text(0.18,0.883,'%.0f m'%(h_levels[0]),fontsize=8)
fig.text(0.18,0.655,'%.0f m'%(h_levels[1]),fontsize=8)
fig.text(0.18,0.525,'%.0f m'%(h_levels[2]),fontsize=8)
fig.text(0.18,0.43,'%.0f m'%(h_levels[3]),fontsize=8)
fig.text(0.18,0.355,'%.0f m'%(h_levels[4]),fontsize=8)
fig.text(0.18,0.295,'%.0f m'%(h_levels[5]),fontsize=8)
fig.text(0.18,0.244,'%.0f m'%(h_levels[6]),fontsize=8)
fig.text(0.18,0.20,'%.0f m'%(h_levels[7]),fontsize=8)
fig.text(0.18,0.16,'%.0f m'%(h_levels[8]),fontsize=8)
fig.text(0.18,0.125,'%.0f m'%(h_levels[9]),fontsize=8)

# Add the relevant special lines
da_temps = units.Quantity(np.linspace(-40,150,20), 'degC')
skew.plot_dry_adiabats(t0=da_temps,colors='red',alpha=0.3)
ma_temps = units.Quantity(np.concatenate((np.linspace(-40,0,5),np.linspace(4,36,9))),'degC')
ma_pres_levels = units.Quantity(np.linspace(max_pres,220.), 'mbar')
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
skew.ax.set_xlim(min_temp,max_temp)
skew.ax.set_ylim(max_pres,min_pres)

# Create a hodograph
if hodo:
    ax_hod = inset_axes(skew.ax, '25%', '25%', loc=1)
    h = Hodograph(ax_hod, component_range=80.)
    h.add_grid(increment=20)
    h.plot_colormapped(u, v, hght)

# Show the plot
plt.show()

# Save the plot
plt.savefig(outpath+'/'+out_fname)

# Ftp to field catalog
