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
There is an option to create a wet bulb temperature plot as well.

For testing:
./plot_skewt.py --inpath /home/disk/bob/impacts/upperair/sbu/20220517 --infile upperair.SBU_sonde.202205171809.RadarTruck.nc --outpath /tmp --fmt SBUnc_mobile --hodo false --wb false --vlim

"""

import os
import sys
sys.path.append('/home/disk/meso-home/brodzik/python/metpy/skewt')

import argparse
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
from ftplib import FTP
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from wetbulb import plot_wb
from wetbulb import calc_wetbulb_temp

KM2FEET = 3280.84
degC = '\N{DEGREE CELSIUS}'
min_pres = 100.
max_pres = 1050.
min_temp = -40.
max_temp = 40.
barb_spacing = 25
debug = False

# Parse input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inpath', action='store', dest='inpath', default='.')
parser.add_argument('--infile', action="store", dest="infile", required=True, default=None)
parser.add_argument('--outpath', action='store', dest='outpath', default='.')
parser.add_argument('--fmt', action="store", dest="fmt", required=True, default='', help='Albany|CMICH|CSU|HsinChu|MUnc|MUtxt_ws|NCSU|Oswego|Purdue|SBUnc|SBUnc_mobile|UIUCnc|UNCA|UND|UNDws|UWYO|Valpo')
parser.add_argument('--hodo', action='store', dest='hodo', default=False, help='True or False (default=False)')
parser.add_argument('--wb', action='store', dest='wb', default=False, help='True or False (default=False)')
parser.add_argument('--vlim', action='store', dest='vlim', default=7, help='vertical extent in km (default=7)')
parser.add_argument('--lat', action='store', dest='lat', default=0.0, help='site latitude (overridden if value in input file)')
parser.add_argument('--lon', action='store', dest='lon', default=0.0, help='site longitude (overridden if value in input file)')
pargs = parser.parse_args()

try:
    if pargs.hodo.lower() == 'false':
        pargs.hodo = False
    else:
        if debug:
            print('   Improper hodo flag. Use \"False\" to turn off the hodograph plot.')
        pargs.hodo = True
except:
    pass

try:
    if pargs.wb.lower() == 'false':
        pargs.wb = False
    else:
        if debug:
            print('   Improper wb flag. Use \"False\" to turn off the hodograph plot.')
        pargs.wb = True
except:
    pass

if debug:
    print('TOP OF plot_skewt:')
    print('   inpath  =',pargs.inpath)
    print('   infile  =',pargs.infile)
    print('   outpath =',pargs.outpath)
    print('   fmt     =',pargs.fmt)
    print('   hodo    =',pargs.hodo)
    print('   wb      =',pargs.wb)
    print('   vlim    =',pargs.vlim)
    print('   lat     =',pargs.lat)
    print('   lon     =',pargs.lon)

###########################################
# Read sounding data into dataframe with column names as indicated
# col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']

(df,out_fname,figtitle) = read_infile(pargs.inpath,pargs.infile,pargs.fmt,pargs.lat,pargs.lon)
if debug:
    print('FOR SKEWT:')
    print('   out_fname =',out_fname)
    print('   figtitle  =',figtitle)
if df.empty:
      sys.exit('Unable to create plot for file = '+pargs.inpath+'/'+pargs.infile)

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
if debug:
    print('Decide whether to add hodograph: hodo =',pargs.hodo)
if pargs.hodo:
    if debug:
        print('   Will create hodograph')
    ax_hod = inset_axes(skew.ax, '25%', '25%', loc=1)
    h = Hodograph(ax_hod, component_range=80.)
    h.add_grid(increment=20)
    h.plot_colormapped(u, v, hght)

# Show the plot
#plt.show()

# Save the plot
plt.savefig(pargs.outpath+'/'+out_fname)

# Plot wetbulb plot if desired
if pargs.wb:
    (base,ext) = os.path.splitext(out_fname)
    (category,platform,datetime,product) = base.split('.')
    if 'skewT' not in product:
        product = product+'_Wet_Bulb'
    else:
        product = product.replace('skewT','Wet_Bulb')
    wb_out_fname = category+'.'+platform+'.'+datetime+'.'+product+ext
        
    plot_wb(df, pargs.outpath, wb_out_fname, figtitle, pargs.vlim)


