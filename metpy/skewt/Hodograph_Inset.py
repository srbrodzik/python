# Copyright (c) 2016 MetPy Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause

"""
===============
Hodograph Inset
===============

Layout a Skew-T plot with a hodograph inset into the plot.
"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, Hodograph, SkewT
from metpy.units import units
from calculations import calcs

degC = '$^{\circ}$C'

###########################################
# Upper air data can be obtained using the siphon package, but for this example we will use
# some of MetPy's sample data.

col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']

#df = pd.read_fwf(get_test_data('may4_sounding.txt', as_file_obj=False),
#                 skiprows=5, usecols=[0, 1, 2, 3, 6, 7], names=col_names)
df = pd.read_fwf('input/may4_sounding.txt',
                 skiprows=5, usecols=[0, 1, 2, 3, 6, 7], names=col_names)

# Drop any rows with all NaN values for T, Td, winds
df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'),
               how='all').reset_index(drop=True)

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

# Create a new figure. The dimensions here give a good aspect ratio
fig = plt.figure(figsize=(10, 9))
#fig.suptitle('Sample SkewT', fontsize=16)
fig.text(0.4,0.9,'MetPy SkewT', fontsize=18)
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
fig.text(0.86,0.20,'Surface')
fig.text(0.86,0.18,'P: %.1f hPa'%df['pressure'][0])
fig.text(0.86,0.16,'H: %.0f m'%(df['height'][0]))
fig.text(0.86,0.14,'T: %.1f %s'%(df['temperature'][0], degC))
fig.text(0.86,0.12,'T$_{D}$: %.1f %s'%(df['dewpoint'][0], degC))

#add_metpy_logo(fig, 115, 100)

# Grid for plots
skew = SkewT(fig, rotation=45)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p, T, 'r')
skew.plot(p, Td, 'g')
skew.plot_barbs(p, u, v)
skew.ax.set_ylim(1000, 100)

# Add the relevant special lines
skew.plot_dry_adiabats(colors='black',alpha=0.3)
skew.plot_moist_adiabats(colors='green',alpha=0.3)
skew.plot_mixing_lines(colors='purple',alpha=0.3)

# Good bounds for aspect ratio
skew.ax.set_xlim(-40, 40)

# Create a hodograph
ax_hod = inset_axes(skew.ax, '25%', '25%', loc=1)
h = Hodograph(ax_hod, component_range=80.)
h.add_grid(increment=20)
h.plot_colormapped(u, v, hght)

# Show the plot
plt.show()
