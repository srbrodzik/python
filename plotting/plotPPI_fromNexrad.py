#!/usr/bin/python3

# Plot nexrad archived data using NWS method of displaying reflectivity

import pyart
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib import colors
#from colorspacious import cspace_converter
#from collections import OrderedDict

#----------------------------------------------------------------------------
# In case anyone ever wants to use a color scale with
# discrete color bins (mimics how the National Weather Service displays
# reflectivity), here it is:

# Segmented dBZ colormap
dbz_colors=[
    (0.39, 0.39, 0.39),
    (0.39, 0.39, 0.39),
    (0.81, 1, 1),
    (0.80, 0.61, 0.80),
    (0.60, 0.41, 0.60),
    (0.40, 0.21, 0.38),
    (0.80, 0.80, 0.61),
    (0.60, 0.60, 0.41),
    (0.39, 0.39, 0.39),
    (0.16, 0.91, 0.90),
    (0.09, 0.63, 0.95),
 
    (0.05, 0.13, 0.94),
    (0.16, 0.98, 0.18),
    (0.12, 0.76, 0.13),
    (0.07, 0.55, 0.08),
    (0.99, 0.96, 0.22),
    (0.89, 0.73, 0.16),
    (0.98, 0.58, 0.15),
    (0.98, 0.05, 0.11),
    (0.82, 0.04, 0.08),
    (0.73, 0.03, 0.07),
    (0.96, 0.16, 0.98),
    (0.59, 0.35, 0.77),
    (0.99, 0.99, 0.99),
    (0.99, 0.99, 0.99),
    (0.99, 0.99, 0.99),
    (0.99, 0.99, 0.99),
    (0.99, 0.99, 0.99)
]

n_dbz_colors = len(dbz_colors)

dbz_cbar_limits = [
    0.0, 0.03571, 0.07143, 0.1429, 0.1786, 0.2143, 0.25,
    0.2857, 0.3214, 0.3571, 0.3929, 0.4286, 0.4643, 0.5,
    0.5357, 0.5714, 0.6071, 0.6429, 0.6786, 0.7143, 0.75,
    0.7857, 0.8214, 0.8571, 0.8929, 0.9286, 0.9643, 1
]

n_dbz_limits = len(dbz_cbar_limits)

dbz_ticks = [
    -40, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
    50, 55, 60, 65, 70, 75, 100
]

dbz_ticklabs = [
    '-40', '-30', '-25', '-20', '-15', '-10', '-5', '0', '5', '10', '15',
    '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75',
    '100'
]

dbz_cmap = colors.LinearSegmentedColormap.from_list('test', dbz_colors,
                                                    n_dbz_colors)
#----------------------------------------------------------------------------

indir = '/home/disk/meso-home/brodzik/nexrad/KGJX'
infile = 'KGJX20210802_060045_V06.ar2v'
input_filename = indir+'/'+infile

radar = pyart.io.read_nexrad_archive(
    input_filename,
    include_fields=('reflectivity','velocity')
)
display = pyart.graph.RadarMapDisplay(radar)
dbz_figure = plt.figure(figsize=[12, 10], clear=True)
dbz_ax = dbz_figure.add_subplot(111)
display.plot_ppi(
    'reflectivity', sweep=1,
    ax=dbz_ax,
    title="Radar Base Reflectivity - Sweep 1",
    colorbar_flag=True,
    colorbar_label='dBZ',
    vmin=-40., vmax=100.,
    ticks=dbz_ticks, ticklabs=dbz_ticklabs,
    cmap=dbz_cmap,
    #cmap='pyart_Carbone17',
    ## embelish=True
)
plt.savefig(input_filename + '_BREF.png', dpi=200, transparent=False)

"""
supported library values for cmap are:
'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r',
'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r',
'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges',
'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r',
'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r',
'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu',
'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn',
'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3',
'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu',
'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot',
'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r',
'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r',
'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r',
'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r',
'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow',
'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r',
'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot',
'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma',
'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink',
'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'pyart_BlueBrown10',
'pyart_BlueBrown10_r', 'pyart_BlueBrown11', 'pyart_BlueBrown11_r',
'pyart_BrBu10', 'pyart_BrBu10_r', 'pyart_BrBu12', 'pyart_BrBu12_r',
'pyart_Bu10', 'pyart_Bu10_r', 'pyart_Bu7', 'pyart_Bu7_r', 'pyart_BuDOr12',
'pyart_BuDOr12_r', 'pyart_BuDOr18', 'pyart_BuDOr18_r', 'pyart_BuDRd12',
'pyart_BuDRd12_r', 'pyart_BuDRd18', 'pyart_BuDRd18_r', 'pyart_BuGr14',
'pyart_BuGr14_r', 'pyart_BuGy8', 'pyart_BuGy8_r', 'pyart_BuOr10',
'pyart_BuOr10_r', 'pyart_BuOr12', 'pyart_BuOr12_r', 'pyart_BuOr8',
'pyart_BuOr8_r', 'pyart_BuOrR14', 'pyart_BuOrR14_r', 'pyart_Carbone11',
'pyart_Carbone11_r', 'pyart_Carbone17', 'pyart_Carbone17_r', 'pyart_Carbone42',
'pyart_Carbone42_r', 'pyart_Cat12', 'pyart_Cat12_r', 'pyart_EWilson17',
'pyart_EWilson17_r', 'pyart_GrMg16', 'pyart_GrMg16_r', 'pyart_Gray5',
'pyart_Gray5_r', 'pyart_Gray9', 'pyart_Gray9_r', 'pyart_HomeyerRainbow',
'pyart_HomeyerRainbow_r', 'pyart_LangRainbow12', 'pyart_LangRainbow12_r',
'pyart_NWSRef', 'pyart_NWSRef_r', 'pyart_NWSVel', 'pyart_NWSVel_r',
'pyart_NWS_SPW', 'pyart_NWS_SPW_r', 'pyart_PD17', 'pyart_PD17_r',
'pyart_RRate11', 'pyart_RRate11_r', 'pyart_RdYlBu11b', 'pyart_RdYlBu11b_r',
'pyart_RefDiff', 'pyart_RefDiff_r', 'pyart_SCook18', 'pyart_SCook18_r',
'pyart_StepSeq25', 'pyart_StepSeq25_r', 'pyart_SymGray12', 'pyart_SymGray12_r',
'pyart_Theodore16', 'pyart_Theodore16_r', 'pyart_Wild25', 'pyart_Wild25_r',
'pyart_balance', 'pyart_balance_r', 'rainbow', 'rainbow_r', 'seismic',
'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r',
'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain',
'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted',
'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r'
"""
