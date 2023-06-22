#!/net/denali/home/brfuchs/epd-7.2-2-rh3-x86_64/
# Brody Fuchs, CSU, May 2015
# brfuchs@atmos.colostate.edu

# Code to now realtime plot radar and lightning
# with NMQ output grib files, the damn format keeps
# changing!! Good thing here is that flash counting
# is already done in another realtime code so we just
# need to read in that data

from __future__ import division

import os
import sys
import numpy as np
import datetime 
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from time import gmtime, strftime, sleep
from matplotlib import colors
import rt_tools
import glob
import time
from cmap import make_cmap
from scipy.ndimage.measurements import label
import mmmpy # this version just skips the conversion to netcdf and reads from binary, saves time!!
import matplotlib.patheffects as PE
import lightning_tools as LT
import gzip

#import cartopy.crs as ccrs


max_ht = 3.0 # km MSL
lma_flag = 1
cbar_flag = 0
lobe_flag = 1
n_images = 12
thresh1 = 35. # dBZ threshold for cells
min_area_pts = 15. # minimum number of points to be called a cell
sdist = 15. # maximum a flash can be from a cell to be counted in km
re = 6374.0 # earth radius


den_cvals = pop_color_vals = [(255,255,255), (0,200,0), (0,0,255), (205,20,127), (255,0,0), \
		(255,128,0), (255,208,0), (255,255,255)] 
# blue, purple, red, orange, yellow, white
#den_cvals = [(255,255,255), (255,0,0)]
den_cmap = make_cmap(den_cvals, bit = True)

nc_dir='/var/home/brfuchs/public_html/realtime/ncdummy'
radar_dir = '/net/andes/data3/mosaic3d/ldm/data/mrefl_ops'
f_dir = '/var/home/brfuchs/public_html/realtime/temp_flash_count'
cur_dir = '/var/home/brfuchs/public_html/realtime/current'
case_dir = '/var/home/brfuchs/public_html/realtime/cases'
base_dir = '/var/home/brfuchs/public_html'
dz_thresh = -1.0
rt_dir = '/var/home/brfuchs/public_html/realtime'
l_dir = '/data2/colma/data'
n_flash = 0
flash_flag = 1
src_alpha = 1
flash_alpha = 1
size = 20
f_size = 20
delta = 6.0
cell_flag = 0
radar_flag = 0 # if radar data actually available, set to 0 first

horizontal_res = 2.0

rt_flag = True

date_fmt = '%Y%m%d-%H%M%S'

# get cell number otherwise exit
if len(sys.argv) > 1:
    for i, arg in enumerate(sys.argv):
        if arg == 'time':
            date_str = sys.argv[i+1]
	    radar_time = datetime.datetime.strptime(date_str, date_fmt)
	    rt_flag = False
else:
   # this is realtime stuff

    radar_time = datetime.datetime.utcnow()



print 'realtime flag: {}'.format(rt_flag)

print radar_time
year = radar_time.year
month = radar_time.month
date = radar_time.day
hour = radar_time.hour
minute = radar_time.minute


dummy = np.fmod(minute, delta)
minute_f_s = np.int(minute - dummy)
print minute_f_s


date_str = '%02d%02d%02d'%(year-2000,month,date)


# start of stitching, cropping file
    
params_file = open('%s/latlon.txt'%rt_dir).readlines()
params = [np.float(i.split('\t')[1]) for i in params_file]
latmin, latmax, lonmin, lonmax = params[0], params[1], params[2], params[3]

# histogram bins
latbins = np.arange(latmin, latmax, 0.01)
lonbins = np.arange(lonmin, lonmax, 0.01)
altbins = np.arange(0, 15, 0.2)

# NOW HAVE CORRECT TILES
# MAKING FINAL ARRAYS TO BE FED DATA

dznewlat = latmax-0.01*np.arange((latmax-latmin)*100+1)
dznewlon = lonmin+0.01*np.arange((lonmax-lonmin)*100+1)

radarcbar=['PeachPuff','Aqua','DodgerBlue','MediumBlue','Lime', \
		'LimeGreen','Green','Yellow','Gold','Orange','OrangeRed', \
		'Red','Crimson','Fuchsia','Indigo','DarkCyan','White']

os.system('which python') # wanna check which python is being used here
os.system('echo $PATH')
lowdbz=5
highdbz=80
deldbz=5
rangedbz=(highdbz-lowdbz)
radar_cmap=colors.ListedColormap(radarcbar)
radar_bounds=np.linspace(lowdbz,highdbz,2+rangedbz/deldbz)
radar_norm=colors.BoundaryNorm(radar_bounds, radar_cmap.N)
# getting current time

# NOW READ IN THE TILES WITH THE NEW WAY
filename = '%s/%d%02d%02d/MergedRefl_%d%02d%02d_%02d%02d*.grib2.gz'\
				%(radar_dir, year, month, date, year, month, date, hour, minute_f_s)



print filename
os.system('rm %s/*.gz'%nc_dir)




done = False
for k in range(10):
    if done == False:
	f = sorted(glob.glob(filename))[::2]
	#print f

	if len(f) >= 5: # files are there

	    print 'FOUND FILES'
	# now loop thru the files
	    for fi in range(len(f)):
		f_base = os.path.basename(f[fi])
		grib2_loc = f_base.find('grib2')
		this_ht = np.float(f_base[grib2_loc-6: grib2_loc-1])
		if this_ht <= max_ht:

	    	    os.system('cp %s %s'%(f[fi], nc_dir))
#	    os.system('gunzip %s/%s.gz'%(nc_dir, filename))

	    	    tile1 = mmmpy.MosaicTile()
	    	    tile1.read_mosaic_grib('%s/%s'%(nc_dir, os.path.basename(f[fi])), wgrib2_path='/usr/local64/bin/', \
			keep_nc=False, nc_path='./ncdummy', wgrib2_name='wgrib2', verbose=False, \
			latrange=[latmin, latmax], lonrange=[lonmin, lonmax])
		# the new way to read in the grib files that are output, pretty fast since you can pre-subset
		# the data with the latrange and lonrange before using wgrib2 command

		    if fi is 0: dz = tile1.mrefl3d # make array if first time thru loop
		    else: dz = np.vstack((dz, tile1.mrefl3d)) # else stack on top of existing array

		else: break # if we're above the max_ht, break out of the reading-in-files loop
	    
	    done = True

	else:
	    print 'STILL HAVENT FOUND FILES'
	    time.sleep(10)	
#dz = tile1.mrefl3d.squeeze() # getting the reflectivity data

h_mult = 1

try:
    cr = np.max(dz, axis = 0)[::h_mult, ::h_mult]
    radar_flag = 1
    dlat = tile1.LatGridSpacing*h_mult
    dlon = tile1.LonGridSpacing*h_mult
except NameError:
    print 'FILES ARE NOT FOUND FOR SOME REASON, STUPID NMQ, STILL GOING AHEAD WITH PLOTTING LMA DATA'
    dlat, dlon = 0.01*h_mult, 0.01*h_mult # if no radar, just assume dlat and dlon are 0.01, safe assumption

print 'RADAR FLAG: %d'%radar_flag

#cell_flag = 1 # demanding this to be 0 so don't contour and such

#f_dir gets to /current

directory='%s/'%f_dir

latcen = (latmin+latmax)/2.0
loncen = (lonmin+lonmax)/2.0


if lma_flag:


    # lma_files are files to be copied over and read in
    # ISSUES WITH INCREMENTING MINUTES BACKWARD AT THE BEGINNING OF THE HOUR
    lt = [datetime.datetime(year,month,date,hour,minute_f_s) - datetime.timedelta(minutes = k) for k in range(np.int(delta))] 


		# lma file times

    lma_files = np.array(['%s/%d%02d%02d_%02d%02d00_0060_flashes.txt' \
		%(f_dir,k.year,k.month,k.day,\
		k.hour,k.minute) for k in lt])


    print 'possible lma files: {}'.format(lma_files)

    src_files = np.array(['%s/%d%02d%02d_%02d%02d00_0060_sources.txt' \
		%(f_dir,k.year,k.month,k.day,\
		k.hour,k.minute) for k in lt])


    lma_check = np.array(sorted(glob.glob('%s/*%s*flashes*.txt'%(f_dir,date_str))))
    #print lma_check
    good_files = np.intersect1d(lma_files, lma_check) # files that are there to grab

    src_check = np.array(sorted(glob.glob('%s/*%s*sources*.txt'%(f_dir,date_str))))
    good_src_files = np.intersect1d(src_files, src_check)
    print good_src_files

    end_log_file = '%s/logfiles/%s.log'%(rt_dir, os.path.basename(good_files[-1])[11:13])
    logfilesize = os.stat(end_log_file).st_size
    print 'log file size: {}'.format(logfilesize) # if last LMA file is not finished, log file will be 0 bytes, so wait

    if (good_files.shape[0] < delta) or (logfilesize == 0): # if files not found, wait 1 min then check again
        print 'WAITING FOR LMA FILES'
	print 'LMA FILES PRESENT: ', good_files
        sleep(40) # wait 90 seconds
    #    lma_check = np.array(sorted(glob.glob('%s/%s/*/LYLOUT*.dat.gz'%(l_dir,date_str)))) # checking for files again
    	lma_check = np.array(sorted(glob.glob('%s/*%s*flashes*.txt'%(f_dir,date_str))))
        good_files = np.intersect1d(lma_files, lma_check) # files that are there to grab

    n_files = good_files.shape[0]
    print 'delta: {}, nfiles: {}'.format(delta, n_files)

    flash = {'lat': np.zeros(0), 'lon': np.zeros(0), 'alt': np.zeros(0), 'area': np.zeros(0)}
    src = {'lat': np.zeros(0), 'lon': np.zeros(0), 'alt': np.zeros(0)}
    
	# now proceed to read in the flash files
    for ifi in range(good_files.shape[0]):
	try:
	    blah = np.genfromtxt(good_files[ifi], skip_header = 1)
	    temp = {'lat': blah[:, 1], 'lon': blah[:, 2], 'alt': blah[:, 3], 'area': blah[:, 7]}
	    flash = LT.dict_concatenate(flash, temp) 
	    blah_src = np.genfromtxt(good_src_files[ifi], skip_header = 1)
	    temp_src = {'lat': blah_src[:,1] , 'lon': blah_src[:,2], 'alt': blah_src[:,3]}
	    src = LT.dict_concatenate(src, temp_src) 


	except Exception:
	    pass

    print 'number of flashes: {}'.format(flash['lat'].shape[0])


else:
    lma_time = np.zeros(0)
    n_flash = 0


################## START OF JUST ONE PLOT WITH JUST RADAR ###################
figscale = 1.1

DPI=72
figa = plt.figure(figsize=(figscale*820.0/DPI, figscale*760.0/DPI),dpi=DPI) # figure size and creation
ax_rad = plt.subplot2grid((6,6), (0,0), colspan=5, rowspan=5)
ax_lat = plt.subplot2grid((6,6), (5,0), colspan=5, rowspan=1)
ax_lat.set_ylim(0,15)
ax_lat.set_xlim(lonmin, lonmax)

ax_lon = plt.subplot2grid((6,6), (0,5), colspan=1, rowspan=5)
ax_lon.set_xlim(0,15)
ax_lon.set_ylim(latmin, latmax)

ax_hist = plt.subplot2grid((6,6), (5,5), colspan=1, rowspan=1)
ax_hist.set_ylim(0,15)

plt.subplots_adjust(bottom=0.03, top=0.93, left=0.05, right=0.95, hspace=0.14, wspace=0.14)



m = Basemap(projection='merc',llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax,resolution='i', ax = ax_rad)
#    map.bluemarble(scale=1)    
m.drawstates(color='k',linewidth=1.0)   
m.drawrivers(color='b',linewidth=0.5)
m.drawcoastlines(color='k', linewidth=1.0)
m.drawcounties(color='0.6', linewidth=0.5)
m.drawmapscale(lonmax-0.50,latmin+0.25,loncen,latcen,100,units='km')
m.drawmapscale(lonmin+0.2,latmin+0.25,loncen,latcen,25,units='km')
m.drawmeridians(np.arange(np.floor(lonmin), np.ceil(lonmax), 0.5), labels=[False, False, False, True])
m.drawparallels(np.arange(np.floor(latmin), np.ceil(latmax), 0.5), labels=[True, False, False, False])


#o,a = np.meshgrid(tile1.Longitude,tile1.Latitude)
#x,y=map(o,a)



# do radar stuff if we have it
if radar_flag:
    bad = cr < 0.
    cr[bad] = 0.

    cr_dummy = np.copy(cr)
    dummy_bad = cr_dummy < thresh1
    cr_dummy[dummy_bad] = 0.
    cell1, nfeat = label(cr_dummy) # doing cell identification
    print 'number of original features: {}'.format(nfeat)

    nlat = dz.shape[1]
    nlon = dz.shape[2]

    if nfeat > 1:
        cell_flag = 1
        uniq = np.unique(cell1)[1:] # leave out first one cuz it's 0

        for icell in uniq: # loop thru each cell to see size
            this_cell = np.where(cell1 == icell)
            n_pts = this_cell[0].shape[0]
    #        print icell, n_pts
            if n_pts < min_area_pts: # if cell too small, set cell value back to 0
                cell1[this_cell] = 0

        sub1 = np.unique(cell1)[1:]
        print 'features after filtering: {}'.format(sub1)
    x, y = m(tile1.Longitude,tile1.Latitude)

    # PLOTTING RADAR CONTOUR FILLS
    dummy1=ax_rad.contourf(x,y,cr,levels=np.linspace(lowdbz,highdbz,1+rangedbz/deldbz),cmap=radar_cmap,norm=radar_norm,alpha=0.9)
    # IF CELLS ARE FOUND, PLOT THEM CONTOUR STYLE
    if cell_flag == 1:
        try:
            ax_rad.contour(x, y, cell1, levels = [1, np.inf], colors = 'black', linewidths = 1.2)
        except ValueError:
            pass


#cities=['Ft Collins','Denver','Limon','Sterling','Cheyenne','Ft Morgan','Estes Park','Idaho Spgs','Greeley']
#citylons=[-105.075,-105,-103.67,-103.2,-104.82,-103.8,-105.55,-105.51,-104.73]
#citylats=[40.56,39.75,39.25,40.62,41.15,40.25,40.37,39.74,40.42]
cit = open('%s/city_latlon.txt'%rt_dir).readlines()
cities = []
#cities = np.array([i.split()[3].split(',')[0] for i in cit if (len(i.split())==4) or i.split()[3].split(',')[0]+' '+i.split()[4].split(',')[0] if (len(i.split())>4)])
for i in cit:
    length = len(i.split())
    if length == 4:
	cities.append(i.split()[3].split(',')[0])
    else:
	cities.append(i.split()[3].split(',')[0]+' '+i.split()[4].split(',')[0])

cities = np.array(cities)
#cities = np.array([i.split()[3].split(',')[0] for i in cit if len(i]) 
citylats = np.array([float(i.split()[1]) for i in cit])
citylons = np.array([-1.0*float(i.split()[2]) for i in cit])
g_lat = rt_tools.between(citylats, latmin+10*dlat, latmax-10*dlat) 
g_lon = rt_tools.between(citylons, lonmin+10*dlon, lonmax-40*dlon) 
# wont plot cities within 40 km of bdry, text overflow
g = np.intersect1d(g_lat[0], g_lon[0]) # where both lon and lat in range
cities = cities[g]
citylats = citylats[g]
citylons = citylons[g]
mtlat=[40.25,40.54]
mtlon=[-105.61,-105.2]

#berklat = 40.5810802
#berklon = -105.1298763
#xberk, yberk = m(berklon, berklat)
#ax_rad.text(xberk, yberk, 'B', horizontalalignment = 'center', \
#	verticalalignment = 'center', color = 'DarkOrange')
#ax_rad.scatter(xberk-400 ,yberk, s = 140, facecolors = 'none') 

# now go grab the first LMA file to check which stations are active and plot
# active stations in green and non-active in red
# start by making station dictionary
# maybe just read from lma status page? then can color yellow 
# if not active for more than a day

print 'lma_flag: %d'%lma_flag
stations = {}
# station name will be dictionary keys
#if lma_flag == 1:
#    print 'files: '%lma_files
try:

    last_lma_file = sorted(glob.glob('%s/%02d%02d%02d/%02d/LYLOUT*'%(l_dir, year-2000, month, date, hour)))[-1]
    stat = gzip.open(last_lma_file).readlines() # file is zipped so can use gzip.open to look at it normally w/o unzipping
    for line in stat:
        if line.find('Number of active stations:') != -1:
            n_stn = int(line.split()[-1])
        if line.find('Active stations:') != -1:
            act_stat = line.split()[2:]
            break

    print 'n_stn: %d'%n_stn
    print 'active stations: ', act_stat

    stat_list = open('%s/stations.txt'%rt_dir).readlines()
    s_id = [i.split()[0] for i in stat_list] # station ID letters to be plotted
    s_lat = [float(i.split()[2]) for i in stat_list]
    s_lon = [float(i.split()[3]) for i in stat_list] #station lats/lons to be plotted
    xst, yst = m(s_lon, s_lat)

    #good_stat = [i for i in s_id if i in act_stat]


    for il, line in enumerate(stat_list):
        if line.split()[0] in act_stat:
            m.scatter(xst[il], yst[il], s = 100, c = 'LimeGreen')
        else:
            m.scatter(xst[il], yst[il], s = 100, c = 'red')
except (AttributeError,IndexError, IOError):
    stat = open('%s/stations.txt'%rt_dir).readlines()
    s_lat = [float(i.split()[2]) for i in stat]
    s_lon = [float(i.split()[3]) for i in stat] #station lats/lons to be plotted
    xst,yst=m(s_lon,s_lat)
    m.scatter(xst,yst, s=100, facecolors = 'none', edgecolors = 'black') # plotting stations    
    n_stn = 0
    print 'NO LMA FILES?'



xmt,ymt=m(mtlon,mtlat)
#x1,y1=m(g_lons,g_lats)
x1,y1 = m(citylons, citylats)
if cbar_flag == 1:
    try:
    	cb = plt.colorbar(dummy1,ticks=np.linspace(lowdbz,highdbz,1+rangedbz/deldbz),shrink=0.7,pad=0.03, fraction = 0.07)
    	cb.set_label('dBZ')
    except NameError:
    	pass
m.plot(x1,y1,'k.')
m.plot(xmt,ymt,'r^')
#    for name,xpt,ypt in zip(marker,xmark,ymark):
#        m.plot(xpt,ypt,name)
#for name,xpt,ypt in zip(g_cities,x1,y1):
for name,xpt,ypt in zip(cities,x1,y1):

    ax_rad.text(xpt,ypt,name)
# plotting DD lobes

#lobe_lats = [40.591, 40.726, 40.0367, 40.196, 40.6526, 40.9453, 40.9326, 41.0903,40.445] # last one is the lma
#lobe_lons = [-105.044, -104.308, -105.162, -104.020, -105.333, -104.110, -105.004, -104.517,-104.639]
#lobe_rads = [47.7, 47.7, 73.7, 73.7, 80.0, 80.0, 32.2, 32.2, 150.0]
#l_colors = ['black', 'black', 'orange', 'orange', 'green', 'green', 'red', 'red', 'yellow']

#lobe_lats = [40.0367, 40.196, 40.6526, 40.9453, 40.445] # last one is the lma
#lobe_lons = [-105.162, -104.020, -105.333, -104.110, -104.639]
#lobe_rads = [73.7, 73.7, 80.0, 80.0, 100.0]
#lobe_colors = ['orange', 'orange', 'green', 'green', 'fuchsia']

lobe_lats = [40.50, 40.50] # last one is the lma
lobe_lons = [-104.4, -104.4]
lobe_rads = [100.0, 200.0]
lobe_colors = ['fuchsia', 'fuchsia']



if lobe_flag == 1:
#    print 'plotting the DD and other lobes of interest'
#    for i in range(len(lobe_lats)): # plotting circles, DD lobes and LMA range
#	print 'lobe %d'%i
#    	rt_tools.equi(m, lobe_lons[i], lobe_lats[i], lobe_rads[i]/111.0, lw=1.2,color = l_colors[i])#,linestyle = '--')

    radius_deg = 1. / 96.



#    lats = [40.0367, 40.196, 40.6526, 40.9453]
#    lons = 360+np.array([-105.162, -104.02, -105.333, -104.11])

    #lats = [40.196, 40.9453]
    #lons = np.array([-104.02, -104.11])
    #radii = [73.7, 73.7, 80, 80]

    #color = ['#a65628', '#a65628', '#ff7f00', '#ff7f00']

    for lo, la, ra, cl in zip(lobe_lons, lobe_lats, lobe_rads, lobe_colors):

	m.tissot(lo, la, ra * radius_deg, 99, fill=False, linestyle='solid', edgecolor=cl, linewidth=2, zorder=100, alpha=0.8)

#for stn in plt_stns.keys():
#
#    x, y = m(plt_stns[stn][0], plt_stns[stn][1])
#
#    m.scatter(x, y, s=9, color='w', marker='o', edgecolor='none', zorder=1000)



#### THIS IS WHERE FLASH ATTRIBUTION WOULD GO

if cell_flag and flash_flag and radar_flag: # if cells actually present, attribute flashes to them
    flash['feat'] = LT.cell_attribute(flash['lat'], flash['lon'], cell1, latmax, \
                                lonmin, latcen, loncen, nlat, nlon, latdelta = dlat, londelta = dlon, \
                                latincreasing = 0, sdist = 10)

    print 'flash features'
    print 'unique flash features : {}'.format(np.unique(flash['feat']))

    #n_flash_cell = np.zeros(sub1[1:].shape[0]) # dont want '0' cuz that's where no cell is
    n_flash_cell = []
    for s in sub1: # looping thru each cell to find out how many flashes
	print 'sub1 looped cell: {}'.format(s)
        this_flash_cell = np.where(flash['feat'] == s)
       	#print 'this_flash_cell: ', this_flash_cell
    	n_flash_cell.append(np.where(flash['feat'] == s)[0].shape[0])
        this_cell_loc = np.where(cell1 == s)

#       print 'this_cell_loc: ', this_cell_loc
        cell_avg_lat_index = int(np.round(np.average(this_cell_loc[0])))
        cell_avg_lon_index = int(np.round(np.average(this_cell_loc[1])))

        cell_lat = tile1.Latitude[cell_avg_lat_index][0] # tile1.lat and lon are 2d arrays so have to make 1D
        cell_lon = tile1.Longitude[0][cell_avg_lon_index] + 0.07 # shift it a little east of center

    #    print cell_lat, cell_lon, n_flash_cell
        avg_x, avg_y = m(cell_lon, cell_lat)
#       ax_rad.text(avg_x, avg_y, sub1[s]) 
	# print the newest one each time cuz that is getting appended for each cell
	if n_flash_cell[-1] > n_files:
            ax_rad.text(avg_x, avg_y, np.int(np.ceil(n_flash_cell[-1]/good_files.shape[0])), fontsize = 11, \
                color = 'black', horizontalalignment = 'center', path_effects = [PE.withStroke(linewidth = 2.5, foreground = 'w')])
#    print f_feat
    print 'number of flashes in cell: {}'.format(n_flash_cell)


#quit()

#### END OF FLASH ATTRIBUTION

den_norm = 35

# lma_* are the sources that are going to be plotted
if lma_flag == 1:
#    lma_x, lma_y = m(lma_lon, lma_lat)
#    m.scatter(lma_x, lma_y, marker = '.', s=10, c='k', alpha = 0.12) # plot sources very transparent
     # compute histograms of sources to get at densities
    lat_alt_hist, lonedges, altedges = np.histogram2d(src['lon'], src['alt']/1000.0, bins = (lonbins, altbins))
    print 'max of LMA lat src density : {}'.format(lat_alt_hist.max())

    #lat_alt_hist /= 1.2*lat_alt_hist.max() # basically normalizing the values
    if lat_alt_hist.max() > 0: # actually any sources?? then plot
    	lat_alt_norm = np.max([lat_alt_hist.max(), den_norm])
	lat_alt_hist = lat_alt_hist.astype(float)/float(lat_alt_norm)
        ax_lat.pcolormesh(lonedges, altedges, lat_alt_hist.T, cmap=den_cmap, vmin=0, vmax=1)
        ax_lat.set_xlim(lonmin, lonmax)

	


        lon_alt_hist, altedges, latedges = np.histogram2d(src['alt']/1000.0, src['lat'], bins = (altbins, latbins))

    	#lat_alt_hist /= 1.2*lat_alt_hist.max()
    	print 'max of LMA lon src density : {}'.format(lon_alt_hist.max())
    	lon_alt_norm = np.max([lon_alt_hist.max(), den_norm])
	lon_alt_hist = lon_alt_hist.astype(float)/float(lon_alt_norm)
        ax_lon.pcolormesh(altedges, latedges, lon_alt_hist.T, cmap=den_cmap, vmin=0, vmax=1) # vmax fn of resolution

        alt_hist, edges1d = np.histogram(src['alt']/1000.0, bins = altbins)
	print 'plotting ax_hist'
	alt_hist = alt_hist.astype(float)/float(alt_hist.sum())
        ax_hist.plot(alt_hist, altbins[:-1], 'black', linewidth = 2)    
	ax_hist.set_xlim(0, 0.1)
	ax_hist.set_xticks([])

    ax_lon.set_ylim(latmin, latmax)
    ax_lon.yaxis.tick_right()
    ax_lon.xaxis.tick_top()

    plt.setp(ax_hist.get_xticklabels(), visible=False)
#    plt.setp(ax_hist.get_yticklabels(), visible=False)
    ax_hist.yaxis.tick_right()

    line_alpha = 0.25

    ax_lon.axvline(x=4, color='black', linewidth=2, alpha=line_alpha)
    ax_lon.axvline(x=8, color='black', linewidth=2, alpha=line_alpha)
    ax_lon.axvline(x=12, color='black', linewidth=2, alpha=line_alpha)

    ax_lat.spines['top'].set_visible(False)

    ax_lat.axhline(y=4, color='black', linewidth=2, alpha=line_alpha)
    ax_lat.axhline(y=8, color='black', linewidth=2, alpha=line_alpha)
    ax_lat.axhline(y=12, color='black', linewidth=2, alpha=line_alpha)




    #print 'lat_alt_hist max', lat_alt_hist.max()
    #print 'flash areas', 90+2.0*flash['area']

    if flash_flag == 1:
	print 'PLOTTING FLASHES'
    	flash_x, flash_y = m(flash['lon'],flash['lat'])
    	m.scatter(flash_x, flash_y, marker= '.', s= 60+3.0*flash['area'], alpha=0.35, edgecolors='black', facecolors='none')
    	n_flash = np.float(flash['lat'].shape[0])*delta/n_files
	n_src = np.float(src['lat'].shape[0])*delta/n_files
	#print flash['lon']

    else: 
    	n_src = 0
	n_files = good_files.shape[0]

#n_src = 0

plt.suptitle('Composite dBZ at %d/%02d/%02d - %02d:%02d, %d LMA sources, %d LMA flashes\n%d active stations, %d LMA files\t Storm flash rates are min$^{-1}$'
	%(year,month,date,hour,minute_f_s,n_src,n_flash,n_stn,n_files), fontsize = 16)

if rt_flag:
    plt.savefig('%s/%d_%02d%02d_%02d%02d_%02d_%02d_%03d_%03d.png'%(cur_dir,year,month,date,
	hour,minute_f_s,np.abs(latmin),np.abs(latmax),np.abs(lonmin),np.abs(lonmax)),dpi=DPI)
    os.system('pngcrush -ow -q -reduce -brute %s/%d_%02d%02d_%02d%02d_%02d_%02d_%03d_%03d.png'%(cur_dir,year,month,date,
        hour,minute_f_s,np.abs(latmin),np.abs(latmax),np.abs(lonmin),np.abs(lonmax)))



    ######### BEGIN NEW WAY TO FILE MANIPUATE ##############

    # start by moving oldest file in the directory
    # to appropriate archive directory
    # need to read only the most recent file in
    # make plot then write png file to realtime directory
    # then go into radarmet and delete radar.txt file with all pngs to plot
    # and recreate it with new list of files
    cur_files = sorted(glob.glob('%s/*.png'%cur_dir), key = os.path.getmtime) # sorted by time, oldest first

    print len(cur_files)
    latest_file = cur_files[-1]
    if len(cur_files) >= n_images+1: # only move files from current directory if there are enough in there already

	for im in range(0, len(cur_files) - n_images): # loop thru the files that are more than n_images and move them

	    tf = cur_files[im]
	    tf_base = os.path.basename(tf)
	    tf_year = tf_base[0:4]
	    tf_day = tf_base[5:9]
	    
	    os.system('mv %s %s/%s/%s%s'%(tf, rt_dir, tf_year, tf_year, tf_day))

    new_files = sorted(glob.glob('%s/*.png'%cur_dir), key = os.path.getmtime) # sorted by time, oldest first
    files_string = ' '.join(new_files)

    # convert valid files into an animated gif for phones that dont have flash
    #os.system('convert -delay 50 -loop 0 %s /var/home/brfuchs/public_html/anim.gif'%files_string)
    # copy the most recent file to latest.png for still images
    os.system('cp %s /var/home/brfuchs/public_html/latest.png'%(latest_file))
    # NOW GET FILES IN CURRENT DIRECTORY TO MAKE A RTRADAR.TXT FILE TO BE READ BY FLANIS
    #txt_files = sorted(glob.glob('%s/*.png'%f_dir), key = os.path.getmtime)
    new_bases = [os.path.basename(_) for _ in new_files]
    outfile = open('/var/home/brfuchs/public_html/rtradar.txt','w')
    for ifi in new_bases:
	outfile.write('./realtime/current/%s\n'%ifi)

    # UPDATE CURRENT DIRECTORY'S index.html for archived viewing
    os.system("make_index %s/current 'Current Images'"%(rt_dir))
    # UPDATE TODAY'S DIRECTORY to add new archived image
    os.system("make_index %s/%d/%d%02d%02d '%d%02d%02d'"%(rt_dir,year,year,month,date,year,month,date))

    # ************ ALSO MAKE AN ANIMATED GIF FOR MOBILE VIEWING *************
    os.system('convert -delay 50 -loop 0 %s/current/*.png %s/anim.gif'%(rt_dir, base_dir))



else:
    plt.savefig('%s/%d_%02d%02d_%02d%02d_%02d_%02d_%03d_%03d.png'%(case_dir,year,month,date,
	hour,minute_f_s,np.abs(latmin),np.abs(latmax),np.abs(lonmin),np.abs(lonmax)),dpi=DPI)


print datetime.datetime.utcnow() - radar_time
