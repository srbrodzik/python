# Brody Fuchs, April 2018
# plotting timeseries of rain rate (from PPI) and echo top height (from RHI)

import numpy as np 
import os
import argparse
import time
import Config
import glob
import datetime
import string_tools
from csuram import RadarData, RadarConfig, analysis_tools
import matplotlib.pyplot as plt 
from collections import OrderedDict
import matplotlib.dates as mdates
import cPickle as pickle 
import scipy.stats as stats
import plot_tools as ptools
import matplotlib.patches as mpatches
import matplotlib.dates as mdates


myFmt = mdates.DateFormatter('%H:%M')


base_path = os.path.dirname(os.path.realpath('__file__'))

# parse the command line arguments
parser = argparse.ArgumentParser(description='Put in a file to be processed')


parser.add_argument('--date', action="store", dest="date", default=None)


pargs = parser.parse_args()


dt_fmt = '%Y%m%d-%H%M%S'
prefix = 'SEA'
DPI = 150

data_path = '%s/radar_stats'%(base_path)


rr_files = sorted(glob.glob('%s/*%s*rr.log'%(data_path, pargs.date)))
dbzht_files = sorted(glob.glob('%s/*%s*dbzhts.log'%(data_path, pargs.date)))

colors = {0: 'gray', 10: 'blue', 20: 'green', 30: '#ffcc00', 40: 'orange', 50: 'red'}





for rrf in rr_files:
	rr = []
	rr_times = []
	dt_string = os.path.basename(rrf)[:8]
	date_obj = datetime.datetime.strptime(dt_string, '%Y%m%d')
	#	print date_obj
	print 'Generating rain rate figure for %s'%(dt_string)

	rrdata = np.genfromtxt(rrf, dtype=None, skip_header=0)
	rr_times.extend([datetime.datetime.strptime(_, dt_fmt) for _ in rrdata['f0']])
	rr.extend(rrdata['f1'])

	rr = np.array(rr)
	rr_times = np.array(rr_times)

	fig, ax = plt.subplots(1,1)

	ax.plot(rr_times, rr)
	ax.grid(True)

	ax.set_ylabel('Max RR (mm/hr)')
	ax.set_ylim(0, 200)
	ax.set_xlim(date_obj, date_obj+datetime.timedelta(days=1))

	fig.suptitle('SEAPOL max rain rate timeseries (mm/hr) for %s'%(dt_string), fontsize=14)
	ax.xaxis.set_major_formatter(myFmt)


	plt.tight_layout()
	fig.subplots_adjust(top=0.90)

	ptools.save_figure(fig, '%s/figures/radar_stats/%s'%(base_path, dt_string), '%s_maxrr_timeseries.png'%(prefix), dpi=DPI)







for dbf in dbzht_files:
	#dbzhts = OrderedDict(0: [], 10: [], 20: [], 30: [], 40: [], 50: [])
	dbzhts = OrderedDict( [(0, []), (10, []), (20, []), (30, []), (40, []), (50, [])] )
	dbz_times = []
	db_dt_string = os.path.basename(dbf)[:8]
	db_date_obj = datetime.datetime.strptime(db_dt_string, '%Y%m%d')


	print 'Generating echo top height figure for %s'%(db_dt_string)



	dbzdata = np.genfromtxt(dbf, dtype=None, skip_header=0)
	dbz_times.extend([datetime.datetime.strptime(_, dt_fmt) for _ in dbzdata['f0']])
	dbzhts[0].extend(dbzdata['f1'])
	dbzhts[10].extend(dbzdata['f2'])
	dbzhts[20].extend(dbzdata['f3'])
	dbzhts[30].extend(dbzdata['f4'])
	dbzhts[40].extend(dbzdata['f5'])
	dbzhts[50].extend(dbzdata['f6'])

	for _k in dbzhts.keys():
	    dbzhts[_k] = np.array(dbzhts[_k])


	fig, ax = plt.subplots(1,1)

	for k in dbzhts.keys():
	    dtb = ax.bar(dbz_times, dbzhts[k], facecolor=colors[k], edgecolor='none', width=0.003)

	ax.grid(True)
	ax.set_ylabel('Altitude (km MSL)')
	ax.set_ylim(0,15)
	ax.set_xlim(db_date_obj, db_date_obj+datetime.timedelta(days=1))
	#ax.legend(dtb)
	handles = []
	for idb in range(len(dbzhts.keys())):
		#print idb, colors[idb], dbzhts.keys()[idb], dbzhts.keys()[idb+1]
		handles.append(mpatches.Patch(color = colors[dbzhts.keys()[idb]], label = '%d dBZ'%(dbzhts.keys()[idb])))

    	#analysis_tools.multiple_axes_one_legend(ax, handles, prop = {'size': 10}, ncol = 4, loc = 'upper right', \
       	#bbox_to_anchor = (0.65, -0.15), fancybox = False)
	ax.legend(handles=handles, prop = {'size': 6}, fancybox = False, loc='upper center', ncol=6, bbox_to_anchor = (0.45, -0.075))


	ax.xaxis.set_major_formatter(myFmt)



	plt.tight_layout()

	fig.subplots_adjust(top=0.92, bottom=0.13)
	fig.suptitle('RHI-based echo top height timeseries for %s'%(db_dt_string))

	ptools.save_figure(fig, '%s/figures/radar_stats/%s'%(base_path, db_dt_string), '%s_dbzht_timeseries.png'%(prefix), dpi=DPI)


