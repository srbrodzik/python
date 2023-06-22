# Brody Fuchs

# trying to get echo top heights from RHIs

import numpy as np 
import os
import pyart
import argparse
import time
import Config
import glob
import datetime
import string_tools
from csuram import RadarData, RadarConfig
import matplotlib.pyplot as plt 
from collections import OrderedDict
import matplotlib.dates as mdates
import cPickle as pickle 
import scipy.stats as stats

def beam_height(el, r):
	return r*np.sin(np.radians(el)) + r**2/(2.0*6370)



base_path = '/home/rmet/SPURS2/realtime_radar'
time_thresh = 7200.0


parser = argparse.ArgumentParser(description='Put in a file to be processed')


parser.add_argument('--load', action="store", dest="load", type=bool, default=False)
parser.add_argument('--config', action="store", dest="config", default=None)

pargs = parser.parse_args()

if pargs.config is None:
    cfg = Config.Config('%s/realtime_config.yaml'%(base_path))
    ycfg = cfg.v
else:
    cfg = Config.Config('%s/%s'%(base_path, pargs.config))
    ycfg = cfg.v

if pargs.load:
    print 'trying to load save file'

    load_check = glob.glob('%s/echodata.p'%(base_path))
    if load_check:
		load_data = pickle.load( open('%s/echodata.p'%(base_path), 'rb') )
		times = load_data['times']
		hts = load_data['hts']
		nfiles = load_data['nfiles']

    else:
		print 'not able to load data'
		times = []
		hts = OrderedDict( [('10', []), ('20', []), ('30', []), ('40', []), ('50', [])] )
		nfiles = []
		pass
else:
	print 'making a new file'
	times = []
	hts = OrderedDict( [('10', []), ('20', []), ('30', []), ('40', []), ('50', [])] )
	nfiles = []




cfg = RadarConfig.RadarConfig(dz=ycfg['dz_name'], zdr=ycfg['zdr_name'], kdp=ycfg['kdp_name'], 
                                rho=ycfg['rho_name'], hid='HID')

radar_path = '%s/realtime_cfradial/rhi'%(base_path)

radar_files = np.array(sorted(glob.glob('%s/*.nc'%(radar_path))))


cape = np.load('CAPE.npz')
cape_vals = cape['cape']
cape_time_str = cape['times']
cape_times = np.array([datetime.datetime.strptime(_, '%Y-%m-%d_%H%M') for _ in cape_time_str])

cape_ptop = cape['p_top']
bad_sounding = cape_ptop >= 100.0
cape_vals[bad_sounding] = np.nan


radar_times = np.array( [ datetime.datetime.strptime(string_tools.parse_time_string(os.path.basename(_))[:15], '%Y%m%d_%H%M%S') 
								for _ in radar_files] )
PRF = 1200.0
mur = 3.0e8/(2.0*PRF)


colors = {'10': 'blue', '20': 'green', '30': '#ffcc00', '40': 'orange', '50': 'red'}

for i, ct in enumerate(cape_times):

	if ct not in times: # only do this if this time has not already been logged


		times.append(ct)
		pass
		rt_diffs = np.array([np.abs((ct-_).total_seconds()) for _ in radar_times])

		if rt_diffs.min() >= time_thresh:
			print 'No RHI near this sounding time'
			for __k in hts.keys():
				hts[__k].append(np.nan)
		else:
			print 'sounding time: {}'.format(ct)
			#print 'close RHIs: {}'.format(radar_files[rt_diffs<=time_thresh])

			valid_rhi_files = radar_files[rt_diffs<=time_thresh]
			print '%d RHIs close to sounding time'%(len(valid_rhi_files))

			#times.append(ct)
			nfiles.append(len(valid_rhi_files))



			for k in hts.keys():
				hts[k].append(0)

			#print 'heights: {}'.format(hts)
			#print 'heights keys: {}'.format(hts.keys())



			for vrf in valid_rhi_files:

				#good_radar_file = radar_files[np.argmin(rt_diffs)]
				radar = pyart.io.read(vrf, file_field_names=True)
				gate_spacing = mur/radar.ngates

			# we'll need to go thru each sweep and look for where values are above 30 dBZ


				for s in range(radar.nsweeps):
					dbz_data = radar.get_field(s, 'DZQC')

					for h in hts.keys():
						#print '{} dBZ'.format(h)
						wh_ht = np.where(dbz_data >= float(h))
						if len(wh_ht[0]):
							xyz = radar.get_gate_x_y_z(s)
							beam_range = np.sqrt(xyz[0]**2 + xyz[1]**2)
							zvals = xyz[2]
							z_wh = zvals[wh_ht]

							z_rep_val = np.percentile(z_wh/1000.0, 95)
						#	print '{} dBZ, z value: {}, current max z val: {}'.format(h, z_rep_val, hts[h][-1])

							if z_rep_val > hts[h][-1]:
								hts[h][-1] = z_rep_val
						#		print 'reassigning the max z val: {}'.format(hts[h][-1])
								#time.sleep(0.05)

							else:
						#		print 'not reassigning'
								pass




for _k in hts.keys():
	hts[_k] = np.array(hts[_k])

times = np.array(times)
nfiles = np.array(nfiles)


valid_rhis = np.logical_not(np.isnan(hts['10']))
valid_rhis_and_cape = np.logical_and(valid_rhis, np.logical_not(bad_sounding))


fig, ax = plt.subplots(1,1, figsize=(10,8))
for ih in hts.keys():
	ax.plot(times, hts[ih], color=colors[ih], label='%s dBZ'%(ih), linewidth=2)
	ax.scatter(times, hts[ih], color=colors[ih], s=nfiles**1.5)


ax2 = ax.twinx()
ax2.plot(cape_times, cape_vals, label='CAPE', color='black', linewidth=2)
ax2.set_ylabel('CAPE (J/kg)')



for label in ax.get_xticklabels():
        label.set_rotation(15) 

date_formatter = mdates.DateFormatter('%m/%d %HZ')
# Set the major tick formatter to use your date formatter.
ax.xaxis.set_major_formatter(date_formatter)
# This simply rotates the x-axis tick labels slightly so they fit nicely.
fig.autofmt_xdate()

ax.legend(loc='upper right')
ax2.legend(loc='upper left')
ax.set_xlim(times.min(), times.max())
ax.set_ylabel('Height')
ax.grid(True)
plt.tight_layout()
plt.savefig('echotopheights.png', dpi=120)

###############################################################
###############################################################

valid_capes = cape_vals[valid_rhis_and_cape]



# regressions down here



fig, ax = plt.subplots(1,1, figsize=(8,6))

for ih, h in enumerate(hts.keys()):

	valid_hts = hts[h][valid_rhis_and_cape]

	slope, intercept, r_value, p_value, std_err = stats.linregress(valid_capes, valid_hts)
	dummy_array = np.array([valid_capes.min(), valid_capes.max()])

	ax.scatter(valid_capes, valid_hts, edgecolors='none', c=colors[h], s=50)
	ax.plot(dummy_array, slope*dummy_array+intercept, color=colors[h], label='%s dBZ fit'%(h))


	ax.text(100, 12-ih/2.0, 'slope: %.2e; R$^2$: %.2f'%(slope, r_value**2), color=colors[h])


ax.set_xlim(0, 3000)
ax.set_ylim(-0.5, 14)
ax.grid(True)

ax.set_xlabel('CAPE (J/kg)')
ax.set_ylabel('Height (km)')

plt.tight_layout()

plt.savefig('cape_h30_scatter.png', dpi=120)



##############################################################

echo_data = {'times': times, 'nfiles': nfiles, 'hts': hts}


# save the data
pickle.dump( echo_data, open('%s/echodata.p'%(base_path), 'wb'))