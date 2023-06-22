# Brody Fuchs
# Oct 2017

# This is code to direct the pipeline of code all the way from the original radar data UF file
# that comes in from the server to the final product

# ********* UPDATE ************
# Feb 2018: This is a copy of the older driver script, but needs some new info based on 
# the fact that we're now using DROPS, and have to deal with some of the annoyances that 
# are brought on by that.

# This is a driver of QC'ing data that will be used after the fact. And will take in times
# to process (I think)

import numpy as np
import os
import glob
import subprocess as sp
import time
import datetime
import argparse
import sys
import Config
import gentools as gtools 
import plot_tools as ptools 
import radar_tools as radtools 
from csuram import analysis_tools as atools

# Feb '18: this is updated to reflect the new flow of the code
# Step 0(?): rsync radar file to radar_data directory
# Step 1: Run qc_radar.py --file=$filename --config=$configfile on the raw radar file to do some basic QCing of the data
#			It uses CSU radartools to read in the data and do some basic checks on the data
#			Like turning the 'UNKNOWN_82' variable into the correlation coefficent
#			And fixing the fact that the RHIs have a huge bug with missing the range factor in the radar equation.
#			It also writes out the file in cfradial format
#
# Step 2: uf_nc_convert_drops.py --file=$filename --config=$configfile on the cfradial file output from the qc_radar.py script
# 			This converts the cfradial file to a UF format, then runs DROPS on the UF file, then converts that output file
#			back to a cfradial file so that it can be used by scripts downstream of it.
#			It also does some moving of files around and changing names because Radxconvert and DROPS do some weird stuff
# 			with directories and filenames
# 
# Step 3: plot_drops_radar.py --file=$filename --config=$configfile has the code that used to be in qc_radar.py that took care of the plotting
#			So this gets run after DROPS and after the file has been converted back to a cfradial
#
# Step 4: plot_volume.py --file=$filename --config=$configfile --var=$varname will take in one variable and plot all sweeps of the data.
#			Originally used to just plot all PPI tilts at once 
#
#
# Step 5: batch_radx_seapol.py --file=$filename --realtime=True is called to use RadxConver to grid the QCed
# 			cfradial file, it is then output into realtime_gridded
# 			The code also renames the file with the same timestamp as the cfradial file
#
# Step 6: process_gridded_radar.py --realtime=True is called to read in the gridded and QCed radar data
# 			HID and rainrate calculations are made here.
#
#


# parse_time_string(filename)


def process_pipeline(fname):

    fbase = os.path.basename(fname)
    a = gtools.parse_time_string(fbase)
    print 'a: {}'.format(a)

    date_string = a[:8]
    print 'date string: {}'.format(date_string)


#    already_check_ppi = sorted(glob.glob('%s/ppi/%s*%s*%s'%(path2, prefix2, a, ftype2)))
#    already_check_rhi = sorted(glob.glob('%s/rhi/%s*%s*%s'%(path2, prefix2, a, ftype2)))
#
#    p_ppi = substring_in_string_list(a, already_check_ppi)
#    p_rhi = substring_in_string_list(a, already_check_rhi)

#    if already_check_rhi:
#    	print 'This RHI has already been processed'
#	keep_going = False
#    elif already_check_ppi:
#    	print 'This PPI has already been processed'
#	keep_going = False
#    else:
#    	print 'This file has not been processed'
#	keep_going = True


    if True:

	qc_cmd = 'python %s/qc_radar.py --file=%s --config=%s'%(base_path, fname, cfg_base)
	print 'RUNNING: %s'%qc_cmd
	raw_check = os.system(qc_cmd)

	# Looks as if os.system waits until the process is done before moving on
	# which is perfect for our purposes
	# now need to check if this is an RHI.


	# okay so after that is done, we need to just do a check to see if the file is there
	#check2 = sorted(glob.glob('%s/%s*.%s'%(path2, prefix2, ftype2)))



#	check_qc = glob.glob('%s/ppi/%s*%s*%s'%(qc_data_path, qc_prefix, a, qc_ftype))
	check_qc = glob.glob('%s/*/%s*%s*%s'%(qc_data_path, qc_prefix, a, qc_ftype))

    	print 'qc check: {}'.format(check_qc)

	# ****************************************************************************
	# ** NEXT IS THE CONVERSION TO UF AND DROPS AND THEN BACK TO A NETCDF FILE
	# ****************************************************************************

	if len(check_qc):
	# this has to be done regardless of PPI or RHI
	    uf_nc_cmd = 'python %s/uf_nc_convert_drops.py --file=%s --config=%s'%(base_path, check_qc[0], cfg_base)
	    print 'RUNNING: %s'%uf_nc_cmd
	    os.system(uf_nc_cmd)

	else:
	    pass



	check_drops_general = glob.glob('%s/%s*%s*%s'%(drops_data_path, drops_prefix, a, drops_ftype))

	print 'DROPS check: {}'.format(check_drops_general)
	if len(check_drops_general):
	    drops_plot_cmd = 'python %s/plot_drops_radar.py --file=%s --config=%s'%(base_path, check_drops_general[0], cfg_base)	
	    print 'RUNNING: %s'%drops_plot_cmd
	    os.system(drops_plot_cmd)




	check_drops_ppi = glob.glob('%s/%s*%s*ppi*%s'%(drops_data_path, drops_prefix, a, drops_ftype))

	print 'PPI drops check: {}'.format(check_drops_ppi)


	# this will return a string no matter what, so if no files there
	# will return an empty string, so do a check on the len

	if len(check_drops_ppi): # this will fail if it's an RHI
	    print 'Found QCed file in %s'%(drops_data_path)

	# now need to plot???
	# All of the subsequent code is based on the presence of the drops file

	# *******************
	    pv_cmd = 'python %s/plot_volume.py --file=%s --config=%s --var=dz_name'%(base_path, check_drops_ppi[0], cfg_base)
	    print 'RUNNING: %s'%pv_cmd
	    os.system(pv_cmd)

	# now need to do the gridding

	    grid_cmd = 'python %s/batch_radx_seapol.py --file=%s --config=%s'%(base_path, check_drops_ppi[0], cfg_base)
	    print 'RUNNING: %s'%grid_cmd
	    os.system(grid_cmd)


	    if pargs.process:
	    # now need to check to see if the file is in path3 now
	    # Check for a gridded file
	    #check_qc = glob.glob('%s/ppi/%s*%s*%s'%(qc_data_path, qc_prefix, a, qc_ftype))
		check_grid = glob.glob('%s/%s*%s*%s'%(grid_data_path, grid_prefix, a, grid_ftype))

		if len(check_grid):
	    
		    process_cmd = 'python -W ignore %s/process_gridded_radar.py --file=%s --realtime=True --config=%s'%(base_path, check_grid[0], cfg_base)
		    print 'RUNNING: %s'%process_cmd
		    os.system(process_cmd)


	    if pargs.rain:
		check_process = glob.glob('%s/%s*%s*%s'%(rain_data_path, rain_prefix, a, rain_ftype))

		if len(check_process):
		    pass
		    rain_cmd = 'python %s/rainfall_products.py --file %s --config %s'%(base_path, check_process[0], cfg_base)
		    print 'RUNNING: %s'%rain_cmd
		    os.system(rain_cmd)


	else:
	    #print 'Did NOT find a file in %s'%(path2)
	    # now check for an RHI with the file base, if so, call the code that will pop up a viewer
	#    rhi_qc_check = glob.glob('%s/rhi/%s*%s*%s'%(drops_data_path, drops_prefix, a, drops_ftype))
	    check_drops_rhi = glob.glob('%s/%s*%s*rhi*%s'%(drops_data_path, drops_prefix, a, drops_ftype))
	    if len(check_drops_rhi):
	    	# call the viewer script
	    	#print 'Displaying the RHI sweeps'
#	    	rdisp_cmd = 'bash %s/pop_up_rhi.sh &'%(base_path)
#	    	print 'display command: {}'.format(rdisp_cmd)
#	    	os.system(rdisp_cmd)
#	    	#os.system('./%s/pop_up_rhi.sh'%(base_path))
	    	pass

	# Regardless, run the plot_radar_stats code?
	if pargs.stats:
	    stats_cmd = 'python %s/plot_radar_stats.py --date %s'%(base_path, date_string)
	    print 'running %s'%(stats_cmd)
	    os.system(stats_cmd)





    else:
    	pass



start_dt = datetime.datetime.utcnow()
start = time.time()
print '-------------------------------------------------'
print start_dt

dt_fmt = '%Y%m%d_%H%M%S'
tdiff_thresh = 480. # this is in seconds
pause_time = 15 # time to sleep in seconds

base_path = os.path.dirname(os.path.realpath('__file__'))

parser = argparse.ArgumentParser(description='Put in a file to be processed')

parser.add_argument('--realtime', action="store", dest="realtime", type=int, default=0)
parser.add_argument('--testtime', action="store", dest="testtime")
parser.add_argument('--directory', action="store", dest="directory", default=None)
parser.add_argument('--process', action="store", dest="process", type=int, default=1)
parser.add_argument('--force', action="store", dest="force", type=int, default=0)
parser.add_argument('--rain', action="store", dest="rain", type=int, default=1)
parser.add_argument('--filetype', action="store", dest="filetype", default=None)
parser.add_argument('--config', action="store", dest="config", default=None)
parser.add_argument('--start', action="store", dest="start", default=None)
parser.add_argument('--end', action="store", dest="end", default=None)
parser.add_argument('--stats', action="store", dest="stats", type=int, default=0)




pargs = parser.parse_args()


if pargs.config is None:
    cfg_file = '%s/realtime_config.yaml'%(base_path)
    cfg = Config.Config(cfg_file)
    cfg = cfg.v
else:
    cfg_file = '%s/%s'%(base_path, pargs.config)
    cfg = Config.Config(cfg_file)
    cfg = cfg.v


cfg_base = os.path.basename(cfg_file)

#print 'config file: {}'.format(cfg_file)


# realtime will be the default mode here, but also have a test mode with the test data
if pargs.testtime is not None:
	utctime = datetime.datetime.strptime(pargs.testtime, dt_fmt)
	if pargs.directory is not None:
		suffix = '/%s'%(pargs.directory)
	else:
		suffix = ''
	end_time = utctime
	start_time = end_time - datetime.timedelta(minutes=1)

elif pargs.realtime:
	pass
	utctime = datetime.datetime.utcnow()
	dirname = utctime.strftime('%Y%m%d')
	suffix = '/spurs2/%s'%(dirname)
	end_time = utctime
	start_time = end_time - datetime.timedelta(seconds=tdiff_thresh)

elif pargs.start is not None and pargs.end is not None:
    try:
    	start_time = datetime.datetime.strptime(pargs.start, dt_fmt)
    	end_time = datetime.datetime.strptime(pargs.end, dt_fmt)
    except Exception:
    	print 'You need to specify start and end times with the format: %s'%(dt_fmt)
    	sys.exit()

else:
	print 'Need to specify a test time if not operating in realtime'
	sys.exit()


# # realtime will be the default mode here, but also have a test mode with the test data
# if pargs.testtime is not None:
#     utctime = datetime.datetime.strptime(pargs.testtime, dt_fmt)
#     if pargs.directory is not None:
# 	suffix = '/%s'%(pargs.directory)
#     else:
# 	suffix = ''

# elif pargs.realtime:
#     pass
#     utctime = datetime.datetime.utcnow()
#     dirname = utctime.strftime('%Y%m%d')
#     suffix = '/spurs2/%s'%(dirname)


# else:
#     print 'Need to specify a test time if not operating in realtime'
#     sys.exit()


# try:
#     start_time = datetime.datetime.strptime(pargs.start, dt_fmt)
#     end_time = datetime.datetime.strptime(pargs.end, dt_fmt)
# except Exception:
#     print 'You need to specify start and end times with the format: %s'%(dt_fmt)
#     sys.exit()



#tdiff_thresh = 480. # this is in seconds
#pause_time = 15 # time to sleep in seconds



raw_data_path = cfg['raw_data_path']
raw_ftype = cfg['raw_ftype']

qc_data_path = cfg['qc_data_path']
qc_ftype = cfg['qc_ftype']
qc_prefix = cfg['qc_prefix']

drops_data_path = cfg['drops_data_path']
drops_ftype = cfg['drops_ftype']
drops_prefix = cfg['drops_prefix']

grid_data_path = cfg['grid_data_path']
grid_ftype = cfg['grid_ftype']
grid_prefix = cfg['grid_prefix']

rain_data_path = cfg['rain_data_path']
rain_ftype = cfg['rain_ftype']
rain_prefix = cfg['rain_prefix']


#path1 = '%s/%s%s'%(base_path, cfg['path1'], suffix)
#prefix1 = cfg['prefix1']

#path2 = '%s/%s'%(base_path, cfg['path2'])
#ftype2 = cfg['ftype2']
#prefix2 = cfg['prefix2']

#path3 = '%s/%s'%(base_path, cfg['path3'])
#ftype3 = cfg['ftype3']
#prefix3 = cfg['prefix3']

#path4 = '%s/%s'%(base_path, cfg['path4'])

raw_files = np.array(sorted(glob.glob('%s/????????/*%s'% (raw_data_path, raw_ftype) )))
#file_sizes = np.array([os.path.getsize(_)/1.0e6 for _ in raw_files])

print 'number of possible files: {}'.format(len(raw_files))


if len(raw_files):
	# if some files within the timeframe, then loop thru them and check if they've been processed.
	# if they have, then do nothing. If they haven't, process them

    file_time_strings = np.array( [gtools.parse_time_string(os.path.basename(_)) for _ in raw_files] ) 
    file_times = np.array( [datetime.datetime.strptime(_, dt_fmt) for _ in file_time_strings] )

    # now need to figure out which files are between the beginning and the end

    valid_file_ind = np.where( (file_times >= start_time) & (file_times <= end_time) )

    if len(valid_file_ind) > 0: # It should equal 1
		# loop thru here

    	valid_files = raw_files[valid_file_ind]
	print 'temporally valid files: {}'.format(valid_files)

	for vf in valid_files:
	    #vf = raw_files[valid_file_ind[-1]]
	    print 'Testing file: {}'.format(vf)
			#vbase = os.path.basename(vf)[3:]
			#print 'vbase: {}'.format(vbase)
	    process_pipeline(vf)
			


    else:
	print 'Exiting'
	sys.exit()


	# If you've made it this far, you have a file. Now run the qc_radar on it
	#print 'You made it here!'
	#print 'Running QC on {}'.format(vf)
	#process_pipeline(fname)


print 'All done in %.1f seconds'%(time.time() - start)
























