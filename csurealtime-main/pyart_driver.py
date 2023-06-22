# Brody Fuchs
# Oct 2017

# This is code to direct the pipeline of code all the way from the original radar data UF file
# that comes in from the server to the final product

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

# Step 0(?): SCP radar file to radar_data directory
# Step 1: Run qc_radar.py --realtime=True on the raw radar file to do some basic QCing of the data
#			It uses CSU radartools and some combination of standard deviation of phase
# 			Plus despeckling and an insect filter. It also filters rho and ZDR as well
# 			Finally, it writes out a cfradial file into realtime_cfradial
# Step 2: batch_radx_seapol.py --realtime=True is called to use RadxConver to grid the QCed
# 			cfradial file, it is then output into realtime_gridded
# 			The code also renames the file with the same timestamp as the cfradial file
# Step 3: fixed_process_radar.py --realtime=True is called to read in the gridded and QCed radar data
# 			HID and rainrate calculations are made here.


def process_pipeline(fname):


    fbase = os.path.basename(fname)
    a = gtools.parse_time_string(fbase)

    date_string = a[:8]

    print 'a: {}, date string: {}'.format(a, date_string)

    if True:

	qc_cmd = 'python -W ignore %s/qc_radar.py --file=%s --config=%s'%(base_path, fname, cfg_base)
	print 'RUNNING: %s'%qc_cmd
	raw_check = os.system(qc_cmd)


	#dummy1 = os.system('python %s/qc_radar.py --file=%s --config=%s'%(base_path, fname, cfg_base))

	# Looks as if os.system waits until the process is done before moving on
	# which is perfect for our purposes
	# now need to check if this is an RHI.


#	check_qc = glob.glob('%s/ppi/%s*%s*%s'%(qc_data_path, qc_prefix, a, qc_ftype))
	check_qc = glob.glob('%s/*/%s*%s*%s'%(qc_data_path, qc_prefix, a, qc_ftype))

    	print 'qc check: {}'.format(check_qc)


	if len(check_qc):
	    plot_cmd = 'python -W ignore %s/plot_polar_radar.py --file=%s --config=%s'%(base_path, check_qc[0], cfg_base)	
	    print 'RUNNING: %s'%plot_cmd
	    os.system(plot_cmd)


	    # now need to do a check for PPI file or not

	blah = '%s/%s*%s*ppi*%s'%(qc_data_path, qc_prefix, a, qc_ftype)

	check_polar_ppi = glob.glob('%s/ppi/%s*%s*ppi*%s'%(qc_data_path, qc_prefix, a, qc_ftype))

	print 'PPI polar check: {}'.format(check_polar_ppi)


	if len(check_polar_ppi):

	    pv_cmd = 'python -W ignore %s/plot_volume.py --file=%s --config=%s --var=dz_name'%(base_path, check_polar_ppi[0], cfg_base)
	    print 'RUNNING: %s'%pv_cmd
	    os.system(pv_cmd)
	    #dummy_2 = os.system('python %s/plot_volume.py --file=%s --config=%s --var=dz_name'%(base_path, check2[0], cfg_base))

	# if we found a file, then run that thru the batch_radx_seapol.py
	# script to grid the data with RadxConvert
	    #dummy2 = os.system('python %s/batch_radx_seapol.py --file=%s --config=%s'%(base_path, check2[0], cfg_base))

	    grid_cmd = 'python -W ignore %s/batch_radx_seapol.py --file=%s --config=%s'%(base_path, check_polar_ppi[0], cfg_base)
	    print 'RUNNING: %s'%grid_cmd
	    os.system(grid_cmd)


	    if pargs.process:
	    # now need to check to see if the file is in path3 now
	    # Check for a gridded file
	    #check_qc = glob.glob('%s/ppi/%s*%s*%s'%(qc_data_path, qc_prefix, a, qc_ftype))
		check_grid = glob.glob('%s/%s*%s*%s'%(grid_data_path, grid_prefix, a, grid_ftype))
		#print '%s/%s*%s*%s'%(grid_data_path, grid_prefix, a, grid_ftype)
		#print 'check grid length: {}'.format(len(check_grid))

		if len(check_grid):
	    
		    process_cmd = 'python -W ignore %s/process_gridded_radar.py --file=%s --realtime=1 --config=%s'%(\
						base_path, check_grid[0], cfg_base)
		    print 'RUNNING: %s'%process_cmd
		    os.system(process_cmd)

	    print 'pargs.rain: {}'.format(pargs.rain)

	    if pargs.rain:
		check_process = glob.glob('%s/%s*%s*%s'%(rain_data_path, rain_prefix, a, rain_ftype))
		print '%s/%s*%s*%s'%(rain_data_path, rain_prefix, a, rain_ftype)
		print 'len check process: {}'.format(len(check_process))


		if len(check_process):
		    pass
		    rain_cmd = 'python -W ignore %s/rainfall_products.py --file %s --config %s'%(base_path, check_process[0], cfg_base)
		    print 'RUNNING: %s'%rain_cmd
		    os.system(rain_cmd)



	else: # It must be an RHI
	    if pargs.show_rhi:
	    # now check for an RHI with the file base, if so, call the code that will pop up a viewer
		rhi_qc_check = glob.glob('%s/rhi/%s*%s*rhi*%s'%(qc_data_path, qc_prefix, a, qc_ftype))
	    	if len(rhi_qc_check):
	    # 	# call the viewer script
	     	    rdisp_cmd = 'bash %s/pop_up_rhi.sh &'%(base_path)
	     	    print 'display command: {}'.format(rdisp_cmd)
	     	    os.system(rdisp_cmd)



    else:
    	pass



start_dt = datetime.datetime.utcnow()
start = time.time()
print '-------------------------------------------------'
print start_dt

dt_fmt = '%Y%m%d_%H%M%S'


base_path = os.path.dirname(os.path.realpath(__file__))
print 'driver base path: {}'.format(base_path)

parser = argparse.ArgumentParser(description='Put in a file to be processed')

parser.add_argument('--file', action="store", dest="file", default=None)
parser.add_argument('--realtime', action="store", dest="realtime", type=int, default=0)
parser.add_argument('--testtime', action="store", dest="testtime")
parser.add_argument('--directory', action="store", dest="directory", default=None)
parser.add_argument('--process', action="store", dest="process", type=int, default=1)
parser.add_argument('--force', action="store", dest="force", type=int, default=0)
parser.add_argument('--filetype', action="store", dest="filetype", default=None)
parser.add_argument('--config', action="store", dest="config", default=None)
parser.add_argument('--rain', action="store", dest="rain", type=int, default=1)
parser.add_argument('--start', action="store", dest="start", default=None)
parser.add_argument('--end', action="store", dest="end", default=None)
parser.add_argument('--show_rhi', action="store", dest="show_rhi", type=int, default=0)
parser.add_argument('--file_path', action="store", dest="file_path", default=None)
parser.add_argument('--overwrite', action="store", dest="overwrite", type=int, default=0)


pargs = parser.parse_args()


if pargs.config is None:
    cfg_file = '%s/realtime_config.yaml'%(base_path)
    cfg = Config.Config(cfg_file)
    cfg = cfg.v
else:
    cfg_file = '%s/%s'%(base_path, pargs.config)
    cfg = Config.Config(cfg_file)
    cfg = cfg.v



#tdiff_thresh = 480. # this is in seconds
#pause_time = 15 # time to sleep in seconds

tdiff_thresh = cfg['file_tdiff_thresh']
pause_time = cfg['pause_time']

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
#    dirname = utctime.strftime('%Y%m%d')
#    suffix = '/spurs2/%s'%(dirname)
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






raw_data_path = cfg['raw_data_path']
raw_ftype = cfg['raw_ftype']

qc_data_path = cfg['qc_data_path']
qc_ftype = cfg['qc_ftype']
qc_prefix = cfg['qc_prefix']

grid_data_path = cfg['grid_data_path']
grid_ftype = cfg['grid_ftype']
grid_prefix = cfg['grid_prefix']

rain_data_path = cfg['rain_data_path']
rain_ftype = cfg['rain_ftype']
rain_prefix = cfg['rain_prefix']

raw_files = np.array(sorted(glob.glob('%s/????????/*%s'% (raw_data_path, raw_ftype) )))
file_sizes = np.array([os.path.getsize(_)/1.0e6 for _ in raw_files])


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
	    print 'Testing file: {}'.format(vf)

	    if not pargs.overwrite:

#		print 'CHECKING IF FILE EXISTS'
	    	vbase = os.path.basename(vf)[3:]
#		print '%s/*/*%s*'%(qc_data_path, vbase)
	    	already_exist = glob.glob('%s/*/*%s*'%(qc_data_path, vbase))
		print 'already exist? {}'.format(len(already_exist))

		time.sleep(3)
		if len(already_exist) == 0:
		    print 'PROCESSING FILE'
		    process_pipeline(vf)
		else:
		    print '%s has already been processed so will not overwrite'%(vbase)

	    else:
	    	process_pipeline(vf)
			


    else:
	print 'Exiting'
	sys.exit()


	# If you've made it this far, you have a file. Now run the qc_radar on it
	#print 'You made it here!'
	#print 'Running QC on {}'.format(vf)
	#process_pipeline(fname)


print 'All done in %.1f seconds'%(time.time() - start)
























