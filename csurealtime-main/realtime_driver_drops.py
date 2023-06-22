# Brody Fuchs
# Oct 2017

# This is code to direct the pipeline of code all the way from the original radar data UF file
# that comes in from the server to the final product

# ********* UPDATE ************
# Feb 2018: This is a copy of the older driver script, but needs some new info based on 
# the fact that we're now using DROPS, and have to deal with some of the annoyances that 
# are brought on by that.

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
# Step 6: fixed_process_radar.py --realtime=True is called to read in the gridded and QCed radar data
# 			HID and rainrate calculations are made here.
#
#


def process_pipeline(fname):

    fbase = os.path.basename(fname)
    a = parse_time_string(fbase)

    already_check_ppi = sorted(glob.glob('%s/ppi/%s*%s*%s'%(path2, prefix2, a, ftype2)))
    already_check_rhi = sorted(glob.glob('%s/rhi/%s*%s*%s'%(path2, prefix2, a, ftype2)))

    p_ppi = substring_in_string_list(a, already_check_ppi)
    p_rhi = substring_in_string_list(a, already_check_rhi)

    if already_check_rhi:
    	print 'This RHI has already been processed'
	keep_going = False
    elif already_check_ppi:
    	print 'This PPI has already been processed'
	keep_going = False
    else:
    	print 'This file has not been processed'
	keep_going = True


    if keep_going:


	dummy1 = os.system('python %s/qc_radar.py --file=%s --config=%s'%(base_path, fname, cfg_base))

	# Looks as if os.system waits until the process is done before moving on
	# which is perfect for our purposes
	# now need to check if this is an RHI.



	# okay so after that is done, we need to just do a check to see if the file is there
	#check2 = sorted(glob.glob('%s/%s*.%s'%(path2, prefix2, ftype2)))
	check2 = glob.glob('%s/ppi/%s*%s*%s'%(path2, prefix2, a, ftype2))


    	print 'check2: {}'.format(check2)


	# this will return a string no matter what, so if no files there
	# will return an empty string, so do a check on the len

	if len(check2): # this will fail if it's an RHI
	    print 'Found QCed file in %s'%(path2)

	else:
	    #print 'Did NOT find a file in %s'%(path2)
	    # now check for an RHI with the file base, if so, call the code that will pop up a viewer
	    rhi_qc_check = glob.glob('%s/rhi/%s*%s*%s'%(path2, prefix2, a, ftype2))
	    if len(rhi_qc_check):
	    	# call the viewer script
	    	#print 'Displaying the RHI sweeps'
	    	rdisp_cmd = 'bash %s/pop_up_rhi.sh &'%(base_path)
	    	print 'display command: {}'.format(rdisp_cmd)
	    	os.system(rdisp_cmd)
	    	#os.system('./%s/pop_up_rhi.sh'%(base_path))
	    	pass


	    #sys.exit()


	dummy_2 = os.system('python %s/plot_volume.py --file=%s --config=%s --var=dz_name'%(base_path, check2[0], cfg_base))

	# if we found a file, then run that thru the batch_radx_seapol.py
	# script to grid the data with RadxConvert
	dummy2 = os.system('python %s/batch_radx_seapol.py --file=%s --config=%s'%(base_path, check2[0], cfg_base))


	if pargs.process:
	# now need to check to see if the file is in path3 now
	    check3 = glob.glob('%s/%s*%s*%s'%(path3, prefix3, a, ftype3))

	# now that I have that file, push it thru the fixed grid 
	# radar processing
		#print 'Running the rain/HID processing now'
	    dummy3 = os.system('python -W ignore %s/fixed_process_radar.py --file=%s --realtime=True --config=%s'%(
									base_path, check3[0], cfg_base))
    else:
    	pass



start_dt = datetime.datetime.utcnow()
start = time.time()
print '-------------------------------------------------'
print start_dt

dt_fmt = '%Y%m%d_%H%M%S'

base_path = os.path.dirname(os.path.realpath('__file__'))

parser = argparse.ArgumentParser(description='Put in a file to be processed')

#parser.add_argument('--noarg', action="store_true", default=False)
#parser.add_argument('--file', action="store", dest="file")
parser.add_argument('--realtime', action="store", dest="realtime", type=bool, default=True)
parser.add_argument('--testtime', action="store", dest="testtime")
parser.add_argument('--directory', action="store", dest="directory", default=None)
parser.add_argument('--process', action="store", dest="process", type=bool, default=True)
parser.add_argument('--force', action="store", dest="process", type=bool, default=False)
parser.add_argument('--filetype', action="store", dest="filetype", default=None)
parser.add_argument('--config', action="store", dest="config", default=None)


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

elif pargs.realtime:
	pass
	utctime = datetime.datetime.utcnow()
	dirname = utctime.strftime('%Y%m%d')
	suffix = '/spurs2/%s'%(dirname)


else:
	print 'Need to specify a test time if not operating in realtime'
	sys.exit()



tdiff_thresh = 480. # this is in seconds
pause_time = 15 # time to sleep in seconds




path1 = '%s/%s%s'%(base_path, cfg['path1'], suffix)
ftype1 = cfg['ftype1']
prefix1 = cfg['prefix1']

path2 = '%s/%s'%(base_path, cfg['path2'])
ftype2 = cfg['ftype2']
prefix2 = cfg['prefix2']

path3 = '%s/%s'%(base_path, cfg['path3'])
ftype3 = cfg['ftype3']
prefix3 = cfg['prefix3']

path4 = '%s/%s'%(base_path, cfg['path4'])

raw_files = np.array(sorted(glob.glob('%s/*%s'% (path1, ftype1) )))
file_sizes = np.array([os.path.getsize(_)/1.0e6 for _ in raw_files])


if len(raw_files):
	# if some files within the timeframe, then loop thru them and check if they've been processed.
	# if they have, then do nothing. If they haven't, process them

    valid_file_ind = get_valid_files(raw_files, utctime)
	#print 'valid files: {}'.format(raw_files[valid_file_ind])

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
























