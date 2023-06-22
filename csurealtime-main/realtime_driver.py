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


def get_valid_files(filenames, t):
	file_bases = np.array([os.path.basename(_) for _ in filenames])
	#print file_bases[:5]
	file_times = np.array([datetime.datetime.strptime(_[len(prefix1):-1*(len(ftype1)+1)], dt_fmt) for _ in file_bases])

	file_time_diffs = np.array([(t - _).total_seconds() for _ in file_times])

	# these are the differences between the current (or specified time) and the file times in seconds
	# We want to see if there is a file that is greater than 0, but is also less than some threshold value

	vfi = np.where( (file_time_diffs >= 0) & (file_time_diffs < tdiff_thresh))[0]
	return vfi



def is_process_running():
	pass

		# p1 = sp.Popen(['ls', '-th', '%s/%s'%(outpath, stripped_date)], stdout=sp.PIPE)
		# p2 = sp.Popen(['head', '-1'], stdin=p1.stdout, stdout=sp.PIPE)
		# p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
		# out, err = p2.communicate()


# def grab_recent_file(filepath):
# 	p1 = sp.Popen(['ls', '-th', '%s'%(filepath)], stdout=sp.PIPE)
# 	p2 = sp.Popen(['head', '-1'], stdin=p1.stdout, stdout=sp.PIPE)
# 	p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
# 	out, err = p2.communicate()

# 	return out.strip()


def grab_files_time_sorted(path_string):
	searchedfile = glob.glob(path_string)
	files = sorted(searchedfile, key=lambda file: os.path.getctime(file))
	return files

def grab_recent_file(path_string):
	files = grab_files_time_sorted(path_string)
	return files[-1]

start_dt = datetime.datetime.utcnow()
start = time.time()
print '-------------------------------------------------'
print start_dt

dt_fmt = '%Y%m%d_%H%M%S'

base_path = '/home/rmet/SPURS2/realtime_radar'

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

print 'config file: {}'.format(cfg_file)


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



tdiff_thresh = 240. # this is in seconds
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

raw_files = sorted(glob.glob('%s/*%s'% (path1, ftype1) ))

if len(raw_files):
	# if some files are in the directory, then check the time of the most recent file

	valid_file_ind = get_valid_files(raw_files, utctime)
	#print 'valid files: {}'.format(raw_files[valid_file_ind])

	if len(valid_file_ind) > 0: # It should equal 1
		vf = raw_files[valid_file_ind[-1]]
		print 'There is a valid file: {}'.format(vf)
		vbase = os.path.basename(vf)[3:]
		print 'vbase: {}'.format(vbase)

	elif (len(valid_file_ind) == 0) and (pargs.realtime):
		# if it's running in realtime mode but no valid files are present
		print 'No valid files yet'
		done = False
		for i in range(5):
			print 'Attempt %d, waiting %d seconds...'%(i+1, pause_time)
			time.sleep(pause_time)

			raw_files = sorted(glob.glob('%s/*%s'% (path1, ftype1) ))
			valid_file_ind = get_valid_files(raw_files, utctime)
			if len(valid_file_ind):
				done = True
				break

		if done:
			vf = raw_files[valid_file_ind[0]]
			print 'Got the file: {}'.format(vf)
			pass
		else:
			print 'Waited but could not find the file, exiting'
			sys.exit()

	else:
		print 'Exiting'
		sys.exit()


	# If you've made it this far, you have a file. Now run the qc_radar on it
	#print 'You made it here!'
	#print 'Running QC on {}'.format(vf)
	already_check = glob.glob('%s/ppi/%s*%s'%(path2, prefix2, ftype2))


	dummy1 = os.system('python %s/qc_radar.py --file=%s --config=%s'%(base_path, vf, cfg_base))

	# Looks as if os.system waits until the process is done before moving on
	# which is perfect for our purposes


	# okay so after that is done, we need to just do a check to see if the file is there
	#check2 = sorted(glob.glob('%s/%s*.%s'%(path2, prefix2, ftype2)))
	check2 = grab_recent_file('%s/ppi/%s*%s'%(path2, prefix2, ftype2))


    	print 'check2: {}'.format(check2)


	# this will return a string no matter what, so if no files there
	# will return an empty string, so do a check on the len
	if len(check2):
		print 'Found QCed file in %s'%(path2)

	else:
		print 'Did NOT find a file in %s'%(path2)
		sys.exit()

	dummy_2 = os.system('python %s/plot_volume.py --file=%s --config=%s --var=dz_name'%(base_path, check2, cfg_base))

	# if we found a file, then run that thru the batch_radx_seapol.py
	# script to grid the data with RadxConvert
	dummy2 = os.system('python %s/batch_radx_seapol.py --file=%s --config=%s'%(base_path, check2, cfg_base))


	if pargs.process:
	# now need to check to see if the file is in path3 now
		check3 = grab_recent_file('%s/%s*%s'%(path3, prefix3, ftype3))

	# now that I have that file, push it thru the fixed grid 
	# radar processing
		#print 'Running the rain/HID processing now'
		dummy3 = os.system('python -W ignore %s/fixed_process_radar.py --file=%s --realtime=True --config=%s'%(
									base_path, check3, cfg_base))


print 'All done in %.1f seconds'%(time.time() - start)
























