# This script will take in diff times and simulate realtime using the driver script

import os
import numpy as np 

date_string = '20170901'
times = ['000030', '000400', '000800', '001200', '001600', '002000', '002400', '002800', '003200', '003600', '004000']
config_file = 'uf_config.yaml'

for t in times:
	os.system('python realtime_driver.py --testtime=%s_%s --config=%s'%(date_string, t, config_file))