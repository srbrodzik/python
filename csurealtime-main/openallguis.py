# this code opens the GUIs used for the realtime processing

import os

base_path = '/home/rmet/SPURS2/realtime_radar'

os.system('python %s/tklooper_update.py --path=%s/figures/raw_ppi --n_images=12 &' %(base_path, base_path))

os.system('python %s/tklooper_update.py --path=%s/figures/rain --n_images=12 &' %(base_path, base_path))

os.system('python %s/tklooper_update.py --path=%s/figures/large_dbz_vel_ppi --n_images=12 &' %(base_path, base_path))

os.system('python %s/processing_status.py &'%(base_path))