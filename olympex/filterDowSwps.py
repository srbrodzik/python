import os
import shutil
import logging as log

log.basicConfig(level=log.INFO,
                format='%(levelname)s %(message)s')

# input directory
swpDir = '/home/disk/bob/olympex/raw/dow/low/20151112/'

# junk sweep directory
junkSwpDir = swpDir + 'junkSweeps/'

# check to make sure junk sweep directory exists; if not, create it
if not os.path.exists(junkSwpDir):
    os.makedirs(junkSwpDir)

# move all sweeps with negative angles to junk sweep dir
for sfile in os.listdir(swpDir):
    if sfile.find('.-0.') >= 0:
        shutil.move(swpDir+sfile,junkSwpDir)
        log.info('file = {} moved to {}'.format(sfile,junkSwpDir))

# move rest of sweeps from volume numbers in junk sweep dir to junk sweep dir
for jfile in os.listdir(junkSwpDir):
    index = jfile.find('_v')
    suffix = jfile[index:]
    log.info('moving sweeps from {} to {}'.format(suffix,junkSwpDir))
    for sfile in os.listdir(swpDir):
        if (sfile.endswith(suffix)):
            shutil.move(swpDir+sfile,junkSwpDir)

# move first sweep of certain rhi volumes to junk sweep dir
for sfile in os.listdir(swpDir):
    if sfile.find('53.3_RHI') >= 0:
        shutil.move(swpDir+sfile,junkSwpDir)
        log.info('file = {} moved to {}'.format(sfile,junkSwpDir))
    elif sfile.find('53.2_RHI') >= 0:
        shutil.move(swpDir+sfile,junkSwpDir)
        log.info('file = {} moved to {}'.format(sfile,junkSwpDir))
