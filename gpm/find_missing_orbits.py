#!/usr/bin/python3

# Find missing orbits in downloaded and subsetted GPM 2A-Ku files
# File naming convention is this:  2A-UW.GPM.Ku.V8-20180723.20200511-S100918-E114150.035232.V06A.HDF5

import os

indir = '/home/disk/bob/gpm/hdf_subsetted/2Ku/2023/06'

actual_orbits = []
for file in os.listdir(indir):
    if file.endswith('HDF5'):
        parts = file.split('.')
        orbit = int(parts[5])
        actual_orbits.append(orbit)
num_actual_orbits = len(actual_orbits)

first_orbit = sorted(actual_orbits)[0]
last_orbit  = sorted(actual_orbits)[-1]
num_expected_orbits = last_orbit - first_orbit + 1

if num_actual_orbits != num_expected_orbits:
    expected_orbits = range(first_orbit,last_orbit+1)
    missing_orbits = sorted(list(set(expected_orbits) - set(actual_orbits)))
    print('These orbits are missing:')
    print(missing_orbits)
    #f'Missing orbits are: {missing_orbits}'
else:
    print('No orbits are missing')
