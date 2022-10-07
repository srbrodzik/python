#!/usr/bin/python3

import os
import netCDF4 as nc4
import numpy as np

indir = '/home/disk/bob/gpm/nam_ku/classify/ex_data_v07'
months = ['01','02','03','04','05','06',
          '07','08','09','10','11','12']
years = ['2015','2016','2017','2018','2019']

minDiameter = 99
maxDiameter = -99

fout = open(indir+'/dm_info','w')

for iyear in years:
    for imonth in months:
        print('year = {} and month = {}: '.format(iyear,imonth))
        minMonth = 99
        maxMonth = -99
        for fname in os.listdir(indir+'/'+iyear+'/'+imonth):
            if fname.endswith('nc'):

                #print('{}'.format(fname))
                #fout.write('fname = {}\t'.format(fname))
 
                # open and read dsd_dm from input file
                ncid = nc4.Dataset(indir+'/'+iyear+'/'+imonth+'/'+fname,'r')
                dm_id = ncid.variables['dsd_dm']
                dm = dm_id[:]
                ncid.close()

                # find min and max dsd_dm
                minDm = np.min(dm)
                maxDm = np.max(dm)
                #fout.write('min = {:.2f}\tmax = {:.2f}\n'.format(minDm,maxDm))
                
                if minDm < minMonth:
                    minMonth = minDm
                if maxDm > maxMonth:
                    maxMonth = maxDm

                if minDm < minDiameter:
                    minDiameter = minDm
                if maxDm > maxDiameter:
                    maxDiameter = maxDm

        print('min = {:.2f} and max = {:.2f}\n'.format(minMonth,maxMonth))
        fout.write('{}{}: min = {:.2f} and month = {:.2f}\n'.format(iyear,imonth,minMonth,maxMonth))

print('min/max Dm = {:.2f}/{:.2f}'.format(minDiameter,maxDiameter))
fout.write('\nALL YEARS AND MONTHS: min = {:.2f} and max = {:.2f}\n'.format(minDiameter,maxDiameter))

fout.close()
