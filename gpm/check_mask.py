import os
import netCDF4 as nc4
import numpy as np

#User inputs
indir = '/home/disk/bob/gpm/sam_ku/classify/ex_data_v05'
years = ['2018']
months = ['06']

num_years = len(years)
num_months = len(months)

for yy in range(0,num_years):
    for mm in range(0,num_months):
        for fname in os.listdir(indir+'/'+years[yy]+'/'+months[mm]):
            print 'fname = ',fname
            ncid = nc4.Dataset(indir+'/'+years[yy]+'/'+months[mm]+'/'+fname,'r')
            bsr_mask = np.array(ncid.variables['bsr_mask_str'])
            dcc_mask = np.array(ncid.variables['dcc_mask_str'])
            dwc_mask = np.array(ncid.variables['dwc_mask_str'])
            wcc_mask = np.array(ncid.variables['wcc_mask_str'])
            storm_mask = np.array(ncid.variables['storm_mask_str'])

            unique_bsr, counts_bsr = np.unique(bsr_mask, return_counts=True)
            unique_dcc, counts_dcc = np.unique(dcc_mask, return_counts=True)
            unique_dwc, counts_dwc = np.unique(dwc_mask, return_counts=True)
            unique_wcc, counts_wcc = np.unique(wcc_mask, return_counts=True)
            unique_sto, counts_sto = np.unique(storm_mask, return_counts=True)
            print 'STR - BSR:',unique_bsr,' DCC:',unique_dcc,' DWC:',unique_dwc,' WCC:',unique_wcc,' STO:',unique_sto

            #bsr_mask = np.array(ncid.variables['bsr_mask_mod'])
            #dcc_mask = np.array(ncid.variables['dcc_mask_mod'])
            #dwc_mask = np.array(ncid.variables['dwc_mask_mod'])
            #wcc_mask = np.array(ncid.variables['wcc_mask_mod'])
            #storm_mask = np.array(ncid.variables['storm_mask_mod'])
            #ncid.close()

            #unique_bsr, counts_bsr = np.unique(bsr_mask, return_counts=True)
            #unique_dcc, counts_dcc = np.unique(dcc_mask, return_counts=True)
            #unique_dwc, counts_dwc = np.unique(dwc_mask, return_counts=True)
            #unique_wcc, counts_wcc = np.unique(wcc_mask, return_counts=True)
            #unique_sto, counts_sto = np.unique(storm_mask, return_counts=True)
            #print 'MOD - BSR:',unique_bsr,' DCC:',unique_dcc,' DWC:',unique_dwc,' WCC:',unique_wcc,' STO:',unique_sto

