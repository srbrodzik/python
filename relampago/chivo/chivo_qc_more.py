# -*- coding: utf-8 -*-
"""
This code applies the DZ_qc mask to HydroClass and RainRate

Created June 28, 2021

@author: mrocq
"""

#import the goods
import numpy as np
import pyart
import os

#Adds a newly created field to the radar object.
def add_field_to_radar_object(field, radar, field_name='FH', units='unitless', 
                              long_name='Hydrometeor ID', standard_name='Hydrometeor ID',
                              dz_field='ZC'):

    # should actually pass in fill_value
    fill_value = -32768

    # THIS WHOLE SECTION SEEMS REDUNDANT - WE ALREADY KNOW WE PASSED IN A MASKED ARRAY
    # BASED ON DZ_qc MASK
    # ensures that field is a masked array
    masked_field = np.ma.asanyarray(field)
    # fill new mask looking for fill_values
    masked_field.mask = masked_field == fill_value
    if hasattr(radar.fields[dz_field]['data'], 'mask'):
        setattr(masked_field, 'mask', 
                np.logical_or(masked_field.mask, radar.fields[dz_field]['data'].mask))
        
    field_dict = {'data': masked_field,
                  'units': units,
                  'long_name': long_name,
                  'standard_name': standard_name,
                  '_FillValue': fill_value}
    radar.add_field(field_name, field_dict, replace_existing=True)
    return radar

#Apply masked field from filtered reflectivity to other parameters
def radar_qc_more(file):
       
    radar = pyart.io.read_cfradial(file, file_field_names=True)
    dbz = radar.fields['DZ_qc']['data']
    dbz_masked = np.ma.getmask(dbz)
    hydro = radar.fields['HydroClass']['data']
    hcq = radar.fields['HCQuality']['data']
    
    bad = -32768
    # can use this instead
    # bad = dbz.fill_value
    hydro_qc = 1.0 * hydro
    hydro_qc[dbz_masked] = bad
    hcq_qc = 1.0 * hcq
    hcq_qc[dbz_masked] = bad

       
    #Add new QCed hydroclass field to radar file
    radar_qced = add_field_to_radar_object(hydro_qc, radar, field_name='HydroClass_qc', units='', long_name='Hydrometeor Class (QCed)',
                                   standard_name='HydroClass (QCed)', dz_field='corrected_reflectivity')
    
    #Add new QCed hcquality field to radar file
    radar_qced = add_field_to_radar_object(hcq_qc, radar, field_name='HCQuality_qc', units='%', long_name='HCQuality (QCed)',
                                   standard_name='HCQuality (QCed)', dz_field='corrected_reflectivity')
   
    
    return (radar_qced)

    
    
    
    
    
