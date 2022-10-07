from sharppy.sharptab import params
from datetime import datetime
import numpy as np
import os

## written by Greg Blumberg - CIMMS
## and
## Kelton Halbert - University of Oklahoma
## wblumberg@ou.edu
## keltonhalbert@ou.edu

def get_mean_pwv(station):
    '''
    Function to get the mean precipitable water vapor (inches) values
    for every month of the year, based on the station passed through.
    
    Parameters
    ----------
    Staion: Can be the 4 letter station ID (i.e. KOUN), 3 letter station 
    ID (i.e. OUN), or the 5 digit WMO ID (i.e. 72357).
    
    Returns
    -------
    mean_pwv: 1D Array of float values corresponding to the mean PWV
    Jan - December.
    
    Written by Greg Blumberg
    and Kelton Halbert.
    
    '''
    ## a rudimentary check to see what type of station identifier
    ## was passed through

    if station == None:
        return np.ma.masked

    if len(station) == 4:
        id_index = 0
        station = station.upper()
    elif len(station) == 3:
        id_index = 2
        station = station.lower()
    elif len(station) == 5:
        id_index = 1
    else:
        #print "Invalid station ID"
        return np.ma.masked

    ## open the file, release the kraken!
    ## get the arrays of station IDs
    pwv_means = np.loadtxt(os.path.dirname( __file__) + '/PW-mean-inches.txt', skiprows=0, dtype=str, delimiter=',')
    sites = pwv_means[1:, id_index]
    ## get the index of the user supplied station - add 1 to account for the header
    
    station_idx = np.where( sites == station )[0]
    if len(station_idx) == 0:
        mean_pwv = np.ma.masked
    else:
        station_idx = station_idx + 1
        ## get the PWV
        mean_pwv = pwv_means[station_idx, 3:][0].astype(np.float)
    return mean_pwv

def get_stdev_pwv(station):
    '''
        Function to get the standard deviation from the mean precipitable water 
        vapor (inches) values for every month of the year, based on the station 
        passed through.
        
        Parameters
        ----------
        Staion: Can be the 4 letter station ID (i.e. KOUN), 3 letter station
        ID (i.e. OUN), or the 5 digit WMO ID (i.e. 72357).
        
        Returns
        -------
        stdev_pwv: 1D Array of float values corresponding to the standard deviation
        from the mean PWV
        Jan - December.
        
        Written by Greg Blumberg
        and Kelton Halbert.
        
        '''
    ## a rudimentary check to see what type of station identifier
    ## was passed through
    if station == None:
        return 0
    if len(station) == 4:
        id_index = 0
        station = station.upper()
    elif len(station) == 3:
        id_index = 2
        station = station.lower()
    elif len(station) == 5:
        id_index = 1
    else:
        #print "Invalid station ID"
        return
    ## open the file, release the kraken!
    ## get the arrays of station IDs
    pwv_stdevs = np.loadtxt(os.path.dirname( __file__) + '/PW-stdev-inches.txt', skiprows=0, dtype=str, delimiter=',')
    sites = pwv_stdevs[1:, id_index]
    ## get the index of the user supplied station - add 1 to account for the header
    station_idx = np.where( sites == station )[0]
    if len(station_idx) == 0:
        stdev_pwv = np.ma.masked
    else:
        station_idx = station_idx + 1
        ## get the PWV
        stdev_pwv = pwv_stdevs[station_idx, 3:][0].astype(np.float)
    return stdev_pwv

def pwv_climo(prof, station, month=None):
    '''
    month is an integer from 1-12
    station_id_3 is the station ID (lower case)
    prof is the profile object

    This function uses the PWV climatology databases provided by Matt Bunkers (NWS/UNR)
    and accepts a SHARPPY profile object, the station name, and the month
    and returns to the user a number indicating where in the distribution the profile's PWV
    value lies.  This function is used in SHARPPY to help provide the user climatological
    context of the PWV index they are viewing.
    
    x can equal 0, 1, 2, or 3
    If the returned value is x, the PWV lies outside +x standard deviations of the mean
    If the returned value is -x, the PWV lies outside -x standard deviations of the mean
    If the returned value is 0, the PWV lies within 1 standard deviation of the mean
    
    Written by Greg Blumberg
    and Kelton Halbert.
    '''
    if not month:
        month = datetime.now().month

    # Calculate the PWV up to 300 mb so it's consistent with the PWV Climo
    pwv_300 = params.precip_water(prof, pbot=None, ptop=300)
    # pwv_300 needs to be in inches (if it isn't already)

    # Load in the PWV mean and standard deviations
    pwv_means = get_mean_pwv(station)
    pwv_stds = get_stdev_pwv(station)
    if pwv_means is np.ma.masked:
        return 0
    elif pwv_means is None:
	return 0
    month_mean = float(pwv_means[month-1])
    month_std = float(pwv_stds[month-1])

    sigma_1 = (month_mean - month_std, month_mean + month_std)
    sigma_2 = (month_mean - (2.*month_std), month_mean + (2.*month_std))
    sigma_3 = (month_mean - (3.*month_std), month_mean + (3.*month_std))

    if len(np.where(pwv_300 > sigma_3[1])[0]) > 0:
        # Means the PWV value is outside +3 sigma of the distribution
        flag = 3
    elif len(np.where(pwv_300 < sigma_3[0])[0]) > 0:
        # Means the PWV value is outside -3 sigma of the distribution
        flag = -3
    elif len(np.where(pwv_300 > sigma_2[1])[0]) > 0:
        # Means the PWV value is outside +2 sigma of the distribution
        flag = 2
    elif len(np.where(pwv_300 < sigma_2[0])[0]) > 0:
        # Means the PWV value is outside -2 sigma of the distribution
        flag = -2
    elif len(np.where(pwv_300 > sigma_1[1])[0]) > 0:
        # Means the PWV value is outside +1 sigma of the distribution
        flag = 1
    elif len(np.where(pwv_300 < sigma_1[0])[0]) > 0:
        # Means the PWV value is outside -1 sigma of the distribution
        flag = -1
    else:
        # Means that the PWV value is within 1 sigma of the distribution
        flag = 0

    return flag

