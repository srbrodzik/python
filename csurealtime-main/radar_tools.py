import pyart
import numpy as np 
import glob
import os



def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    TJL - taken from old Py-ART version
    Calculate Cartesian coordinate from radar coordinates
    Parameters
    ----------
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    az : array
        Azimuth angle of the radar in degrees.
    ele : array
        Elevation angle of the radar in degrees.
    Returns
    -------
    x, y, z : array
        Cartesian coordinates in meters from the radar.
    Notes
    -----
    The calculation for Cartesian coordinate is adapted from equations
    2.28(b) and 2.28(c) of Doviak and Zrnic [1]_ assuming a
    standard atmosphere (4/3 Earth's radius model).
    .. math::
        z = \\sqrt{r^2+R^2+r*R*sin(\\theta_e)} - R
        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})
        x = s * sin(\\theta_a)
        y = s * cos(\\theta_a)
    Where r is the distance from the radar to the center of the gate,
    :math:\\theta_a is the azimuth angle, :math:\\theta_e is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).
    References
    ----------
    .. [1] Doviak and Zrnic, Doppler Radar and Weather Observations, Second
        Edition, 1993, p. 21.
    """
    theta_e = ele * np.pi / 180.0  # elevation angle in radians.
    theta_a = az * np.pi / 180.0  # azimuth angle in radians.
    R = 6371.0 * 1000.0 * 4.0 / 3.0  # effective radius of earth in meters.
    r = rng * 1000.0  # distances to gates in meters.

    z = (r ** 2 + R ** 2 + 2.0 * r * R * np.sin(theta_e)) ** 0.5 - R
    s = R * np.arcsin(r * np.cos(theta_e) / (R + z))  # arc length in m.
    x = s * np.sin(theta_a)
    y = s * np.cos(theta_a)
    return x, y, z


def get_z_from_radar(radar):
    """Input radar object, return z from radar (km, 2D)"""
    azimuth_1D = radar.azimuth['data']
    elevation_1D = radar.elevation['data']
    srange_1D = radar.range['data']
    sr_2d, az_2d = np.meshgrid(srange_1D, azimuth_1D)
    el_2d = np.meshgrid(srange_1D, elevation_1D)[1]
    xx, yy, zz = radar_coords_to_cart(sr_2d/1000.0, az_2d, el_2d)
    return zz + radar.altitude['data']


def check_sounding_for_montonic(snd_T, snd_z):
    """
    So the sounding interpolation doesn't fail, force the sounding to behave
    monotonically so that z always increases. This eliminates data from
    descending balloons.
    """
#    snd_T = sounding.soundingdata['temp']  # In old SkewT, was sounding.data
#    snd_z = sounding.soundingdata['hght']  # In old SkewT, was sounding.data
    dummy_z = []
    dummy_T = []
#    if not snd_T.mask[0]: #May cause issue for specific soundings
    if True:
        dummy_z.append(snd_z[0])
        dummy_T.append(snd_T[0])
        for i, height in enumerate(snd_z):
            if i > 0:
                if snd_z[i] > snd_z[i-1]:
                    dummy_z.append(snd_z[i])
                    dummy_T.append(snd_T[i])
        snd_z = np.array(dummy_z)
        snd_T = np.array(dummy_T)
    return snd_T, snd_z


def interpolate_sounding_to_radar(snd_T, snd_z, radar):
    """Takes sounding data and interpolates it to every radar gate."""
    radar_z = get_z_from_radar(radar)
    radar_T = None
    snd_T, snd_z = check_sounding_for_montonic(snd_T, snd_z)
    shape = np.shape(radar_z)
    rad_z1d = radar_z.ravel()
    rad_T1d = np.interp(rad_z1d, snd_z, snd_T)
    return np.reshape(rad_T1d, shape), radar_z


def extract_unmasked_data(radar, field, bad=-32768):
    """Simplify getting unmasked radar fields from Py-ART"""
    return radar.fields[field]['data'].filled(fill_value=bad)


def add_field_to_radar_object(field, radar, field_name='FH', units='unitless', 
                              long_name='Hydrometeor ID', standard_name='Hydrometeor ID',
                              dz_field='ZC'):
    """
    Adds a newly created field to the Py-ART radar object. If reflectivity is a masked array,
    make the new field masked the same as reflectivity.
    """
    fill_value = -32768
    masked_field = np.ma.asanyarray(field)
    masked_field.mask = masked_field == fill_value
    if hasattr(radar.fields[dz_field]['data'], 'mask'):
        setattr(masked_field, 'mask', 
                np.logical_or(masked_field.mask, radar.fields[dz_field]['data'].mask))
        fill_value = radar.fields[dz_field]['_FillValue']
    field_dict = {'data': masked_field,
                  'units': units,
                  'long_name': long_name,
                  'standard_name': standard_name,
                  '_FillValue': fill_value}
    radar.add_field(field_name, field_dict, replace_existing=True)
    return radar






def slant_height_indices(x, y, z, ranges, heights):

    xmesh, ymesh = np.meshgrid(x, y)
    dist = np.sqrt(xmesh**2 + ymesh**2)

    z1 = np.argmin(np.abs(z-heights[0]))
    z2 = np.argmin(np.abs(z-heights[1]))
    z3 = np.argmin(np.abs(z-heights[2]))
    z4 = np.argmin(np.abs(z-heights[3]))
    z5 = np.argmin(np.abs(z-heights[4]))

    z_ind_array = np.zeros(xmesh.shape, int)

    z_ind_array[(dist>=ranges[0]) & (dist<=ranges[1])] = z1
    z_ind_array[(dist>ranges[1]) & (dist<=ranges[2])] = z2
    z_ind_array[(dist>ranges[2]) & (dist<=ranges[3])] = z3
    z_ind_array[(dist>ranges[3]) & (dist<=ranges[4])] = z4
    z_ind_array[(dist>ranges[4])] = z5

    return z_ind_array

def radar_slant_values(x, y, z, val_arr, ranges=[0.0, 30.0, 60.0, 80.0, 120.0], heights=[0.5, 1.0, 1.5, 2.0, 2.5]):

    z_ind_array = slant_height_indices(x, y, z, ranges, heights)

    slant_arr = np.zeros(z_ind_array.shape, float)

    for ipt in range(z_ind_array.shape[0]):
        for jpt in range(z_ind_array.shape[1]):
            slant_val = val_arr[z_ind_array[ipt, jpt], ipt, jpt]
            # if (not np.fmod(ipt, 5)) and (not np.fmod(jpt, 5)):
            #   print 'ipt: {}, jpt: {}, dist: {}, zind: {}, dbz value: {}'.format(ipt, jpt, dist[ipt, jpt], z_ind_array[ipt, jpt], dbz_val)
            #   time.sleep(0.05)

            slant_arr[ipt, jpt] = slant_val


    return slant_arr



def composite(arr, bad_val=-9.9, func=np.max):

    comp_arr = func(arr.filled(bad_val), axis=0)
    comp_arr[np.isnan(comp_arr)] = bad_val
    return comp_arr


def rainfall_accumulation(rainarr, tdiffs, minute_thresh):
    #print tdiffs
    good_files = np.where( (tdiffs >= -1.0*minute_thresh*60.0) & (tdiffs <= 0) )
    #print good_files
   # Same with the 60 minutes. We could change these values if needed
   # Here we sum up all the files that fall within the certain timeframes
    accum_rainfall = np.ma.masked_invalid(np.sum(rainarr[good_files]*5.0/(60.0), axis=0))
    return accum_rainfall

    pass

def rainfall_age(rainarr, tdiffs, rainfall_thresh=10.0, default_dtime=20.0):
    age = np.zeros((rainarr.shape[1], rainarr.shape[2]), np.float)
    age = np.ma.masked_where(age == 0.0, age)


   # set to 0's here
   #accrain_array = np.array(rainfall)
    enough_rain = np.where(rainarr/(default_dtime/60.0) >= rainfall_thresh)

   # now need to loop thru each entry in enough_rain
    xy_pairs = zip(enough_rain[1], enough_rain[2])

   # Only go thru unique x,y pairs
    uniq_xy_pairs = list(set(xy_pairs))


    for pt in uniq_xy_pairs:
   # pt is an x,y pair
   # find all the times when it rained more than the threshold value at the given point
        #print pt
        all_pts = np.where( (enough_rain[1] == pt[0]) & (enough_rain[2] == pt[1]) )
        t_inds_thresh = enough_rain[0][all_pts[0]]
        #print t_inds_thresh
   # grab the last time
        last_t = t_inds_thresh[-1] # figure out the time index where it was last raining
        #print last_t
   #print t_inds_thresh
   # figure out what file we're looking at here. This logic is done to allow you to run a file
   # that is not the most current and it will only look backwards relative to that file
        #this_file_ind = np.where(radar_times == file_time)[0][0]
   # Here we just multiply the difference in index by the dtime, which will be set to 4 minutes
        this_age = (tdiffs[last_t])/60.0
        #print this_age
       #this_age = (this_file_ind - last_t)*default_dtime
        if this_age <= rainfall_thresh: # only assign a value if it's under the threshold
               # if we didn't have any thresholds, the map would be covered and unintelligible
            age[pt] = -1.0*this_age

    return age

