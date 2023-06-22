import numpy as np 
import matplotlib.pyplot as plt 
import os
import glob
import gentools as gtools
import matplotlib as mpl



km2nm = 0.54 # 1 km is 0.54 nautical miles
nm2km = 1.0/km2nm # Just to have both conversions


def modify_cmap(orig_cb_string, modify_dict, ncolors=8, new_name='newcmap'):
    
    oldcmap = plt.cm.get_cmap(orig_cb_string, ncolors) #generate a jet map with 10 values
    old_vals = oldcmap(np.arange(ncolors)) #extract those values as an array
    for k, v in modify_dict.iteritems():
	old_vals[k] = v
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(new_name, old_vals)
    return new_cmap


def save_figure(fig, path, fname, **kwargs):
    """
	This writes the figure to the specified path, but the wrinkle here is that it will check
	for the path you've given it, and if it doesn't exist, it will make it. This way if the code
	is running on a new system, it will create the appropriate directories without the user
	having to go in to make the new directories all the time.

    """
    # first check if the specified path already exists.
    path_check = glob.glob(path)
    if len(path_check) == 0:
	print 'making new directory %s because it does not exist'%(path)
        b = os.system('mkdir -v %s'%(path))
        #print b
        if b != 0:

        # means we need to step back and make the parent directory first before making the subdirectory
            os.system('mkdir -v %s'%(os.path.dirname(path)))
            os.system('mkdir -v %s'%(path))
            pass
    else:
	#print 'directory already there, just gonna save the figure'
	pass

    fig.savefig('%s/%s'%(path, fname), **kwargs)



def plot_point(point, angle, length):
     '''
     point - Tuple (x, y)
     angle - Angle you want your end point at in degrees.
     length - Length of the line you want to plot.

     Will plot the line on a 10 x 10 plot.
     '''

     # unpack the first point
     x, y = point

     # find the end point
     endy = length * np.sin(np.radians(angle))
     endx = length * np.cos(np.radians(angle))

     return x+endx, y+endy


def plot_range_rings(a, rlat, rlon, rads, minor_rads=None, conv=nm2km/111.0, text_flag=True):


    for icr in rads:

        icr_conv = icr*conv
	dummycircle = plt.Circle((rlon, rlat), icr_conv, facecolor='none', edgecolor='black', 
	    				linewidth=0.5, alpha=0.5, linestyle='dashed', zorder=9)
	a.add_artist(dummycircle)
	#if icr_conv < lat_width/2.0:
        if text_flag:
	    dummy_text = a.text(rlon, rlat-icr_conv, '%d nm'%(icr), ha='center', va='center', alpha=0.5, zorder=10)

    if minor_rads is not None:        

        for imcr in minor_rads:

	    imcr_conv = imcr*nm2km/111.0
	    dummycircle_minor = plt.Circle((rlon, rlat), imcr_conv, facecolor='none', edgecolor='black', 
	    				linewidth=0.5, alpha=0.25, linestyle='dashed', zorder=9)
	    a.add_artist(dummycircle_minor)
	    #if imcr_conv < lat_width/2.0:
            if text_flag:
                a.text(rlon, rlat-imcr_conv, '%d nm'%(imcr), ha='center', va='center', alpha=0.25, zorder=10)




def plot_azimuths(a, rlat, rlon, azs, minor_azs=None, color='black', linewidth=0.5, alpha=0.5, linestyle='dashed', conv=nm2km/111.0, max_range=120.0):
    for az in azs:
	fx, fy = plot_point((rlon, rlat), az, 200.0/111.0)
	a.plot([rlon, fx], [rlat, fy], color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyle)

    if minor_azs is not None:
    	for maz in minor_azs:
	    mfx, mfy = plot_point((rlon, rlat), maz, 200.0/111.0)
	    a.plot([rlon, mfx], [rlat, mfy], color=color, linewidth=linewidth, alpha=alpha/3.0, linestyle=linestyle)



def plot_azimuth(a, rlat, rlon, azs, minor_azs=None, color='black', linewidth=4, alpha=0.5, linestyle='solid', radius=125.0):
    for az in azs:
        fx, fy = plot_point((rlon, rlat), az, radius)
        hand = a.plot([rlon, fx], [rlat, fy], color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyle)

    return hand



def ax_labels_minutes(ax):
    axt = ax.get_xticks().tolist()
    for j in range(len(axt)):
	axt[j] = gtools.convert_decimal_to_degree(axt[j])
    ax.set_xticklabels(axt)

    axty = ax.get_yticks().tolist()
    #print 'y tick labels: {}'.format(axty)
    for j in range(len(axty)):
	axty[j] = gtools.convert_decimal_to_degree(axty[j])
    ax.set_yticklabels(axty)


def rotate_ticklabs(ax, x=True, y=True, rotation=25):

    if x:
	for label in ax.get_xticklabels():
	    label.set_rotation(rotation)

    if y:
	for label in ax.get_yticklabels():
	    label.set_rotation(rotation)

