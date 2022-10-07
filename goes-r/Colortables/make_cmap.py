# Takes in RGB triples where values of R, G and B range from 0 to 255.
# For each RGB triple there is a corresponding position on the colorbar.
# This function will nterpolate the colors in between the values you choose.
# You must first know the intended range of your colorbar.

# Sample inputs for popular colormap used by NWS
# colors = [(255,255,255), (255,255,0), (140,25,115), (255,125,200),(255,255,255), (0,0,0), \
#           (255,0,0), (255,255,0), (0,255,0), (0,0,115), (0,255,255), (200,200,200), (0,0,0)]
# position = [0, .12, .120001, .18, .180001, .25, .29, .35, .41, .49, .55, .550001, 1]
# my_cmap = make_cmap(colors, position=position, bit=True)

def make_cmap(colors, position=None, bit=False):
    
    import matplotlib as mpl
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap
