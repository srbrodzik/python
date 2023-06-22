import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
from copy import deepcopy


#speccmap = plt.cm.get_cmap("Spectral_r", 5) #generate a jet map with 10 values
#spec_vals = speccmap(np.arange(5)) #extract those values as an array

# test line


ncolors = 12
oldcmap = plt.cm.get_cmap("Spectral_r", ncolors) #generate a jet map with 10 values
vals = oldcmap(np.arange(ncolors)) #extract those values as an array

#vals[0] = [0.0, 0.4, 0.0, 1]
vals[2] = [0.8, 0.8, 0.8, 1] #change the first value
vals[3] = [0.8, 0.8, 0.8, 1] #change the first value
#vals[4] = [0.8, 0.8, 0.8, 1] #change the first value
#vals[2] = [1.0, 1.0, 0.25, 1] #change the first value
#vals[3] = [1.0, 0.5, 0.0, 1]
#vals[4] = [0.7,  0.0,  0.0,  1. ]

newcmap = mpl.colors.LinearSegmentedColormap.from_list("new", vals)



plt.imshow(np.random.rand(18,20)*3 - 1, cmap=newcmap, vmin=-1, vmax=3, interpolation="nearest")
plt.colorbar() 
plt.show()
