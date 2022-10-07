# Calculating the position of a point on a map
# When inverse is False, the input is a point in lon/lat and the output is the
#    point in map coordinates.  When inverse is True, the behavior is the
#    opposite.
# See https://basemaptutorial.readthedocs.io/en/latest/basic_functions.html

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

map = Basemap(projection='aeqd', lon_0 = 10, lat_0 = 50)
plt.show()

print map(10, 50)
print map(20015077.3712, 20015077.3712, inverse=True)

