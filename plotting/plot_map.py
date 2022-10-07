%matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

fig = plt.figure(figsize=(12,8))
m = Basemap(projection='cyl',resolution=None,llcrnrlat=36.2,urcrnrlat=48,llcrnrlon=-88,urcrnrlon=-66.5,)
m.etopo(scale=0.5, alpha=0.5)
#draw_map(m)

plt.show()

