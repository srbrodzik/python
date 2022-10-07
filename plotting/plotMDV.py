import matplotlib.pyplot as pltimport
import pyart

filename = '/home/disk/bob/olympex/mdv/cart/npol/sur_fromCFradial/20151114/101849.mdv'
radar = pyart.io.read_mdv(filename)
display = pyart.graph.RadarMapDisplay(radar)
fig = plt.figure(figsize=[5, 5])
ax = fig.add_subplot(111)
display.plot_ppi_map('filtered_radar_reflectivity', 0, vmin=-16., vmax=64)
plt.show()


# something similar seen here:
# http://arm-doe.github.io/pyart/dev/auto_examples/plotting/plot_ppi_mdv.html
