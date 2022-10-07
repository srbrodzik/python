#!/usr/bin/python3

import sys
import matplotlib.pyplot as plt
import pyart

# SAMPLE INPUTS
# plotPPI_fromCfradial.py cfrad.20180208_192405.000_to_20180208_192658.000_1_SUR.nc DBZH 2

if len(sys.argv) != 4:
    raise SystemExit("Useage: {} {}".format(sys.argv[0], "[file] [field] [sweepNumber]"))
else:
    filename = sys.argv[1]
    field = sys.argv[2]
    sweepNumber = int(sys.argv[3])

print('filename    = ',filename)
print('field       = ',field)
print('sweepNumber = ',sweepNumber)

radar = pyart.io.read_cfradial(filename)
display = pyart.graph.RadarMapDisplay(radar)
display.plot_ppi(field,sweep=sweepNumber)
plt.show()



