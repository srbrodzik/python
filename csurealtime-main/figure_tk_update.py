# import modules that I'm using
import matplotlib
matplotlib.use('TKAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import Tkinter
#from Tkinter import *
import numpy as np
import scipy as sc

import pyart
import argparse
from csuram import RadarConfig
import Config



#Make object for application
class App_Window(Tkinter.Tk):
    def __init__(self, display_object, parent=None):
        Tkinter.Tk.__init__(self, parent)
        self.display_object = display_object
        self.parent = parent
        #self.bind("<Button 1>", self.get_mouse_fig_coords)
        self.initialize()

    def initialize(self):
        button = Tkinter.Button(self, text="Open File", command=self.OnButtonClick).pack(side=Tkinter.TOP)
        self.canvasFig = plt.figure(1)
        self.fig, self.ax = plt.subplots(2, 1, figsize=(6, 8), dpi=100)
        #FigSubPlot = Fig.add_subplot(111)

        self.display_object.plot('DZQC', sweep=0, vmin=cfg.lims[cfg.dz_name][0], vmax=cfg.lims[cfg.dz_name][1], cmap=cfg.cmaps[cfg.dz_name], 
                         colorbar_label='Reflectivity (dBZ)', mask_outside=True, ax=self.ax[0], colorbar_flag=False)
        self.ax[0].set_xlim(0, 120)
        self.ax[0].set_ylim(0, 15)

        self.line1, = self.ax[1].plot([], [], 'r-')
        self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.callbacks.connect('button_press_event', self.get_mouse_axis_coords)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
        self.resizable(True,False)
        self.update()


    def refreshLinePlot(self, x, y):
    	self.line1.set_data(x,y)
	#ax = self.canvas.figure.axes[0]
	self.ax[1].set_xlim(x.min(), x.max())
	self.ax[1].set_ylim(y.min(), y.max())        
	self.canvas.draw()

    def OnButtonClick(self):
        # file is opened here and some data is taken
        # I've just set some arrays here so it will compile alone

        x = np.arange(0, 100)
        y = np.sin(x)

        self.refreshLinePlot(x, y)

    def get_mouse_axis_coords(self, event):
    	if event.inaxes is not None:
    	    self.x, self.y = event.xdata, event.ydata
    	    print('{}, {}'.format(self.x, self.y))
    	else:
    	    print 'Clicked outside'

    def get_mouse_fig_coords(self, event):
    	print dir(event)
    	x, y = event.x, event.y
    	print('{}, {}'.format(x, y))
    


if __name__ == "__main__":

    base_path = '/home/rmet/SPURS2/realtime_radar'
    cfg_file = Config.Config('%s/realtime_config.yaml'%(base_path))
    radarfile = '/run/media/rmet/Seagate Expansion Drive/spurs2_backup/realtime_cfradial/rhi/SEA20171025_220254_rhi.nc'

    ycfg = cfg_file.v
    cfg = RadarConfig.RadarConfig(dz=ycfg['dz_name'], zdr=ycfg['zdr_name'], kdp=ycfg['kdp_name'], 
                                rho=ycfg['rho_name'], hid='HID')


    radar = pyart.io.read(radarfile, file_field_names=True)

    display = pyart.graph.RadarDisplay(radar)

    MainWindow = App_Window(display)
    MainWindow.mainloop()


