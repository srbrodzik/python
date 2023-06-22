# import modules that I'm using
import matplotlib
matplotlib.use('TKAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import Tkinter as tk
import tkFileDialog
#from Tkinter import *
import numpy as np
import scipy as sc
import matplotlib as mpl

import pyart
import argparse
from csuram import RadarConfig
from csuram import analysis_tools as atools
import Config
import os
import plot_tools as ptools

    
smooth_window = 50


#Make object for application
class App_Window(tk.Tk):
    def __init__(self, radar_var, rconfig, parent=None):
        tk.Tk.__init__(self, parent)
        #self.radar = radar_object
        self.parent = parent
        self.radar_var = radar_var
        self.rconfig = rconfig

        #self.geometry('+{}+{}'.format(x, y))
        self.geometry('+100+50')
        self.initialize()

    def initialize(self):
    	self.title('Point and click ASCOPE')

        menubar = tk.Menu(self)
        self.config(menu=menubar)
        
        fileMenu = tk.Menu(menubar)       
        
        submenu = tk.Menu(fileMenu)

        submenu.add_command(label="DBZ", command=self.plotDBZ)
        submenu.add_command(label="ZDR", command=self.plotZDR)
        submenu.add_command(label="Phase", command=self.plotPhase)
        submenu.add_command(label="Kdp", command=self.plotKDP)
        submenu.add_command(label="RhoHV", command=self.plotCorr)
        submenu.add_command(label="Vel", command=self.plotVel)


        fileMenu.add_cascade(label='Field', menu=submenu, underline=0)
        
        fileMenu.add_separator()

        fileMenu.add_command(label='Save Image', command=self.saveImg)
        fileMenu.add_command(label='Save Ray Data', command=self.saveRay)

        
        fileMenu.add_command(label="Exit", underline=0, command=self.onExit)
        menubar.add_cascade(label="File", underline=0, menu=fileMenu)   

        self.filename = tkFileDialog.askopenfilename(initialdir="/home/rmet/projects/realtime_radar/", 
            title="Select file", filetypes=(("netcdf files","*.nc"), ("all files","*")))

        self.radar = pyart.io.read(self.filename, file_field_names=True)
        #print self.radar.azimuth['data'].min(), self.radar.azimuth['data'].max(), self.radar.azimuth['data'].shape
        self.scan_type = self.radar.scan_type

               
        self.canvasFig = plt.figure(1)
        self.fig, self.ax = plt.subplots(2, 1, figsize=(6, 9), dpi=100)


	self.display_object = pyart.graph.RadarDisplay(self.radar)


        self.line1, = self.ax[1].plot([], [], color='red', linewidth=0.7, alpha=1.0)
        self.line2, = self.ax[1].plot([], [], color='black', linewidth=1.0, alpha=0.5)
        self.line3, = self.ax[1].plot([], [], color='green', linewidth=1.0, alpha=0.5)
        self.ax[1].grid(True)


        self.fig.subplots_adjust(top=0.95, bottom=0.07, left=0.10, right=0.95, hspace=0.50)


        self.canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.callbacks.connect('button_press_event', self.ascope_plot)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

	self.plotRadarVariable(self.radar_var)


        self.resizable(True,False)
        self.update()


    def plotRadarVariable(self, variable):

    	self.ax[0].clear()

    	try:
    	    self.axcb.clear()
    	except Exception:
    	    pass

    	self.display_object.plot(variable, sweep=0, vmin=self.rconfig.plot_params[variable]['lims'][0], 
        			vmax=self.rconfig.plot_params[variable]['lims'][1], cmap=self.rconfig.plot_params[variable]['cmap'], 
        			colorbar_label=self.rconfig.plot_params[variable]['longname'], mask_outside=False, ax=self.ax[0], 
        			colorbar_flag=False, colorbar_orient='horizontal')

    	if self.scan_type == 'rhi':

            self.ax[0].set_xlim(0, 120)
            self.ax[0].set_ylim(0, 13)

        elif self.scan_type == 'ppi':

            self.ax[0].set_xlim(-120, 120)
            self.ax[0].set_ylim(-120, 120)
            self.ax[0].set_aspect('equal')



        self.ax[1].set_xlabel('Ray range (km)')
        self.ax[1].set_ylabel('%s %s'%(self.rconfig.plot_params[variable]['longname'], self.rconfig.plot_params[variable]['units']), labelpad=0)
       
        self.ax[1].set_xlim(0, 120)
        self.ax[1].set_ylim(self.rconfig.plot_params[variable]['lims'][0], self.rconfig.plot_params[variable]['lims'][1])



	bds = self.rconfig.plot_params[variable]['lims']
	self.cmap_bounds = np.arange(bds[0], bds[1]+ self.rconfig.plot_params[variable]['delta'], self.rconfig.plot_params[variable]['delta'])
	location = [0.05, 0.50, 0.9, 0.03]
	self.axcb = self.fig.add_axes(location) # x pos, y pos, x width, y width
	self.cb = mpl.colorbar.ColorbarBase(self.axcb, cmap=self.rconfig.plot_params[variable]['cmap'], 
					orientation='horizontal', norm=self.rconfig.plot_params[variable]['norm'])
	self.cb.set_label('%s %s' % (self.rconfig.plot_params[variable]['longname'], self.rconfig.plot_params[variable]['units']), labelpad=-1)



    def plotDBZ(self):
    	self.radar_var = self.rconfig.dz_name
    	self.plotRadarVariable(self.radar_var)

    def plotZDR(self):
    	self.radar_var = self.rconfig.zdr_name
    	self.plotRadarVariable(self.radar_var)

    def plotPhase(self):
    	self.radar_var = self.rconfig.ph_name
    	self.plotRadarVariable(self.radar_var)

    def plotKDP(self):
    	self.radar_var = self.rconfig.kdp_name
    	self.plotRadarVariable(self.radar_var)

    def plotCorr(self):
    	self.radar_var = self.rconfig.rho_name
    	self.plotRadarVariable(self.radar_var)

    def plotVel(self):
    	self.radar_var = self.rconfig.vel_name
    	self.plotRadarVariable(self.radar_var)

    def saveImg(self):
	fname = tkFileDialog.asksaveasfilename(initialdir="/home/rmet",title="Save file", 
			filetypes=(("png files","*.png"), ("all files","*.*")))

	plt.savefig(fname)

    def saveRay(self):
	
	fname = tkFileDialog.asksaveasfilename(initialdir="/home/rmet",title="Save file", 
			filetypes=(("numpy save file","*.npz"), ("all files","*.*")))    	

	np.savez(fname, x=np.ma.filled(self.ray_x, np.nan), y=np.ma.filled(self.ray_y, np.nan))


    def onExit(self):
        
        self.quit()


    def refreshLinePlot(self, x, y, scat_x=None, scat_y=None, title=None, lineplot_angle=None, 
    				radius=None, met_angle=None, smooth=None):
    	
        if hasattr(self, 'scat_handle'):
            self.scat_handle.remove()

        if hasattr(self, 'cursor_scat'):
            try:
            	self.cursor_scat.remove()
            except Exception:
            	pass

    	self.line1.set_data(x,y)
    	self.ray_x = x
    	self.ray_y = y
    	if smooth is not None:
    	    xy_smooth = atools.convolve_smooth(x, y, window='121', width=smooth, mode='extend')[0]
    	    self.line2.set_data(x, xy_smooth)


    	#self.line3.set_data(x, y+20.0*np.log10(x) -30.0)

    	if scat_x is not None and scat_y is not None:
    	    self.scat_handle = self.ax[1].scatter([scat_x], [scat_y], c='black', s=100, alpha=0.5, edgecolors='none', zorder=10)

        if len(self.ax[0].lines):
            try:
            	del self.ax[0].lines[-1]
            except Exception:
            	pass


    	if lineplot_angle is not None:
    	    end_of_ray_x_val = np.cos(np.radians(lineplot_angle))*125.0

            self.raytrace_handle = ptools.plot_azimuth(self.ax[0], 0, 0, [lineplot_angle])

    	    #self.raytrace_handle = self.ax[0].plot(rayx, rayy, color='black', alpha=0.4, linewidth=4)

    	if scat_x is not None and lineplot_angle is not None:
    	    sx, sy = ptools.plot_point((0, 0), lineplot_angle, radius)
    	    self.cursor_scat = self.ax[0].scatter([sx], [sy], 
    	    				s=50, facecolors='none', edgecolors='black', alpha=0.6)
    	    pass


    	if title is not None:
    	    self.ax[1].set_title(title)

    	#plt.tight_layout()
	self.canvas.draw()

    def OnButtonClick(self):
        # file is opened here and some data is taken
        # I've just set some arrays here so it will compile alone

        x = np.arange(0, 100)
        y = np.sin(x)

        self.refreshLinePlot(x, y)



    def get_mouse_axis_coords(self, event):
    	if event.inaxes is not None:
    	    x, y = event.xdata, event.ydata
    	    return x, y
    	else:
    	    print 'Clicked outside'



    def get_mouse_fig_coords(self, event):
    	print dir(event)
    	x, y = event.x, event.y
    	print('{}, {}'.format(x, y))
    


    def ascope_plot(self, event):
    	xax, yax = self.get_mouse_axis_coords(event)

    	cur_range = np.sqrt(xax**2 + yax**2)

    	rad_range = self.radar.range['data']

    	if self.scan_type == 'rhi':
    	    cur_el = np.degrees(np.arctan2(yax, xax))
    	    #if cur_el < 0: cur_el += 360.0

    	    #print 'click --> x: {}, y: {}, range: {}, el: {}'.format(xax, yax, cur_range, cur_el)

    	    el_ind = np.argmin(np.abs(cur_el - self.radar.get_elevation(0)))

    	    dbz_ray = self.radar.get_field(0, self.radar_var)[el_ind]
    	    range_ind = np.argmin(np.abs(cur_range - self.radar.range['data']/1000.0))
    	    dbz_pt = dbz_ray[range_ind]
    	    self.refreshLinePlot(rad_range/1000.0, dbz_ray, scat_x=cur_range, scat_y=dbz_pt, 
    				title='%.1f$^{\circ}$ elevation'%(cur_el), lineplot_angle=cur_el, 
    				radius=cur_range, smooth=smooth_window)

    	elif self.scan_type:
    		# need to plot the angle with the physics angle (0 at east) but get the correct ray index with the
    		# meteorological angle, which is 0 at north and is what the radar uses
    	    cur_az = np.degrees(np.arctan2(yax, xax))
    	    cur_az_met = 90.0 - cur_az
    	    if cur_az_met < 0: cur_az_met += 360.0


    	    az_ind = np.argmin(np.abs(cur_az_met - self.radar.get_azimuth(0)))
    	    dbz_ray = self.radar.get_field(0, self.radar_var)[az_ind]
    	    range_ind = np.argmin(np.abs(cur_range - self.radar.range['data']/1000.0))
    	    dbz_pt = dbz_ray[range_ind]
    	    self.refreshLinePlot(rad_range/1000.0, dbz_ray, scat_x=cur_range, scat_y=dbz_pt, 
    				title='%.1f$^{\circ}$ azimuth'%(cur_az_met), lineplot_angle=cur_az, 
    				radius=cur_range, smooth=smooth_window)









if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Put in a file to be processed')
    parser.add_argument('--var', action="store", dest="var", default='DZQC')
    parser.add_argument('--config', action="store", dest="config", default=None)
    pargs = parser.parse_args()

    print 'pargs.config: {}'.format(pargs.config)

    base_path = os.path.dirname(os.path.realpath('__file__'))
    if pargs.config is None:
    	cfg_file = Config.Config('%s/realtime_config.yaml'%(base_path))
    else:
	cfg_file = Config.Config('%s'%(pargs.config))

    ycfg = cfg_file.v
    cfg = RadarConfig.RadarConfig(dz=ycfg['dz_name'], zdr=ycfg['zdr_name'], kdp=ycfg['kdp_name'], vel=ycfg['vel_name'],
                                #rho=ycfg['rho_name'], hid='HID', vel=ycfg['vel_name'])
                                rho=ycfg['rho_name'], hid='HID')

    print cfg.plot_params.keys()






    MainWindow = App_Window(pargs.var, cfg)
    MainWindow.mainloop()


