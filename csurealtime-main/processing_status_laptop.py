# This script makes a GUI with the status of processing


from itertools import cycle
import glob
import argparse
import datetime
import time
from copy import deepcopy
import os
import numpy as np
import psutil

try:
    # Python2
    import Tkinter as tk
    #from Tkinter.font import Font
except ImportError:
    # Python3
    import tkinter as tk
    #from tkinter.font import Font


base_path = '/home/rmet/projects/realtime_radar'
radar_path = '/shares/radar_data/projects/test'
fig_path = '%s/figures'%(base_path)


output_time_fmt = '%Y/%m/%d %H:%M:%S'
file_time_fmt = '%Y%m%d_%H%M%S'

def file_time(fname):
    #print fname
    try:
        done = False
        for i in range(len(fname)):
            try:
                dummy = int(fname[i])
                number_ind = i
                break
            except Exception, e:
                #print e
                pass

    except Exception, fte:
        print 'Error with getting file name: {}'.format(fte)
        print 'file name: {}'.format(fte)


    ftime_string = fname[number_ind:number_ind+15]
    return datetime.datetime.strptime(ftime_string, file_time_fmt)


def status_info(tdiff, vals=[12.0, 30.0]):
    out = {}
    if tdiff > vals[1]:
        out['textcolor'] = 'white'
        out['bgcolor'] = 'red'
    elif ( (tdiff >= vals[0]) and (tdiff <= vals[1]) ):
        out['textcolor'] = 'white'
        out['bgcolor'] = 'gold'
    else:
        out['textcolor'] = 'white'
        out['bgcolor'] = 'green'            

    return out


class App(tk.Tk):
    '''Tk window/label adjusts to size of image'''
    def __init__(self, x=1800, y=1500, update_cycle=60.0, bg_color='white'):
        # the root will be self
        tk.Tk.__init__(self)
        # set x, y position only
        self.time = datetime.datetime(2000, 1, 1)
        #self.geometry('+{}+{}'.format(x, y))
        self.geometry('550x255+2000+1200')
        self.update_cycle = update_cycle
        self.configure(background=bg_color)
        self.text_bg_color = bg_color

#	myFont = Font(family="Times New Roman", size=12)
	#text = tk.Text(self)
	#text.configure(font=myFont)

	self.option_add("*Font", "helvetica")

        self.title('Realtime radar processing status')

        self.text_string = 'OK'
        self.bg_color = 'green'

        self.toprow = 2

        ###########################################################################

        self.time_label = tk.Label(self, text="Time: %s"%(self.time.strftime(output_time_fmt)), bg=self.text_bg_color)
        self.time_label.grid(row=self.toprow-1, column=0)

        ###########################################################################
        self.label_raw = tk.Label(self, text="Raw data", bg=self.text_bg_color).grid(row=self.toprow, column=0)
        self.raw_file_label = tk.Label(self, text="nothing yet")
        self.raw_file_label.grid(row=self.toprow, column=10)

        ###########################################################################

        self.label_ppi_status = tk.Label(self, text="PPI status", bg=self.text_bg_color).grid(row=self.toprow+1, column=0)

        ###########################################################################


        self.label_qc = tk.Label(self, text="PPI QC data", bg=self.text_bg_color).grid(row=self.toprow+2, column=0)
        self.qc_file_label = tk.Label(self, text="nothing yet")
        self.qc_file_label.grid(row=self.toprow+2, column=10)

        ###########################################################################

	self.label_ppirhi = tk.Label(self, text="Polar PPI figure", bg=self.text_bg_color).grid(row=self.toprow+3, column=0)
	self.polar_label = tk.Label(self, text="nothing yet")
        self.polar_label.grid(row=self.toprow+3, column=10)

        ###########################################################################


        ###########################################################################

	self.label_gridded = tk.Label(self, text="Gridded data", bg=self.text_bg_color).grid(row=self.toprow+4, column=0)

        self.gridded_file_label = tk.Label(self, text="nothing yet")
        self.gridded_file_label.grid(row=self.toprow+4, column=10)

        ###########################################################################

	self.label_rain = tk.Label(self, text="Rainfall", bg=self.text_bg_color).grid(row=self.toprow+6, column=0)
	self.rain_label = tk.Label(self, text="nothing yet")
        self.rain_label.grid(row=self.toprow+6, column=10)

        ###########################################################################


        # This one doesn't have a changing value
        self.label_rhi_status = tk.Label(self, text="RHI status", bg=self.text_bg_color).grid(row=self.toprow+8, column=0)

        ###########################################################################



        self.label_rhi_qc = tk.Label(self, text="RHI QC data", bg=self.text_bg_color).grid(row=self.toprow+10, column=0)
        self.rhi_qc_file_label = tk.Label(self, text="nothing yet")
        self.rhi_qc_file_label.grid(row=self.toprow+10, column=10)

        ###########################################################################



        #self.figures_label = tk.Label(self, text="Figures").grid(row=self.toprow+6, column=0)
        #s = tk.Label(self, text="Raw data").grid(row=self.toprow, column=0)
        
        self.label_rhi = tk.Label(self, text="Polar RHI figure", bg=self.text_bg_color).grid(row=self.toprow+11, column=0)
        self.rhi_label = tk.Label(self, text="nothing yet")
        self.rhi_label.grid(row=self.toprow+11, column=10)

        ###########################################################################



        # This is computer stuff down here


        self.label_computer = tk.Label(self, text="Computer info", bg=self.text_bg_color).grid(row=self.toprow+13, column=0)

        self.label_disk_space = tk.Label(self, text="HD space", bg=self.text_bg_color).grid(row=self.toprow+14, column=0)

        self.label_saved_data = tk.Label(self, text="Saved data file size", bg=self.text_bg_color).grid(row=self.toprow+16, column=0)
        self.label_memory = tk.Label(self, text="Memory usage", bg=self.text_bg_color).grid(row=self.toprow+17, column=0)


        self.dspace_label = tk.Label(self, text="nothing yet")
        self.dspace_label.grid(row=self.toprow+14, column=10)        

        # self.sd_label = tk.Label(self, text="nothing yet")
        # self.sd_label.grid(row=self.toprow+16, column=10)
        self.memory_label = tk.Label(self, text="nothing yet")
        self.memory_label.grid(row=self.toprow+17, column=10)


        
        #print type(self.raw_file_label)

        self.show_values()


    # def get_images(self):

    #     image_files = sorted(glob.glob('%s/*.png'%(self.image_path)))

    #     if len(image_files) > self.nimages:
    #         self.image_files = image_files[-1*self.nimages:]
    #         print 'cutting off images'
    #     else:
    #         self.image_files = deepcopy(image_files)
    #         print 'using all images'
    #         pass    
    #     print self.image_files[-1]

    #     self.pictures = [(tk.PhotoImage(file=image), image) for image in self.image_files]
    #     #print self.pictures
    #     # self.picture_display = tk.Label(self)
    #     self.picture_display.pack()


    def secs_since_start(self):
        now = datetime.datetime.utcnow()
        tdiff = (now-self.time).total_seconds() # here is what will get returned

        return tdiff

        
    def show_values(self):
        # time_check = self.secs_since_start()
        # if time_check >= self.update_cycle:
            #self.picture_display(text='Updating')
            #time.sleep(5)
            #self.get_images()
            #print 'Getting new images'
        self.time = datetime.datetime.utcnow() 
        '''cycle through the images and show them'''
        # next works with Python26 or higher
   

        # ft = file_time(os.path.basename(self.pictures[self.counter][1]))
        # #print ft       
        # time_since_file = (datetime.datetime.utcnow() - ft).total_seconds()/60.0

        # if self.counter == self.nimages-1: # if last file, then run a check
        #     if time_since_file >= 30.0:
        #         self.text_string = 'Files are old'
        #         self.bg_color = 'red'
    
        # raw_check = sorted(glob.glob('%s/%04d????/*%04d*'%(radar_path, self.time.year, self.time.year)))
        # raw_base = os.path.basename(raw_check[-1])
        # raw_time = file_time(raw_base)
        # raw_age = (self.time - raw_time).total_seconds()/60.0
        # raw_status = status_info(raw_age)
        # #raw_file_text = os.path.basename(raw_check[-1])
        # #print 'raw file text: {}'.format(raw_file_text)

        # qc_check = sorted(glob.glob('%s/qc_data/ppi/*%04d*'%(radar_path, self.time.year)))
        # qc_base = os.path.basename(qc_check[-1])
        # qc_time = file_time(qc_base)
        # qc_age = (self.time - qc_time).total_seconds()/60.0
        # qc_status = status_info(qc_age)


        # rhi_qc_check = sorted(glob.glob('%s/qc_data/rhi/*%04d*'%(radar_path, self.time.year)))
        # rhi_qc_base = os.path.basename(rhi_qc_check[-1])
        # rhi_qc_time = file_time(rhi_qc_base)
        # rhi_qc_age = (self.time - rhi_qc_time).total_seconds()/60.0
        # rhi_qc_status = status_info(rhi_qc_age)


        # grid_check = sorted(glob.glob('%s/gridded_data/*%04d*'%(radar_path, self.time.year)))
        # grid_base = os.path.basename(grid_check[-1])
        # grid_time = file_time(grid_base)
        # grid_age = (self.time - grid_time).total_seconds()/60.0
        # grid_status = status_info(grid_age)


        polar_check = sorted(glob.glob('%s/raw_ppi/*/*.png'%(fig_path)))
        polar_base = os.path.basename(polar_check[-1])
        polar_time = file_time(polar_base)
        polar_age = (self.time - polar_time).total_seconds()/60.0
        polar_status = status_info(polar_age)
        

        #print '%s/raw_rhi/*/*.png'%(fig_path)
        rhi_check = sorted(glob.glob('%s/raw_rhi/*/*.png'%(fig_path)))
        rhi_base = os.path.basename(rhi_check[-1])
        rhi_time = file_time(rhi_base)
        rhi_age = (self.time - rhi_time).total_seconds()/60.0
        rhi_status = status_info(rhi_age)

        rain_check = sorted(glob.glob('%s/rain/*/*.png'%(fig_path)))
        rain_base = os.path.basename(rain_check[-1])
        rain_time = file_time(rain_base)
        rain_age = (self.time - rain_time).total_seconds()/60.0
        rain_status = status_info(rain_age)
        
        hdinfo = os.statvfs('/home')
        free_space = (hdinfo.f_bavail * hdinfo.f_frsize)/1e9

        free_space_status = {}
	if free_space < 10.0:
	    free_space_status['textcolor'] = 'white'
	    free_space_status['bgcolor'] = 'red'
	elif ( (free_space >= 10.0) and (free_space <= 20.0) ):
	    free_space_status['textcolor'] = 'white'
	    free_space_status['bgcolor'] = 'gold'
	else:
	    free_space_status['textcolor'] = 'white'
	    free_space_status['bgcolor'] = 'green'            



        # saved_data_size = os.path.getsize('%s/saved_data/radar_times_accum.npz'%(base_path))/1.0e6
        # saved_data_status = status_info(saved_data_size, vals=[5.0, 8.0])


        memory_usage = psutil.virtual_memory().percent
        memory_status = status_info(memory_usage, vals=[70.0, 90.0])




        # self.label_ppirhi = tk.Label(self, text="Polar PPI/RHI").grid(row=self.toprow+8, column=0)
        # self.label_rain = tk.L


        self.time_label.config(text='%s'%(self.time.strftime(output_time_fmt)))
        #self.raw_file_label.config(text='%s ----- %d min old'%(raw_base, raw_age), bg=raw_status['bgcolor'], fg=raw_status['textcolor'])
        #self.qc_file_label.config(text='%s ----- %d min old'%(qc_base, qc_age), bg=qc_status['bgcolor'], fg=qc_status['textcolor'])
        #self.rhi_qc_file_label.config(text='%s ----- %d min old'%(rhi_qc_base, rhi_qc_age), bg=rhi_qc_status['bgcolor'], 
#        												fg=rhi_qc_status['textcolor'])
        

        #self.gridded_file_label.config(text='%s ----- %d min old'%(grid_base, grid_age), bg=grid_status['bgcolor'], fg=grid_status['textcolor'])

        self.polar_label.config(text='%s ----- %d min old'%(polar_base, polar_age), bg=polar_status['bgcolor'], fg=polar_status['textcolor'])
        self.rhi_label.config(text='%s ----- %d min old'%(rhi_base, rhi_age), bg=rhi_status['bgcolor'], fg=rhi_status['textcolor'])
        
        self.rain_label.config(text='%s ----- %d min old'%(rain_base, rain_age), bg=rain_status['bgcolor'], fg=rain_status['textcolor'])

        self.dspace_label.config(text='%.2f GB'%(free_space), bg=free_space_status['bgcolor'], fg=free_space_status['textcolor'])

        #self.sd_label.config(text='%.2f MB'%(saved_data_size), bg=saved_data_status['bgcolor'], fg=saved_data_status['textcolor'])
        self.memory_label.config(text='%d %%'%(memory_usage), bg=memory_status['bgcolor'], fg=memory_status['textcolor'])


        self.after(1000, self.show_values)

    def run(self):
        self.mainloop()
# set milliseconds time between slides
delay = 600

# upper left corner coordinates of app window
x = 100
y = 50

# Above are the defaults, below, we'll read in any arguments from the command line
#parser = argparse.ArgumentParser(description='Put in a file to be processed')
if __name__ == '__main__':

    app = App()
    #app.show_values()
    app.run()
