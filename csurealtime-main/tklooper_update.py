

''' tk_image_slideshow3.py
create a Tkinter image repeating slide show
tested with Python27/33  by  vegaseat  03dec2013

Some modifications by Brody Fuchs, Oct 2017
'''


from itertools import cycle
import glob
import argparse
import datetime
import time
from copy import deepcopy
import os
import numpy as np

try:
    # Python2
    import Tkinter as tk
except ImportError:
    # Python3
    import tkinter as tk

def file_time(fname, prefix='SEA'):
    done = False
    for i in range(len(fname)):
        try:
            dummy = int(fname[i])
            number_ind = i
            break
        except Exception, e:
            #print e
            pass


    ftime_string = fname[number_ind:number_ind+15]
    return datetime.datetime.strptime(ftime_string, '%Y%m%d_%H%M%S')



class App(tk.Tk):
    '''Tk window/label adjusts to size of image'''
    def __init__(self, image_path, x, y, nimages=12, delay=600, loop_delay=2000, update_cycle=20.0):
        # the root will be self
        tk.Tk.__init__(self)
        # set x, y position only
        self.time = datetime.datetime(2000, 1, 1)
        self.geometry('+{}+{}'.format(x, y))
        self.delay = delay
        self.counter = 0
        self.loop_delay = loop_delay
        self.update_cycle = update_cycle
        self.nimages = nimages
        self.image_path = image_path
        #self.get_images()
        self.text_string = 'OK'
        self.bg_color = 'green'



        #self.pictures = cycle((tk.PhotoImage(file=image), image) for image in self.image_files)
        #print self.pictures
        self.picture_display = tk.Label(self)
        self.picture_display.pack()


    def get_images(self):

        image_files = sorted(glob.glob('%s/*/*.png'%(self.image_path)))

        if len(image_files) > self.nimages:
            self.image_files = image_files[-1*self.nimages:]
            print 'cutting off images'
        else:
            self.image_files = deepcopy(image_files)
            print 'using all images'
            pass    
        print self.image_files[-1]

        self.pictures = [(tk.PhotoImage(file=image), image) for image in self.image_files]
        #print self.pictures
        # self.picture_display = tk.Label(self)
        self.picture_display.pack()


    def secs_since_start(self):
        now = datetime.datetime.utcnow()
        tdiff = (now-self.time).total_seconds() # here is what will get returned

        return tdiff

        
    def show_slides(self):
        time_check = self.secs_since_start()
        if time_check >= self.update_cycle:
            #self.picture_display(text='Updating')
            #time.sleep(5)
            self.get_images()
            print 'Getting new images'
            self.time = datetime.datetime.utcnow() 
        '''cycle through the images and show them'''
        # next works with Python26 or higher
        #print datetime.datetime.utcnow()
        #print self.counter
        #img_object, img_name = next(self.pictures)

        ft = file_time(os.path.basename(self.pictures[self.counter][1]))
        #print ft       
        time_since_file = (datetime.datetime.utcnow() - ft).total_seconds()/60.0

        if self.counter == self.nimages-1: # if last file, then run a check
            if time_since_file >= 30.0:
                self.text_string = 'Files are old'
                self.bg_color = 'red'
        
            elif (time_since_file >= 12.0) and (time_since_file < 30.0):
                self.text_string = 'Files are getting old'
                self.bg_color = 'gold'

            else:
                self.text_string = 'OK'
                self.bg_color = 'green'

        self.picture_display.config(image=self.pictures[self.counter][0], text=self.text_string, bg=self.bg_color, fg='white', 
                            font="Helvetica", compound=tk.BOTTOM)
        # shows the image filename, but could be expanded
        # to show an associated description of the image



        self.title('%s (image %d of %d: %d mins old)'%(os.path.basename(self.pictures[self.counter][1]), 
                        self.counter+1, self.nimages, time_since_file))

        
        self.counter += 1
        if self.counter >= self.nimages:
            self.counter = 0 # reset it to 0    
            self.after(self.delay + 1000, self.show_slides)
        else:
            self.after(self.delay, self.show_slides)

    def run(self):
        self.mainloop()
# set milliseconds time between slides
delay = 600
# get a series of gif images you have in the working folder
# or use full path, or set directory to where the images are
#image_files = glob.glob('/home/rmet/SPURS2/realtime_radar/rain_figs/*.png')
# upper left corner coordinates of app window
x = 100
y = 50

# Above are the defaults, below, we'll read in any arguments from the command line
parser = argparse.ArgumentParser(description='Put in a file to be processed')


#parser.add_argument('--file', action="store", dest="file")
#parser.add_argument('--lastfile', action="store", dest="lastfile", type=bool, default=True)
#parser.add_argument('--allfiles', action="store", dest="allfiles", type=bool, default=False)
#parser.add_argument('--realtime', action="store", dest="realtime", type=bool, default=False)
parser.add_argument('--path', action="store", dest="path", default=None)
parser.add_argument('--n_images', action="store", dest="nimages", type=int, default=None)

pargs = parser.parse_args()

#print 'pargs: {}'.format(pargs)

if pargs.path is not None:
    image_path = pargs.path
    #image_files = sorted(glob.glob('%s/*.png'%(image_path)))
else:
    image_path = '/shares/radar_data/projects/test/figures/rain'

# nimages = 12

# if pargs.nimages is not None:
#     nimages = pargs.nimages

# if len(image_files) > nimages:
#     image_files = image_files[-1*nimages:]
#     print 'cutting off images'
# else:
#     print 'using all images'
#     pass    

# print 'nimages: {}'.format(nimages)
#print 'files: {}'.format(image_files)

app = App(image_path, x, y)
app.show_slides()
app.run()
