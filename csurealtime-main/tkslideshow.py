import glob
import os
import sys
import time
from copy import deepcopy

import datetime
import argparse

#import Image
#import ImageTk

import Tkinter as tk

# This script is called by: python tkslideshow.py <filepath>
# This code shouldn't need any changes in here, no matter what computer this is running on
# Maybe the Tk would need to be changed depending on the version of Python

def parse_time_string(fname):
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

        return fname[number_ind:] 

    except Exception, be:
    	print 'error with file parsing: {}'.format(be)


def images(nimages=12):
    im = []
    print 'len sys.argv: {}'.format(len(sys.argv))
    print 'sys.argv: {}'.format(sys.argv[1:])
    if len(sys.argv) > 1:
    #	print 'here'
        for path in sys.argv[1:]:
            print 'here, path: {}'.format(path)
            im.extend(images_for(path))
    else:
        im.extend(images_for(os.getcwd()))
    #print 'images: {}'.format(im)
    return sorted(im)
    #return sorted(im, key=lambda s:s.lower())

def images_for(path, nimages=12):
    print 'path: {}'.format(path)
    
    if os.path.isfile(path):
        return [path]
    i = []
    last_file = sorted(glob.glob('%s/*.png'%(path)))[-1]


    print 'last file: {}'.format(last_file)
    last_base = parse_time_string(os.path.basename(last_file))[:15]
    #print 'last base: {}'.format(last_base)

    matches = sorted(glob.glob("%s/*%s*.png"%(path, last_base)))

    if len(matches) > nimages:
    	matches = matches[-1*nimages:]
      #  print 'cutting off images'
    else:
        #matches = deepcopy(matches)
     #   print 'using all images'
        pass    

    #print 'matches: {}'.format(matches)

    for match in matches:
        if match.lower()[-4:] in ('.jpg', '.png', '.gif'):
            i.append(match)

    return i

class App():
    def __init__(self, x=400, y=0, nimages=12, quit_time=120.0):
        self.root = tk.Tk()
        self.root.pack_propagate(False)
        self.root.config(bg="black", width=1000, height=1200)
        self.root.geometry('+{}+{}'.format(x, y))
        self._fullscreen = False
        #self.root.geometry('550x255+2000+1200')
        #self.f_handler(0)
        self.nimages = nimages
        self._images = images(nimages=self.nimages)
        #print '{} images'.format(len(self._images))
        self._image_pos = -1
        self.quit_time = quit_time
        self.start_time = time.time()
        self.root.title('Pop up image viewer: image %d of %d, %d seconds left'%(1, len(self._images), self.quit_time))

        self.root.bind("<Return>", self.return_handler)
        self.root.bind("<space>", self.space_handler)
        self.root.bind("<Escape>", self.esc_handler)
        self.root.bind("<Left>", self.show_previous_image)
        self.root.bind("<Right>", self.show_next_image)
        self.root.bind("q", self.esc_handler)
        self.root.bind("f", self.f_handler)
        self.root.after(100, self.show_next_image)

        self.label = tk.Label(self.root, image=tk.PhotoImage(file=self._images[0]))
        #self.label.configure(borderwidth=2)
        self.label.pack()

        self.set_timer()
        self.root.mainloop()
   
    slide_show_time = 4
    last_view_time = 0
    paused = False
    image = None

    def f_handler(self, e):
        self._fullscreen = not self._fullscreen
        if self._fullscreen:
            self.root.attributes('-fullscreen', True)
        else:
            self.root.attributes('-fullscreen', False)
            self.root.attributes("-zoomed", True)

    def esc_handler(self, e):
        self.root.destroy()

    def return_handler(self, e):
        self.show_next_image()

    def space_handler(self, _):
        self.paused = not self.paused

    def set_timer(self):
        self.root.after(300, self.update_clock)

    def update_clock(self):
    	tsincestart = time.time() - self.start_time
    	#print 'time since start: {}'.format(tsincestart)
    	#self.root.title('Pop up image viewer: %d seconds left'%(self.quit_time-tsincestart))
    	self.root.title('Pop up image viewer: image %d of %d, %d seconds left'%(self._image_pos+1, 
    					len(self._images), self.quit_time-tsincestart))
        if tsincestart > self.quit_time: sys.exit()
        if time.time() - self.last_view_time > self.slide_show_time and not self.paused:
            self.show_next_image()
        self.set_timer()
        self.check_image_size()

    def show_next_image(self, e=None):
        fname = self.next_image()
        #print fname
        if not fname:
            return
        self.show_image(fname)

    def show_previous_image(self, e=None):
        fname = self.previous_image()
        if not fname:
            return
        self.show_image(fname)

    def show_image(self, fname):
        self.original_image = tk.PhotoImage(fname)
        self.image = None
        self.fit_to_box()
        self.last_view_time = time.time()

    def check_image_size(self):
        if not self.image:
            return
        self.fit_to_box()

    def fit_to_box(self):
#        print 'self.image: {}'.format(self.image)
#        if self.image:
#            print 'image'
#            if self.image.size[0] == self.box_width: return
#            if self.image.size[1] == self.box_height: return
        #print self.box_width
        #print self.box_height

        #width, height = self.original_image.size
        #new_size = scaled_size(width, height, self.box_width, self.box_height)
        
        #self.image = self.original_image.resize(new_size, Image.ANTIALIAS)
        #self.image = self.original_image.resize(new_size)
        #self.label.place(x=self.box_width, y=self.box_height, anchor=tk.CENTER)
        self.image = self.original_image
        tkimage = tk.PhotoImage(file=self.image)
        self.label.configure(image=tkimage)
        self.label.image = tkimage

    @property
    def box_width(self):
        return self.root.winfo_width()

    @property
    def box_height(self):
        return self.root.winfo_height()

    def next_image(self):
        if not self._images: 
            return None
        self._image_pos += 1
        self._image_pos %= len(self._images)
        return self._images[self._image_pos]

    def previous_image(self):
        if not self._images: 
            return None
        self._image_pos -= 1
        return self._images[self._image_pos]

def scaled_size(width, height, box_width, box_height):
    source_ratio = width / float(height)
    box_ratio = box_width / float(box_height)
    if source_ratio < box_ratio:
        return int(box_height/float(height) * width), box_height
    else:
        return box_width, int(box_width/float(width) * height)

def test_scaled_size():
    x = scaled_size(width=1871, height=1223, box_width=1920, box_height=1080)
    assert x == (1652, 1080)
    x = scaled_size(width=100, height=100, box_width=1920, box_height=1080)
    assert x ==(1080, 1080)



if __name__ == '__main__':
    app=App()