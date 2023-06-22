# Just a collection of functions for random string manipulation and reading

import numpy as np 
import sys





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