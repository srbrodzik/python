import os
import subprocess as sp 
import numpy as np
import datetime
import glob


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

def convert_decimal_to_degree(val):
    val_string = str(val)
    vs_split = val_string.split('.')
    decimal_string = vs_split[1]
    if len(decimal_string) == 1:
        decimal_string += '0'
    decimal_float = float(decimal_string)*60.0/100.0
    out_string = "%d$^{\circ}$ %d'"%(int(vs_split[0]), int(decimal_float))
    return out_string


def file_time_location(time1, time_arr):
    out = None  
    for i in range(len(time_arr)-1):
        if time1>time_arr[i] and time1<time_arr[i+1]:
            out = i

    return out


def recent_file(path, before_wildcard='', after_wildcard=''):

    print 'trying to find recent file'
    #p1 = sp.Popen(['ls', '-th', '%s/%s*%s'%(path, before_wildcard, after_wildcard)], stdout=sp.PIPE)
    p1 = sp.Popen(['ls', '-th', '%s/*'%(path)], stdout=sp.PIPE)
    p2 = sp.Popen(['head', '-1'], stdin=p1.stdout, stdout=sp.PIPE)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    out, err = p2.communicate()
    print 'end of trying to find recent file'

    fname = out.strip()
    return fname

def recent_file_glob(path, before_wildcard='', after_wildcard='', index=-1):
    return sorted(glob.glob('%s/%s*%s'%(path, before_wildcard, after_wildcard)), key=os.path.getctime)[index]


def substring_in_string_list(string, str_list):
    out = False
    for sl in str_list:
        if string in sl:
            out = True
            break
    return out


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



def get_valid_files(filenames, t):
    file_bases = np.array([os.path.basename(_) for _ in filenames])
    #print file_bases[:5]
    file_times = np.array([datetime.datetime.strptime(_[len(prefix1):-1*(len(ftype1)+1)], dt_fmt) for _ in file_bases])

    file_time_diffs = np.array([(t - _).total_seconds() for _ in file_times])

    # these are the differences between the current (or specified time) and the file times in seconds
    # We want to see if there is a file that is greater than 0, but is also less than some threshold value

    vfi = np.where( (file_time_diffs >= 0) & (file_time_diffs < tdiff_thresh))[0]
    return vfi



def is_process_running():

    p1 = sp.Popen(['ps', '-eaf'], stdout=sp.PIPE)
    p2 = sp.Popen(['grep', 'realtime_driver.py'], stdin=p1.stdout, stdout=sp.PIPE)
    p3 = sp.Popen(['grep', 'python'], stdin=p2.stdout, stdout=sp.PIPE)
    out, err = p3.communicate()
        # p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        # out, err = p2.communicate()
    return out


def grab_files_time_sorted(path_string):
    searchedfile = glob.glob(path_string)
    files = sorted(searchedfile, key=lambda file: os.path.getctime(file))
    return files

def grab_recent_file(path_string):
    files = grab_files_time_sorted(path_string)
    return files[-1]
