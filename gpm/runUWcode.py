#!/usr/local/python/bin/python -O
"""
:Program:         runUWcode.py
:Author:          Erich Franz Stocker
:System:          PPS special utilities
:Date-Written:    28 September 2014
:Synopsis:
    This mondule contains the code that will run both the UW matlab code for the
    GPM KU L2 analysis and is a modification of a similar program that was 
    used to do analysis for TRMM PR data.

    The following activities are carried out:

    #. determine whether the date was sent in on the command line or is to be
       read from the date file
    #. if required read the date file to get process date
    #. Fix the date for which matlab is to be run
    #. The matLab code is run by day
    #. update the date file to the next day if the date was read from the
      date file rather than sent in on the commandline

    
    Usage: runUWcode.py uw.parm [YYYY MM]

    The current approach uses the matlab program to create daily gzipped netCDF
    files over the specified area of interest.  All this handled by the Matlab
    program.

:copyright: All rights to software owned by NASA/GSFC code 610.2

:Version-change: (30 September 2014) changed the run idl function so it is 
                 called on a month whose last day was processed by the matlab
                 portion.

:Version-change: (1 September 2016) changed to run KU gpm code rather than the
                 TRMM code. Also added DAY as flag to indicate whether 
                 processing starts at beginning of month or some other day in 
                 the month. If the constant flag DAY is **null** or **none**
                 then starts on 1st of month otherwise its value determines
                 the date to start.

:Version-change: (3 September 2016) added feature to allow the year and the
                 month to be processed to be sent in on the commandline. If
                 sent in then the monthly file *uw.date* is neither read or
                 updated. This allow program to be called from runyear.py.

:Version-change: (8 September 2016) removed all vestiges of the IDL run that
                 was done for TRMM but is not required for GPM.
"""
import os
import sys
import logging
import multiprocessing as mp

from dateclass import *
from ppsmpbeoexec import *
from keywordvariable import *

#: module version
VERS='2.2.4'
#: module version date
VERSDATE='8 September 2016'

#:  Day start flag: if DAY not null then the start day is DAY rather than 01
#DAY='08'
DAY='null'

#############################################################
# Setup of the module root logger. This includes the format #
# that is used by others that do not specify their own.     #
#############################################################
FORMAT = '%(asctime)s:[%(levelname)-5s] %(message)s'
logging.basicConfig(format=FORMAT,level=logging.DEBUG) #module logger format
rlog=logging.getLogger('')   #module logger
#rlog.setLevel(logging.INFO)
rlog.setLevel(logging.DEBUG)

class MatLab(KeywordVariable):
    """
    [UWMatLab] class contains the data variables and methods to create the UW
    products. The matLab code and IDL code to create the special products were
    provided by UW and modified by PPS to use compiled code rather than 
    interactive matLab and IDL

    **Many of the key instance variables inherited from KeywordVariable parent**
    """
    def __init__(self,parmIn,year=0,month=0):
        """
        (UWMatLab) constructor inherits from the parent for the instance and
        sets up additional instance variable. Most important job is to call the
        method that sets the starting runDate for the instance.
 
        :param string parmIn: fully qualified name of the input parameter
                              file
        
        :ivar string parmName: copy of the parmIn
        :ivar string year: copy of the year to process
        :ivar string month: copy of the month to process
        :ivar object runDate: Date object representing the start date to run 

        :note: if year and month parameters are default or year sent in as 0
               then it is determined the date file should be used for establishing
               the *runDate* object otherwise if year not zero then year and 
               month sent in are used to create the *runDate* object.

        :raises: general exception, keywordvariable exception or IOException
                 which it passes to instantiator
        """
        self.parmName=parmIn
        self.year=year
        self.month=month
        self.retrievalMethod='' #default value if no inheritance
        try:
            KeywordVariable.__init__(self,self.parmName)
        except KeywordVariableException as e:
            raise

        try:
            self.runDate=self.readDateFile()
        except (IOError,Exception) as e:
            raise

    def readDateFile(self):
        """
        (readDateFile) reads a text file that has YYYY MM in it and
        returns the Date Object for the first date of the month

        #. read the file which contains year and month on a line
        #. parse the line into separate year and month strings
        #. build a date string that has MM-01-YYYYT00:00:01 format and then
           instantiate a Date object from it

        :returns: Date object that has the first day of the year and month
                  which was in the file. Time set to one second past mid-
                  night. But time not important except for incrementing date.

        :raises: IOError if the date file can't be read or general exception
                 for parsing

        :note: it is assumed that if year sent then valid month also sent
        """
        if self.year==0: #no year/month sent so read from file
            try:
                fd=open(self.dateFile,'r')
                dline=fd.readlines()[0]
                fd.close()
                dline=dline.strip() #remove trailing white space
                year,month=dline.split() 
            except IOError as e:
                raise
            except Exception as e:
                raise
        else:  #use the year and month sent
            month=self.month
            year=self.year

        #########################################################
        # If no day in the DAY flag then start on day 1 of the  #
        # month. Otherwise start on the value of DAY of month.  #
        #########################################################
        if DAY.lower() == 'null' or DAY.lower() == 'none':        
            datestr='%02d-01-%04dT00:00:01' %(int(month),int(year))
        else:
            datestr='%02d-%02d-%04dT00:00:01' %(int(month),int(DAY),int(year))

        return Date(datestr)

    def matLab(self,year,month,day):
        """
        (matLab) builds a full matlab command line which would run the matlab,
        provided by UWash, for the specified date. Actually runs a script that
        contains the activation of the MatLab program.
       
        :param string year: year to process
        :param string month: month to process
        :param string day: day to process

        :returns: fully qualifed command line to run the command
        """
        matCmd=self.matExec+' '+year+' '+month+' '+day+' '
        matCmd=matCmd+' '+self.dirOUT

        return  matCmd

    def run(self):
        """
        (run) is the main program that sets up the job and result queues.
        It builds the command line for each day of the process month
        and determines whether the IDL part of the job should be run.

        It then:
           #. fills the job queue
           #. adds a [None,None] entry for each processor to the job queue
           #. job queue entry is a list which has the command as the first
              element and the name of a log file as second element
           #. starts a process to handle an item on the queue for each
              processor
           #. It loops over the day by incrementing a Date object by 
              one day each iteration of the loop. It creates a log file to which
              the output of the matlab programs will be redirected. This is one
              log file per day.
           #. waits for all the jobs to finish

        The key is to load the job queue with the commands that will then
        be handled by the processes. The job queue is what is made available
        to each process for handling until it is empty. The flag value
        None indicates empty. There must be a None for each process created
        to take off the queue

        The result queue contains any messages that were encountered on
        the running of a job on a beowulf node. If not empty then 
        something happened.

        :warning: Currently because of limitations in matlab no more than 15
                  processes can safely be handled.

        :returns: 111 if result queue isn't empty otherwise 0
        """
        self.jobs=PPSmpBeoexec.PPSBeoQueue()
        self.results=PPSmpBeoexec.PPSBeoQueue()
        ############################################
        # adjust for leap year as the daysIn Month #
        # always returns 29 for Feb.               #
        ############################################
        daysToRun=self.runDate.daysInMonth()
        rlog.debug('%s-%d','Days to run',daysToRun)
        rlog.debug('%s-%d','Leap year',self.runDate.isLeapYear())
        rlog.debug('%s-%d','Month number',self.runDate.getMonth())

        rlog.debug(daysToRun)
        for day in range(daysToRun):
            if day ==0:
                oneDay=self.runDate
            else:
                oneDay=self.runDate>>day #increment by number days from start
            dstr=oneDay.isoDate()
            rlog.info('date run '+oneDay.isoDate())
            yrstr=dstr[0:4]
            mostr=dstr[5:7]
            daystr=dstr[8:10]
            matCmd=self.matLab(yrstr,mostr,daystr)
            rlog.info(matCmd)
            matLog=self.logDir+'/uwlog.'+oneDay.isoDate()+'.log'
            matLog=matLog.replace(':','')
            self.jobs.put([matCmd,matLog])  #fill job queue with jobs

        #############################################
        # Put end of job markers for each process   #
        # on the job queue to ensure end of process.#
        #############################################
        for nodes in range(self.numProcessors):
            self.jobs.put([None,None])

        ################################################
        # Allocate a process for the number of allowed #
        # processes as provided in parameter file.     #
        ################################################
        processes=[]
        rtype=self.retrievalMethod
        for nodes in range(self.numProcessors):
            p=PPSmpBeoexec((self.sleepTime,self.jobs,self.results,rtype,))
            processes.append(p)
            p.start()
            time.sleep(self.sleepTime)
        ################################################
        # Wait for the all the processes in the list to#
        # finish, accept a end of job token and rejoin.#
        ################################################
        while True:
            if not mp.active_children():
                break
        
        if self.results.empty():
            return 0
        else:
            for line in self.results.get(False,5):
                print 'Error: '+str(line)
            return 111

    def displayError(self):
        """
        (displayError) displays the contents of the results queue.
        """
        rlog.debug('In display Error')
        cnt=1
        while not self.results.empty():
            rlog.error('('+str(cnt)+') '+str(self.results.get()))
            cnt+=1
         
    def updateDateFile(self):
        """
        (updateDateFile) puts the next year and month into the date file. It
        creates a new date object that adds the daysInMonth to the runDate and
        then just enters the year and month into the date file.
 
        :note: does not update if the user sent year/month on command line
        """ 
        if self.year != 0: #user sent year/month
            return
        newDate=self.runDate>>self.runDate.daysInMonth()
        year="%04d" %(newDate.getYear())
        month="%02d" %(newDate.getMonth())
        try:
            fd=open(self.dateFile,'w')
            fd.write(year+' '+month+'\n')
            fd.close()
        except IOError as e:
            raise

if __name__ == '__main__':
    status=0
    argc=len(sys.argv)
    if argc <2:
        rlog.error('usage: runUWCode.py parameter_file [YYYY MM]')
        status=1
    else:
        if argc == 4:
            year=int(sys.argv[2])
            month=int(sys.argv[3])
        else:
            year=0
            month=0    
        try:
            uw=MatLab(sys.argv[1],year,month)
            status = uw.run()
            if status==111:
                uw.displayError()
            else:
                uw.updateDateFile()
        except KeywordVariableException as e:
            rlog.error(e.getErrMsg()+' '+e.getErrArg())
            status=e.getErrNumber()
        except IOError as e:
            rlog.error(str(e))
            status=2
        except Exception as e:
            rlog.error(str(e))
            status=3

    sys.exit(status)       
