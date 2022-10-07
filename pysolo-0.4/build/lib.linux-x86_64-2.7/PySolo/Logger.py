"""
  This class creates logger for debug.
"""

import logging
import os

class Logger():
    """ create logger """
    TASK_LOGGER = 'NCAR.PySolo.Logger'
    CONSL_LEVEL_RANGE = range(0, 100)
    LOG_FILE = '~/pysolo.log'
    FORMAT_STR = '%(asctime)s %(levelname)s %(message)s'

    def __init__(self,level):
        ## Create logger
        logger = logging.getLogger(self.TASK_LOGGER)
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(self.FORMAT_STR)
        ## Add FileHandler
        fn = os.path.expanduser(self.LOG_FILE) ## locate to user home directory
        fh = logging.FileHandler(fn)
        fh.name = 'File Logger'
        fh.level = logging.WARNING
        fh.formatter = formatter
        logger.addHandler(fh)

        ## Add (optionally) ConsoleHandler
        if level is not None:
            ch = logging.StreamHandler()
            ch.name = 'Console Logger'
            ch.level = level
            ch.formatter = formatter
            logger.addHandler(ch)

        logger.debug('DEBUG')
        logger.info('INFO')
        logger.warning('WARNING')
        logger.critical('CRITICAL')
