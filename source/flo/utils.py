#!/usr/bin/env python
# encoding: utf-8
"""
utils.py

 * DESCRIPTION: This is the collection of utilities for parsing text, running external processes and
 binaries, checking environments and whatever other mundane tasks aren't specific to this project.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2016-04-11.
Copyright (c) 2011 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os
import sys
import string
import traceback
import logging
import time
from glob import glob
from copy import copy
import shutil
from os.path import basename, dirname, abspath, isdir, isfile, exists, join as pjoin
import fileinput

from datetime import datetime

LOG = logging.getLogger(__name__)

def setup_logging(verbosity):
    LOG.debug("Verbosity is {}".format(verbosity))

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[verbosity if verbosity < 4 else 3]

    # Set up the logging
    #console_logFormat = '%(asctime)s : %(name)-12s: %(levelname)-8s %(message)s'
    #console_logFormat = '%(levelname)s:%(name)s:%(msg)s') # [%(filename)s:%(lineno)d]'
    #console_logFormat = '%(asctime)s : %(funcName)s:%(lineno)d:  %(message)s'
    #console_logFormat = '%(message)s'
    console_logFormat = '(%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:  %(message)s'
    #console_logFormat = '%(asctime)s : (%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:' \
        #' %(message)s'

    dateFormat = '%Y-%m-%d %H:%M:%S'
    logging.basicConfig(
        stream=sys.stderr,
        level=level,
        format=console_logFormat,
        datefmt=dateFormat)

def _replaceAll(intputfile, searchExp, replaceExp):
    '''
    Replace all instances of 'searchExp' with 'replaceExp' in 'intputfile'
    '''
    for line in fileinput.input(intputfile, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp, replaceExp)
        sys.stdout.write(line)
    fileinput.close()

def cleanup(work_dir, objs_to_remove):
    """
    cleanup work directiory
    remove evething except the products
    """
    for file_obj in objs_to_remove:
        try:
            if isdir(file_obj):
                LOG.debug('Removing directory: {}'.format(file_obj))
                shutil.rmtree(file_obj)
            elif isfile(file_obj):
                LOG.debug('Removing file: {}'.format(file_obj))
                os.unlink(file_obj)
        except Exception:
            LOG.warn(traceback.format_exc())

def link_files(dest_path, files):
    '''
    Link ancillary files into a destination directory.
    '''
    files_linked = 0
    for src_file in files:
        src = basename(src_file)
        dest_file = pjoin(dest_path, src)
        if not exists(dest_file):
            LOG.debug("Link {0} -> {1}".format(src_file, dest_file))
            os.symlink(src_file, dest_file)
            files_linked += 1
        else:
            LOG.warn('link already exists: {}'.format(dest_file))
            files_linked += 1
    return files_linked

class SipsEnvironment(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

def make_time_stamp_d(timeObj):
    '''
    Returns a timestamp ending in deciseconds
    '''
    dateStamp = timeObj.strftime("%Y-%m-%d")
    # seconds = repr(int(round(timeObj.second + float(timeObj.microsecond) / 1000000.)))
    deciSeconds = int(round(float(timeObj.microsecond) / 100000.))
    deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
    timeStamp = "{}.{}".format(timeObj.strftime("%H:%M:%S"), deciSeconds)
    return "{} {}".format(dateStamp, timeStamp)

def make_time_stamp_m(timeObj):
    """
    Returns a timestamp ending in milliseconds
    """
    dateStamp = timeObj.strftime("%Y-%m-%d")
    # seconds = repr(int(round(timeObj.second + float(timeObj.microsecond) / 1000000.)))
    milliseconds = int(round(float(timeObj.microsecond) / 1000.))
    milliseconds = repr(000 if milliseconds > 999 else milliseconds)
    timeStamp = "{}.{}".format(timeObj.strftime("%H:%M:%S"), str(milliseconds).zfill(3))
    return "{} {}".format(dateStamp, timeStamp)

def execution_time(startTime, endTime):
    '''
    Converts a time duration in seconds to days, hours, minutes etc...
    '''

    time_dict = {}

    delta = endTime - startTime
    days, remainder = divmod(delta, 86400.)
    hours, remainder = divmod(remainder, 3600.)
    minutes, seconds = divmod(remainder, 60.)

    time_dict['delta'] = delta
    time_dict['days'] = int(days)
    time_dict['hours'] = int(hours)
    time_dict['minutes'] = int(minutes)
    time_dict['seconds'] = seconds

    return time_dict

def create_dir(dir):
    '''
    Create a directory
    '''
    returned_dir = copy(dir)
    LOG.debug("We want to create the dir {} ...".format(dir))

    try:
        if returned_dir is not None:
            returned_dir_path = dirname(returned_dir)
            returned_dir_base = basename(returned_dir)
            LOG.debug("returned_dir_path = {}".format(returned_dir_path))
            LOG.debug("returned_dir_base = {}".format(returned_dir_base))
            returned_dir_path = '.' if returned_dir_path=="" else returned_dir_path
            LOG.debug("returned_dir_path = {}".format(returned_dir_path))
            # Check if a directory and has write permissions...
            if not exists(returned_dir) and os.access(returned_dir_path, os.W_OK):
                LOG.debug("Creating directory {} ...".format(returned_dir))
                os.makedirs(returned_dir)
                # Check if the created dir has write permissions
                if not os.access(returned_dir, os.W_OK):
                    msg = "Created dir {} is not writable.".format(returned_dir)
                    raise SipsEnvironment(msg)
            elif exists(returned_dir):
                LOG.debug("Directory {} exists...".format(returned_dir))
                if not(isdir(returned_dir) and os.access(returned_dir, os.W_OK)):
                    msg = "Existing dir {} is not writable.".format(returned_dir)
                    raise SipsEnvironment(msg)
            else:
                raise SipsEnvironment("Cannot create {}".format(returned_dir))
    except SipsEnvironment:
        LOG.debug("Unable to create {}".format(returned_dir))
        LOG.debug(traceback.format_exc())
        returned_dir = None
    except OSError:
        LOG.debug(
            "Unable to create new dir '{}' in {}".format(returned_dir_base, returned_dir_path))
        LOG.debug(traceback.format_exc())
        returned_dir = None
    except Exception:
        LOG.warning("General error for {}".format(returned_dir))
        LOG.debug(traceback.format_exc())
        returned_dir = None

    LOG.debug('Final returned_dir = {}'.format(returned_dir))
    return returned_dir
