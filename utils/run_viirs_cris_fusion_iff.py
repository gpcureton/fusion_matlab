#!/usr/bin/env python
# encoding: utf-8
"""
run_viirs_cris_fusion_iff.py

 * DESCRIPTION: Locally executes the IFF on VIIRS SNPP Fusion L1B files.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2017-09-11.
Copyright (c) 2017 University of Wisconsin SSEC. All rights reserved.
Licensed under GNU GPLv3.
"""

import os
from os.path import basename, dirname, abspath, isdir, isfile, exists, join as pjoin
import sys
import shutil
from glob import glob
from copy import copy

import time
from datetime import datetime, timedelta
from subprocess import check_output, check_call, CalledProcessError
import traceback
import logging

# every module should have a LOG object
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
        stream=sys.stdout, level=level,
        format=console_logFormat,
        datefmt=dateFormat)

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


def create_iff_filename(prefix, platform, dt, output_type='HDF4'):

    origin = 'ssec'
    domain = 'fsn'
    suffix = {'HDF4': '.hdf', 'NETCDF': '.nc'}

    begin_date = dt.strftime('d%Y%m%d')
    begin_time = dt.strftime('t%H%M%S')
    processed_time = datetime.now().strftime('c%Y%m%d%H%M%S')

    fn = '{}{}'.format('_'.join([prefix, platform, begin_date, begin_time,
                                 processed_time, origin, domain]),
                       suffix[output_type])
    return fn

def main(args):
    '''
    Required arguments:
        1. input_dir : location of VGEOM and VL1B inputs
        2. output_type: HDF or NETCDF
        3. output_dir : Directory to place output files
    '''
    verbosity = 2
    setup_logging(verbosity)

    input_dir = args[0]
    output_type = args[1]
    output_dir = args[2]

    LOG.info('input_dir = {}'.format(input_dir))
    LOG.info('output_type = {}'.format(output_type))
    LOG.info('output_dir = {}'.format(output_dir))

    # Return a dictionary of input pairs
    viirs_geolocation_files = glob(pjoin(input_dir, 'VGEOM*.nc'))
    viirs_geolocation_files.sort()
    viirs_cris_fused_files = glob(pjoin(input_dir, 'VL1BM*.nc'))
    viirs_cris_fused_files.sort()

    input_dict = {}

    files_slice = slice(None)
    #files_slice = slice(0,1)

    for fusion_file in viirs_cris_fused_files[files_slice]:
        granule_key = '_'.join(basename(fusion_file).split('_')[2:4])
        input_dict[granule_key] = {}
        input_dict[granule_key]['VL1BM'] = fusion_file

    for geo_file in viirs_geolocation_files[:]:
        granule_key = '_'.join(basename(geo_file).split('_')[2:4])
        if granule_key in input_dict.keys():
            input_dict[granule_key]['VGEOM'] = geo_file

    granule_keys = input_dict.keys()
    granule_keys.sort()

    rc = 0

    # Create the output directory
    current_dir = os.getcwd()
    LOG.debug('Current directory is {}'.format(current_dir))
    create_dir(output_dir)

    # Set the granule length
    granule_length = timedelta(minutes=6.)
    LOG.debug('granule_length is {}'.format(granule_length))

    # Get the current environment settings
    env = os.environ

    # Loop through the granule keys...
    for granule_key in granule_keys:

        # Get the datetime object of the input files...
        start_time = datetime.strptime(granule_key, 'd%Y%m%d_t%H%M%S')
        end_time = start_time + granule_length - timedelta(seconds=1)
        output_file = create_iff_filename('IFFSVM', 'snpp', start_time, output_type)
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DUMMY
        #output_file = basename(glob(pjoin(output_dir, output_file.split('_c')[0]) + '*')[0])
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        LOG.info('granule_key = "{}", start_time, end_time = {}, {}, out_fn = {}'.format(granule_key, start_time, end_time, output_file))

        geo_file = input_dict[granule_key]['VGEOM']
        l1b_file = input_dict[granule_key]['VL1BM']

        # Create the IFF from the NASA VIIRS L1b files
        cmd = 'python -m iff2.iff2 -{5} -o {0} --hdf4 npp viirs-nasa svm {1:%Y%m%d} {1:%H%M%S} {2:%H%M%S} {3} {4}'.format(
                pjoin(output_dir, output_file),
                start_time,
                end_time,
                geo_file,
                l1b_file,
                verbosity * 'v'
                )
        LOG.debug('cmd = "{}"'.format(cmd))

        try:
            rc = check_call([cmd], shell=True, env=env)
            LOG.debug('check_call() return value: {}'.format(rc))
        except CalledProcessError as err:
            rc = err.returncode
            LOG.error("iff2.iff2 returned a value of {}".format(rc))
            return rc

        # Apply the land water mask to the new IFF files
        cmd = 'python -m demlw.ifflw -{1} {0}'.format(
                pjoin(output_dir, output_file),
                verbosity * 'v'
                )
        LOG.debug('cmd = "{}"'.format(cmd))

        # Location of the demlw data location
        demlw_dir = '/mnt/software/flo/iff/landwater'

        try:
            env['DEMLW_DIR'] = demlw_dir
            env['LW_CACHE_FILE'] = pjoin(demlw_dir, 'lw-cache.npy')
            rc = check_call([cmd], shell=True, env=env)
            LOG.debug('check_call() return value: {}'.format(rc))
        except CalledProcessError as err:
            rc = err.returncode
            LOG.error("demlw.ifflw returned a value of {}".format(rc))
            return rc

    #time.sleep(1.)

    return rc


if __name__ == '__main__':
    args = sys.argv[1:]

    usage = \
        """
        Usage: 'run_viirs_cris_fusion_iff.py <options>
        where <options> are:

            1. input_dir : location of VGEOM and VL1B inputs
            2. output_type: HDF4 or NETCDF
            3. output_dir : Directory to place output files
        """

    if len(args) != 3:
        print(usage)

    sys.exit(main(args))
