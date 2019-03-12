#!/usr/bin/env python
# encoding: utf-8
"""
Script to global quicklook images from CrIS/VIIRS Fusion product files.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2018-06-18.
Copyright (c) 2018 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os
from os.path import basename, dirname, curdir, abspath, isdir, isfile, exists, splitext, join as pjoin
import sys
import shutil
import logging
import traceback
import time
from glob import glob
from subprocess import call, check_call
from datetime import datetime

from utils import execution_time, create_dir, satellite_esdt, product_filename

# every module should have a LOG object
LOG = logging.getLogger(__file__)

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

def prepare_env(pkg_root):

    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = ':'.join([pjoin(pkg_root, 'env/lib')])
    env['PATH'] = ':'.join([pjoin(pkg_root, 'env/bin'),
                            '/usr/bin:/bin'])

    LOG.debug("env['PATH'] = \n\t{}".format(env['PATH'].replace(':','\n\t')))
    LOG.debug("env['LD_LIBRARY_PATH'] = \n\t{}".format(env['LD_LIBRARY_PATH'].replace(':','\n\t')))

    return env

def run_fusion_quicklooks(fusion_dir, geo_dir, **kwargs):

    bin_dir = kwargs['bin_dir']
    out_dir = kwargs['out_dir']
    fusion_ql_binary = kwargs['fusion_ql_binary']
    satellite = kwargs['satellite']
    env = kwargs['env']

    fusion_prefix = {'snpp':'VNP', 'noaa20':'VJ1'}[satellite]
    fusion_files = sorted(glob(pjoin(fusion_dir,'{}02FSN*.nc'.format(fusion_prefix))))
    geo_files = sorted(glob(pjoin(geo_dir,'{}03MOD*.nc'.format(fusion_prefix))))

    LOG.info('fusion_prefix: {}'.format(fusion_prefix))
    LOG.info('fusion_glob: {}'.format(pjoin(fusion_dir,'{}02FSN*.nc'.format(fusion_prefix))))
    LOG.info('fusion_files: {}'.format(fusion_files))
    LOG.info('geo_files: \n\t{}'.format(' \n\t'.join(geo_files)))
    LOG.info('fusion_files: \n\t{}'.format(' \n\t'.join(fusion_files)))

    dt = datetime.strptime(fusion_files[0].split('.')[1],'A%Y%j')
    year = dt.utctimetuple().tm_year
    jday = dt.utctimetuple().tm_yday

    # Create the output directory, and change to it.
    current_dir = os.getcwd()
    create_dir(out_dir)
    os.chdir(out_dir)

    rc_fusion_ql = 0

    #run matlab
    cmd = '{}/{} /opt/matlab/2015b/ {} {} {}/ {}/  >> fusion_quicklooks.log'.format(
        bin_dir,
        fusion_ql_binary,
        year,
        jday,
        fusion_dir,
        geo_dir
        )

    # Run the Matlab Fusion code
    LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))
    rc_fusion_ql = check_call([cmd], shell=True, env=env)
    LOG.debug('check_call() return value: {}'.format(rc_fusion_ql))

    # Move matlab file to the output directory
    fusion_ql_files = glob('*.png')
    if len(fusion_ql_files) != 0:
        LOG.info('Found Fusion quicklook files {}.'.format(
            ', '.join([basename(x) for x in fusion_ql_files])
            ))
    else:
        LOG.error('There are no Fusion quicklook files "*.png", aborting')
        rc_fusion_ql = 1
        return rc_fusion_ql

    os.chdir(current_dir)

    return rc_fusion_ql, fusion_ql_files

def main(args):

    pkg_root = abspath(args[0])
    fusion_dir = abspath(args[1])
    geo_dir = abspath(args[2])
    out_dir = abspath(args[3])
    satellite = args[4]

    LOG.info('Package root = {}'.format(pkg_root))
    LOG.info('Fusion dir = {}'.format(fusion_dir))
    LOG.info('Geo dir = {}'.format(geo_dir))
    LOG.info('Output dir = {}'.format(out_dir))
    LOG.info('satellite = {}'.format(satellite))

    py_env_dir = pjoin(pkg_root, 'env')
    bin_dir = pjoin(pkg_root, 'bin')
    #work_dir = os.path.abspath(os.path.curdir)

    LOG.debug('py_env_dir = {}'.format(py_env_dir))

    env = prepare_env(pkg_root)

    kwargs = {}
    kwargs['bin_dir'] = bin_dir
    kwargs['out_dir'] = out_dir
    kwargs['env'] = env

    # Creating the fusion quicklooks...


    kwargs['satellite'] = satellite
    kwargs['fusion_ql_binary'] = 'run_plot_globalVIIRSfusion_fct.sh'

    rc_fusion_ql, fusion_ql_files = run_fusion_quicklooks(fusion_dir, geo_dir, **kwargs)

    # Creating the matlab fusion file failed, exiting...
    if rc_fusion_ql != 0:
        return rc_fusion_ql

    LOG.debug('run_fusion_quicklooks() return value: {}'.format(rc_fusion_ql))
    LOG.debug('run_fusion_quicklooks() generated {}...'.format(fusion_ql_files))

    LOG.debug('python return value = {}'.format(rc_fusion_ql))
    return rc_fusion_ql



if __name__ == '__main__':
    args = sys.argv[1:]

    usage = \
        """
        Usage: 'run_fusion_1gran.sh PACKAGE_ROOT FUSION_DIR GEO_DIR OUTPUT_DIR SATELLITE'
        where:
            PKG_ROOT : Top level package dir (e.g.: 'dist/'
            FUSION_DIR : Directory containing V02FSN files"
            GEO_DIR : Directory containing V03MOD files matching the files in FUSION_DIR"
            OUTPUT_DIR : Directory where fusion quicklooks are generated"
            SATELLITE : The satellite providing the input data (snpp or noaa20)"
        """

    if len(args) != 5:
        print(usage)

    verbosity = 3
    setup_logging(verbosity)

    sys.exit(main(args))
