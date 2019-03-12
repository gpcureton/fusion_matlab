#!/usr/bin/env python
# encoding: utf-8
"""
Script to run fusion code (either for MODIS/AIRS or VIIRS/CrIS)
for one selected imager granule

Note, it is assumed that imager/sounder collocation files are available (one per sounder file)

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2017-08-10.
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

def run_fusion_matlab(geo_files, l1b_files, sounder_files, collo_files, **kwargs):

    bin_dir = kwargs['bin_dir']
    anc_paths = kwargs['anc_paths']
    fusion_binary = kwargs['fusion_binary']
    out_dir = kwargs['out_dir']
    matlab_file_glob = kwargs['matlab_file_glob']
    env = kwargs['env']

    rc_fusion = 0

    # run matlab
    cmd = '{}/{} /opt/matlab/2015b/ {} {} {} {} {}'.format(
        bin_dir,
        fusion_binary,
        ' '.join(geo_files),
        ' '.join(l1b_files),
        ' '.join(sounder_files),
        ' '.join(collo_files),
        ' '.join(anc_paths)
        )
    #cmd = 'touch fusion_output.mat' # Dummy
    #cmd = 'echo "No command"' # Dummy

    # Create the output directory
    current_dir = os.getcwd()
    create_dir(out_dir)

    # Run the Matlab Fusion code
    LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))

    rc_fusion = check_call([cmd], shell=True, env=env)
    LOG.debug('check_call() return value: {}'.format(rc_fusion))


    # Move matlab file to the output directory
    LOG.debug("glob('{}') = {}".format(matlab_file_glob, glob(matlab_file_glob)))
    matlab_file = glob(matlab_file_glob)
    if len(matlab_file) != 0:
        matlab_file = matlab_file[0]
        LOG.info('Found Matlab file "{}", moving to {}...'.format(matlab_file, out_dir))
        if exists(pjoin(out_dir, matlab_file)):
            LOG.info('{} exists, removing...'.format(pjoin(out_dir, matlab_file)))
            os.remove(pjoin(out_dir, matlab_file))
        shutil.move(matlab_file, out_dir)
        matlab_file = glob(pjoin(out_dir, matlab_file))[0]
    else:
        LOG.error('There are no Matlab files "{}" to convert, aborting'.format(matlab_file_glob))
        rc_fusion = 1
        return rc_fusion, None

    return rc_fusion, matlab_file

def convert_matlab_to_netcdf(matlab_file, l1b_file, **kwargs):

    bin_dir = kwargs['bin_dir']
    out_dir = kwargs['out_dir']
    satellite = kwargs['satellite']
    conversion_bin = kwargs['conversion_bin']
    env = kwargs['env']

    rc_fusion = 0

    dt = kwargs['input_dt']
    dt_create = datetime.utcnow()

    # Create the output directory
    current_dir = os.getcwd()

    LOG.debug("l1b_file = {}".format(l1b_file))

    # Copy the un-fused level-1b file to the work directory as a template...
    if satellite=='snpp' or satellite=='noaa20':
        esdt = satellite_esdt('V02FSN', satellite)
        collection = int(basename(l1b_file).split('.')[3])
        fused_l1b_file_new = product_filename(esdt, collection, dt, created=dt_create)
    if satellite=='aqua':
        esdt = satellite_esdt('M02FSN', satellite)
        collection = int(basename(l1b_file).split('.')[3])
        fused_l1b_file_new = product_filename(esdt, collection, dt, created=dt_create)

    unfused_l1b_file = pjoin(out_dir, 'unfused', os.path.basename(fused_l1b_file_new))
    unfused_l1b_dir = os.path.dirname(unfused_l1b_file)
    create_dir(unfused_l1b_dir)

    if os.path.exists(unfused_l1b_file):
        LOG.info('{} exists, removing...'.format(unfused_l1b_file))
        os.remove(unfused_l1b_file)

    LOG.info('Copying {} to {}'.format(l1b_file, unfused_l1b_file))
    shutil.copy(l1b_file, unfused_l1b_file)


    # Removing the fused  file if it exists
    fused_l1b_file = pjoin(out_dir, os.path.basename(unfused_l1b_file))
    if os.path.exists(fused_l1b_file):
        LOG.info('{} exists, removing...'.format(fused_l1b_file))
        os.remove(fused_l1b_file)

    # Convert the Matlab file to the desired format...
    #cmd = 'cp {} {}'.format(unfused_l1b_file, out_dir) # Dummy
    cmd = '{}  {} {} {}'.format(conversion_bin, unfused_l1b_file, matlab_file, out_dir)
    LOG.debug(cmd)
    rc_fusion = check_call([cmd], shell=True, env=env)

    # Determine success...
    fused_l1b_file = glob(fused_l1b_file)
    if len(fused_l1b_file) != 0:
        fused_l1b_file = fused_l1b_file[0]
    else:
        LOG.error('There is no fused file {}, aborting'.format(fused_l1b_file))
        rc_fusion = 1
        return rc_fusion, None

    # Remove the unfused dir...
    LOG.info('Removing the unfused level-1b dir {} ...'.format(unfused_l1b_dir))
    shutil.rmtree(unfused_l1b_dir)

    return rc_fusion, fused_l1b_file

def main(args):

    pkg_root = args[0]
    inst_pair = 'modis_airs' if int(args[1])==1 else 'cris_viirs'
    in_dir = args[2]
    out_dir = args[3]

    LOG.info('Package root = {}'.format(pkg_root))
    LOG.info('Instrument pair = {}'.format(inst_pair))
    LOG.info('Input dir = {}'.format(in_dir))
    LOG.info('Output dir = {}'.format(out_dir))

    py_env_dir = pjoin(pkg_root, 'env')
    bin_dir = pjoin(pkg_root, 'bin')
    anc_dir = pjoin(pkg_root, 'luts')

    LOG.debug('py_env_dir = {}'.format(py_env_dir))

    env = prepare_env(pkg_root)

    kwargs = {}
    kwargs['bin_dir'] = bin_dir
    kwargs['env'] = env
    kwargs['out_dir'] = out_dir

    # Creating the matlab fusion file...

    if inst_pair == 'modis_airs':

        imager = 'modis'
        sounder = 'airs'

        geo_files = sorted(glob(pjoin(in_dir, 'MYD03.*.hdf'))) #imager geolocation file
        l1b_files = sorted(glob(pjoin(in_dir, 'MYD021KM.*.hdf'))) #imager L1B file

        # sounder files (as many as needed to cover the imager file):
        sounder_files = sorted(glob(pjoin(in_dir, 'AIRS.*.hdf')))

        # imager/sounder collocation files (one per sounder file):
        collo_files = sorted(glob(pjoin(in_dir, 'colloc*airs*modis*.nc')))

        kwargs['input_dt'] = [
                datetime.strptime('.'.join(basename(geo_file).split('.')[1:3]), 'A%Y%j.%H%M') for geo_file in geo_files]

        # Various lut files and directories:
        kwargs['anc_paths'] = [pjoin(anc_dir, 'L2.chan_prop.2005.03.01.v9.5.1.txt'),
                               pjoin(anc_dir, 'modis_aqua.srf.nc'),
                               pjoin(anc_dir, 'modis_conv_error_2005.mat')]

        LOG.info('geo_files: \n\t{}'.format(' \n\t'.join(geo_files)))
        LOG.info('l1b_files: \n\t{}'.format(' \n\t'.join(l1b_files)))
        LOG.info('sounder_files: \n\t{}'.format(' \n\t'.join(sounder_files)))
        LOG.info('collo_files: \n\t{}'.format(' \n\t'.join(collo_files)))
        LOG.info('input_dt: \n\t{}'.format(' \n\t'.join([str(dt) for dt in kwargs['input_dt']])))

        kwargs['satellite'] = {'MOD':'terra', 'MYD':'aqua'}[basename(geo_files[0]).split('.')[0][:3]]
        kwargs['fusion_binary'] = 'run_imagersounderfusion_M.sh'
        kwargs['matlab_file_glob'] = 'fusion_output.mat'
        kwargs['conversion_bin'] = pjoin(py_env_dir, 'bin', 'l1b-fusion-modis-airs')

        # The l1b file for which we are making the fusion file...
        l1b_file = l1b_files[0]
        kwargs['input_dt'] = kwargs['input_dt'][0]

    elif inst_pair == 'cris_viirs':

        imager = 'viirs'
        sounder = 'cris'

        geo_files = sorted(glob(pjoin(in_dir, 'VNP03MOD.*.nc'))) #imager geolocation file
        l1b_files = sorted(glob(pjoin(in_dir, 'VNP02MOD.*.nc'))) #imager L1B file

        # sounder files (as many as needed to cover the imager file):
        sounder_files = sorted(glob(pjoin(in_dir, 'SNDR.*.CRIS.*.nc')))

        # imager/sounder collocation files (one per sounder file):
        collo_files = sorted(glob(pjoin(in_dir, 'colloc*cris*viirs*.nc')))

        #kwargs['input_dt'] = datetime.strptime('.'.join(basename(geo_file).split('.')[1:3]), 'A%Y%j.%H%M')
        kwargs['input_dt'] = [
                datetime.strptime('.'.join(basename(geo_file).split('.')[1:3]), 'A%Y%j.%H%M') for geo_file in geo_files]

        # Process the 1st granule in the original list, requirering +-1 context granules.
        gran_idx = 1
        geo_files = geo_files[gran_idx-1:gran_idx+2]
        l1b_files = l1b_files[gran_idx-1:gran_idx+2]
        sounder_files = sounder_files[gran_idx-1:gran_idx+2]
        collo_files = collo_files[gran_idx-1:gran_idx+2]
        kwargs['input_dt'] = kwargs['input_dt'][gran_idx-1:gran_idx+2]

        LOG.info('geo_files: \n\t{}'.format(' \n\t'.join(geo_files)))
        LOG.info('l1b_files: \n\t{}'.format(' \n\t'.join(l1b_files)))
        LOG.info('sounder_files: \n\t{}'.format(' \n\t'.join(sounder_files)))
        LOG.info('collo_files: \n\t{}'.format(' \n\t'.join(collo_files)))
        LOG.info('input_dt: \n\t{}'.format(' \n\t'.join([str(dt) for dt in kwargs['input_dt']])))

        # Various lut files and directories:
        kwargs['anc_paths'] = [pjoin(anc_dir, 'modis_aqua.srf.nc'),
                               pjoin(anc_dir, 'NG_VIIRS_NPP_RSR_filtered_Oct2011_BA/')]

        kwargs['satellite'] = {'VNP':'snpp', 'VJ1':'noaa20'}[basename(geo_files[0]).split('.')[0][:3]]
        kwargs['fusion_binary'] = 'run_imagersounderfusion_V.sh'
        kwargs['matlab_file_glob'] = 'fusion_output.mat'
        kwargs['conversion_bin'] = pjoin(py_env_dir, 'bin', 'l1b-fusion-viirs-cris')

        # The l1b file for which we are making the fusion file...
        gran_idx = 1
        l1b_file = l1b_files[gran_idx]
        kwargs['input_dt'] = kwargs['input_dt'][gran_idx]
        LOG.debug('l1b_file = {}'.format(l1b_file))
        LOG.debug('input_dt = {}'.format(kwargs['input_dt']))

    rc_fusion, matlab_file = run_fusion_matlab(geo_files, l1b_files, sounder_files, collo_files, **kwargs)

    # Creating the matlab fusion file failed, exiting...
    if rc_fusion != 0:
        return rc_fusion

    LOG.debug('run_fusion_matlab() return value: {}'.format(rc_fusion))
    LOG.debug('run_fusion_matlab() generated {}...'.format(matlab_file))

    # Now that we've computed the Matlab file, convert to a NetCDF file...
    rc_fusion, fused_l1b_file = convert_matlab_to_netcdf(matlab_file, l1b_file, **kwargs)

    # Converting the Matlab fusion file to NetCDF failed, exiting...
    if rc_fusion != 0:
        return rc_fusion

    LOG.debug('convert_matlab_to_netcdf() generated {}...'.format(fused_l1b_file))

    LOG.debug('python return value = {}'.format(rc_fusion))
    return rc_fusion



if __name__ == '__main__':
    args = sys.argv[1:]

    usage = \
        """
        Usage: 'run_fusion_1gran.sh PACKAGE_ROOT INSTRUMENT_PAIR INDIR OUTDIR'
        where:
            PKG_ROOT : Top level package dir (e.g.: 'dist/'
            INSTRUMENT_PAIR : 1 (for MODIS/AIRS) or 2 (for VIIRS/CrIS) "
        """

    if len(args) != 4:
        print(usage)

    verbosity = 3
    setup_logging(verbosity)

    sys.exit(main(args))
