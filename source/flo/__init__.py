#!/usr/bin/env python
# encoding: utf-8
"""

Purpose: Run the fusion_matlab package found at:

    https://gitlab.ssec.wisc.edu/eweisz/modis-airs-fusion
    https://gitlab.ssec.wisc.edu/eweisz/viirs-cris-fusion

and the CF-compliance converter found at:

    https://gitlab.ssec.wisc.edu/sounder-imager-fusion/fusion_glue/tree/develop/py/fusion_l1b

For Aqua AIRS/MODIS fusion, example inputs are...

    * NASA MODISS L1B geolocation and radiometric files:
        MYD03.A2015107.1435.006.2015108161923.hdf
        MYD021KM.A2015107.1435.006.2015108162652.hdf

    * NASA AIRS L1B files:
        AIRS.2015.04.17.145.L1B.AIRS_Rad.v5.0.23.0.G15173170714.hdf
        AIRS.2015.04.17.146.L1B.AIRS_Rad.v5.0.23.0.G15173170305.hdf

    * AIRS/MODIS collocation files
        colloc.airs_20150417T1429.modis_20150417T1435.nc
        colloc.airs_20150417T1435.modis_20150417T1435.nc

For Suomi-NPP CrIS/VIIRS fusion, example inputs are...

    * NASA VIIRS L1B geolocation and radiometric files:
        VGEOM_snpp_d20150417_t143600_c20170401181622.nc
        VL1BM_snpp_d20150417_t143600_c20170401181415.nc

    * NASA CrIS L1B files
        SNDR.SNPP.CRIS.20150417T1436.m06.g147.L1B_NSR.std.v01_00_00.T.161217000158.nc

    * CrIS/VIIRS collocation files
        colloc.cris_20150417T1436.viirs_20150417T1436.nc



Copyright (c) 2017 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os
from os.path import basename, join as pjoin
import sys
import string
import shutil
import logging
import traceback
from glob import glob
from datetime import datetime, timedelta
from subprocess import check_output

from flo.computation import Computation
from subprocess import check_call, CalledProcessError
from flo.time import TimeInterval, check_coverage, round_datetime
from flo.util import symlink_inputs_to_working_dir

from flo.dawg import DawgCatalog

from flo.sw.viirs_l1 import ViirsL1b2, ViirsGeo2
from flo.sw.cris_l1b import L1B
from flo.sw.collocation import CrisViirsCollocation

from flo.sw.lib.glutil import delivered_software, support_software
from orbnav_client import Client

from utils import create_dir

dawg_catalog = DawgCatalog('DAWG')

# every module should have a LOG object
LOG = logging.getLogger(__file__)

class FUSION_MATLAB(Computation):

    parameters = ['granule', 'satellite', 'version']

    outputs = ['fused_l1b']

    def find_contexts_gran_length(self, time_interval, satellite, delivery_id):
        '''
        Here we assume that the granule boundaries fall along 6-minute (snpp) or 5-minute (aqua)
        increments, starting at the top of the hour:

        SNPP: [0.,   6.,  12.,  18.,  24.,  30.,  36.,  42.,  48.,  54.]
        AQUA: [0.,   5.,  10.,  15.,  20.,  25.,  30.,  35.,  40.,  45.,  50., 55.]
        '''

        if satellite=='snpp':
            granule_length = timedelta(minutes=6)
        elif satellite=='aqua':
            granule_length = timedelta(minutes=5)
        else:
            return []

        return [{'satellite': satellite, 'version': delivery_id, 'granule': g}
                    for g in time_interval.contained_series(granule_length)]

    def find_contexts_orbnav(self, time_interval, satellite, delivery_id):

        '''
        Return the desired contexts falling within the given time interval.

        time_interval: TimeInterval(datetime_1, datetime_2)
        granule_length: timedelta(seconds=360)
        '''

        LOG.info("Running find_contexts()...")

        #
        # Here we assume that the granule boundaries fall along 6-minute
        # increments, starting at the top of the hour, so that there would be
        # 10 granules per hour at minutes:
        # [0.,   6.,  12.,  18.,  24.,  30.,  36.,  42.,  48.,  54.]
        #

        if satellite=='snpp':
            granule_length = timedelta(minutes=6)
        elif satellite=='aqua':
            granule_length = timedelta(minutes=5)
        else:
            return []

        # Generating List of all possible granules
        allGranuals = {
            g.left for g in time_interval.overlapping_interval_series(granule_length)}

        # Not all contexts are valid for and aerosol product: we only want the
        # day granules. This corresponds to a solar zenith angle < 85 degrees,
        # in this case. We will use OrbNav to determine whether a context
        # granule has any pixels in daylight.

        # Setting up OrbNav Client
        orbnavSatName = {'aqua': 'AQUA', 'snpp': 'SUOMI NPP'}
        orbnavClient = Client()

        try:
            # Get the datetimes when the satellite enters and/or leaves a
            # circle where the solar zenith angle is LESS than 85 degrees.
            orbnavOutput = orbnavClient.suncirc(sat=orbnavSatName[satellite],
                                                start=time_interval.left,
                                                end=time_interval.right,
                                                zenangl='85'
                                                )['data']

            # For each datetime when the satellite crosses the "terminator"
            # defined by solZenAngle=85 degrees, compute a 6-minute time
            # interval centered on the crossing. Any 6-minute granule which
            # has a start time falling in this interval will have daytime
            # pixels.
            half_gran_delta = timedelta(seconds=0.5 * granule_length.seconds)
            dayTimeIntervals = [
                TimeInterval(
                    terminator_intersect_time[0][0] - half_gran_delta,
                    terminator_intersect_time[1][0] + half_gran_delta,
                    right_open=True)
                for terminator_intersect_time in orbnavOutput
            ]

            # Generating list of OrbNav selected granules
            orbnavGranules = set()
            for interval in dayTimeIntervals:
                orbnavGranules = orbnavGranules | {
                    g.left for g in interval.overlapping_interval_series(granule_length)}

            # Finding Overlapping Granuals
            overlappingGranules = sorted(orbnavGranules & allGranuals)
        except:
            overlappingGranules = allGranuals

        return [{'satellite': satellite, 'version': delivery_id, 'granule': g}
                    for g in overlappingGranules]

    def find_contexts(self, time_interval, satellite, delivery_id):
        return self.find_contexts_gran_length(time_interval, satellite, delivery_id)

    def build_task_aqua(self, context, task):
        '''
        Build up a set of inputs for a single context
        '''

        # context aliases
        version = context['version']
        satellite = context['satellite']
        granule = context['granule']

        #
        # Add tasks for the download of the various inputs...
        #
        LOG.debug("context = {}".format(context))

        # Get the level-1b file references from DAWG
        myd03 = dawg_catalog.file('aqua','MOD03', granule, version='c6')
        myd021km = dawg_catalog.file('aqua','MOD021KM', granule, version='c6')
        #airs = dawg_catalog.file(satellite, 'AIRS', granule, version='c6')

        # Generate the CrIS/VIIRS collocation using the CrIS and VIIRS files from DAWG
        #collo_context = {
               #'satellite' : satellite,
               #'cris_granule': granule,
               #'viirs_granule': granule,
               #'cris_sdr_type': 'l1-v1.0rc8',
               #'viirs_sdr_type': 'l1-2.0.2',
               #'version': '0.1.59'
        #}
        #collo = CrisViirsCollocation().dataset('out').product(collo_context)

        # Download the required files
        task.input('geo', myd03)
        task.input('l1b', myd021km)
        #task.input('sounder',  airs)
        #task.input('collo', collo)

    def build_task_snpp(self, context, task):
        '''
        Build up a set of inputs for a single context
        '''


        # context aliases
        version = context['version']
        satellite = context['satellite']
        granule = context['granule']

        #
        # Add tasks for the download of the various inputs...
        #
        LOG.debug("context = {}".format(context))

        # Get the level-1b file references from DAWG
        vgeom = dawg_catalog.file(satellite, 'VGEOM', granule, version='2.0.2')
        vl1b = dawg_catalog.file(satellite, 'VL1BM', granule, version='2.0.2')
        cris = dawg_catalog.file(satellite, 'CL1B', granule, version='v1.0rc8')

        # Generate the CrIS/VIIRS collocation using the CrIS and VIIRS files from DAWG
        collo_context = {
               'satellite' : satellite,
               'cris_granule': granule,
               'viirs_granule': granule,
               'cris_sdr_type': 'l1-v1.0rc8',
               'viirs_sdr_type': 'l1-2.0.2',
               'version': '0.1.59'
        }
        collo = CrisViirsCollocation().dataset('out').product(collo_context)

        # Download the required files, and add them to the task inputs.
        task.input('geo', vgeom)
        task.input('l1b', vl1b)
        task.input('sounder',  cris)
        task.input('collo', collo)

    def build_task(self, context, task):

        satellite = context['satellite']

        if satellite == 'snpp':
            self.build_task_snpp(context, task)
        elif satellite == 'aqua':
            raise NotImplementedError("Fusion_matlab support for aqua hasn't been implemented yet.")
            self.build_task_aqua(context, task)
        else:
            raise ValueError('Invalid satellite: {}'.format(satellite))

    def run_fusion_matlab(self, geo_file, l1b_file, sounder_files, collo_files, **kwargs):
        '''
        Run the Matlab fusion binary on the input level-1b files to generate the *.mat output file.
        '''

        bin_dir = kwargs['bin_dir']
        anc_dir = kwargs['anc_dir']
        fusion_binary = kwargs['fusion_binary']
        out_dir = kwargs['out_dir']
        matlab_file_glob = kwargs['matlab_file_glob']
        env = kwargs['env']

        rc_fusion = 0

        #run matlab
        cmd = '{}/{} {} {} {} {} {} {}'.format(
            bin_dir,
            fusion_binary,
            support_software.lookup('matlab/2015b').path,
            geo_file,
            l1b_file,
            ' '.join(sounder_files),
            ' '.join(collo_files),
            anc_dir
            )

        # Create the output directory
        current_dir = os.getcwd()
        LOG.debug('The current directory is {}'.format(current_dir))
        create_dir(out_dir)

        # Run the Matlab Fusion code
        try:
            LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))
            rc_fusion = check_call([cmd], shell=True, env=env)
            #<<< Dummy <<<<<<<<<
            #shutil.copy('/mnt/sdata/geoffc/fusion_matlab/test_data/outputs/fusion_viirs_20150417_t143600.mat',
                        #pjoin(current_dir, 'fusion_viirs_20150417_t143600.mat'))
            #rc_fusion = 0
            #>>> Dummy >>>>>>>>>
        except CalledProcessError as err:
            rc_fusion = err.returncode
            LOG.error("Matlab binary {} returned a value of {}".format(fusion_binary, rc_fusion))
            return rc_fusion, []

        # Move matlab file to the output directory
        matlab_file = glob(matlab_file_glob)
        if len(matlab_file) != 0:
            matlab_file = matlab_file[0]
            LOG.debug('Found Matlab file "{}", moving to {}...'.format(matlab_file, out_dir))
            shutil.move(matlab_file, out_dir)
            matlab_file = glob(os.path.join(out_dir, matlab_file))[0]
        else:
            LOG.error('There are no Matlab files "{}" to convert, aborting'.format(matlab_file_glob))
            rc_fusion = 1
            return rc_fusion, []

        return rc_fusion, matlab_file


    def convert_matlab_to_netcdf(self, matlab_file, l1b_file, **kwargs):
        '''
        Transcode the Matlab *.mat file into a CF-compliant HDF4 (aqua) or NetCDF4 (snpp) file.
        '''

        py_interp = kwargs['py_interp']
        bin_dir = kwargs['bin_dir']
        anc_dir = kwargs['anc_dir']
        out_dir = kwargs['out_dir']
        matlab_file_dt_str = kwargs['matlab_file_dt_str']
        matlab_file_dt_filespec = kwargs['matlab_file_dt_filespec']
        conversion_bin = kwargs['conversion_bin']
        env = kwargs['env']

        rc_fusion = 0

        dt = datetime.strptime(os.path.basename(matlab_file), matlab_file_dt_filespec)
        dt_string = dt.strftime(matlab_file_dt_str)

        LOG.debug('dt_string = {}'.format(dt_string))

        # Create the output directory
        current_dir = os.getcwd()

        # Copy the un-fused level-1b file to the work directory as a template...
        unfused_l1b_file = os.path.join(out_dir, 'unfused', os.path.basename(l1b_file))
        unfused_l1b_dir = os.path.dirname(unfused_l1b_file)
        create_dir(unfused_l1b_dir)

        if os.path.exists(unfused_l1b_file):
            LOG.debug('{} exists, removing...'.format(unfused_l1b_file))
            os.remove(unfused_l1b_file)

        LOG.debug('Copying {} to {}'.format(l1b_file, unfused_l1b_file))
        shutil.copy(l1b_file, unfused_l1b_file)

        # Removing the fused file if it exists
        fused_l1b_file = os.path.join(out_dir, os.path.basename(unfused_l1b_file))
        if os.path.exists(fused_l1b_file):
            LOG.debug('{} exists, removing...'.format(fused_l1b_file))
            os.remove(fused_l1b_file)

        cmd = '{} {}  {} {} {}'.format(
                py_interp,
                conversion_bin,
                unfused_l1b_file,
                matlab_file,
                out_dir
                )

        # Convert the Matlab file to the desired format...
        try:
            LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))
            rc_fusion = check_call([cmd], shell=True, env=env)
            #<<< Dummy <<<<<<<<<
            #shutil.copy('/mnt/sdata/geoffc/fusion_matlab/test_data/outputs/VL1BM_snpp_d20150417_t143600_c20170401181415.nc',
                        #pjoin(out_dir, 'VL1BM_snpp_d20150417_t143600_c20170401181415.nc'))
            #rc_fusion = 0
            #>>> Dummy >>>>>>>>>
        except CalledProcessError as err:
            rc_fusion = err.returncode
            LOG.error("CF converter {} returned a value of {}".format(conversion_bin, rc_fusion))
            return rc_fusion, []

        # Determine success...
        fused_l1b_file = glob(fused_l1b_file)
        if len(fused_l1b_file) != 0:
            fused_l1b_file = fused_l1b_file[0]
        else:
            LOG.error('There is no fused file {}, aborting'.format(fused_l1b_file))
            rc_fusion = 1
            return rc_fusion, []

        # Remove the unfused dir...
        LOG.debug('Removing the unfused level-1b dir {} ...'.format(unfused_l1b_dir))
        shutil.rmtree(unfused_l1b_dir)

        # Move the final fused file to the work directory
        LOG.debug('Found final fused output file "{}", moving to {}...'.format(fused_l1b_file, current_dir))
        shutil.move(fused_l1b_file, current_dir)
        fused_l1b_file = glob(pjoin(current_dir, basename(fused_l1b_file)))[0]

        # Remove the fused_outputs directory
        LOG.debug('Removing the fused_outputs dir {} ...'.format(out_dir))
        shutil.rmtree(out_dir)

        return rc_fusion, fused_l1b_file

    def prepare_env(self, dist_root, inputs, context):
        LOG.debug("Running prepare_env()...")

        LOG.debug("package_root = {}".format(self.package_root))
        LOG.debug("dist_root = {}".format(dist_root))

        env = dict(os.environ)
        envroot = pjoin(dist_root, 'env')

        LOG.debug("envroot = {}".format(envroot))

        env['LD_LIBRARY_PATH'] = ':'.join([pjoin(envroot, 'lib'),
                                           pjoin(dist_root, 'lib')])
        env['PATH'] = ':'.join([pjoin(envroot, 'bin'),
                                pjoin(dist_root, 'bin'),
                                '/usr/bin:/bin'])


        LOG.debug("env['PATH'] = \n\t{}".format(env['PATH'].replace(':','\n\t')))
        LOG.debug("env['LD_LIBRARY_PATH'] = \n\t{}".format(env['LD_LIBRARY_PATH'].replace(':','\n\t')))

        return env

    def run_task(self, inputs, context):

        LOG.debug("Running run_task()...")

        for key in context.keys():
            LOG.debug("run_task() context['{}'] = {}".format(key, context[key]))

        granule = context['granule']
        satellite = context['satellite']
        delivery_id = context['version']

        # Get the location of the binary package
        delivery = delivered_software.lookup('fusion_matlab', delivery_id=delivery_id)
        dist_root = pjoin(delivery.path, 'dist')
        envroot = pjoin(dist_root, 'env')

        # Get the required  environment variables
        env = self.prepare_env(dist_root, inputs, context)

        # What is the path of the python interpreter
        py_interp = "{}/bin/python".format(envroot)
        LOG.debug("py_interp = '{}'".format(py_interp))

        # Where are we running the package
        work_dir = os.path.abspath(os.path.curdir)
        LOG.debug("working dir = {}".format(work_dir))

        # What are out inputs?
        for input in inputs.keys():
            inputs_dir = os.path.dirname(inputs[input])
            LOG.debug("inputs['{}'] = {}".format(input,inputs[input]))
        LOG.debug("Inputs dir = {}".format(inputs_dir))

        #py_env_dir = pjoin(pkg_root, 'env')
        bin_dir = pjoin(dist_root, 'bin')
        anc_dir = pjoin(dist_root, 'luts')
        out_dir = pjoin(work_dir, 'fused_outputs')

        # Setup the require keyword arguments for the fusion_matlab package
        kwargs = {}
        kwargs['py_interp'] = py_interp
        kwargs['bin_dir'] = bin_dir
        kwargs['anc_dir'] = anc_dir
        kwargs['env'] = env
        kwargs['out_dir'] = out_dir

        if satellite=='snpp':
            kwargs['fusion_binary'] = 'run_imagersounderfusion_V.sh'
            kwargs['matlab_file_glob'] = 'fusion_viirs_*_t*.mat'
            kwargs['matlab_file_dt_str'] = '%Y%m%d_t%H%M%S'
            kwargs['matlab_file_dt_filespec'] = 'fusion_viirs_{}.mat'.format(kwargs['matlab_file_dt_str'])
            kwargs['conversion_bin'] = os.path.join(envroot, 'bin', 'l1b-fusion-viirs-cris')
        elif satellite=='aqua':
            kwargs['fusion_binary'] = 'run_imagersounderfusion_M.sh'
            kwargs['matlab_file_glob'] = 'fusion_modis_*.*.mat'
            kwargs['matlab_file_dt_str'] = '%Y%j.%H%M'
            kwargs['matlab_file_dt_filespec'] = 'fusion_modis_{}.mat'.format(kwargs['matlab_file_dt_str'])
            kwargs['conversion_bin'] = os.path.join(envroot, 'bin', 'l1b-fusion-modis-airs')
        else:
            return {}

        # Run the fusion_matlab package
        rc_fusion, matlab_file = self.run_fusion_matlab(
                                                   inputs['geo'],
                                                   inputs['l1b'],
                                                   [inputs['sounder']],
                                                   [inputs['collo']],
                                                   **kwargs
                                                  )

        LOG.debug('run_fusion_matlab() return value: {}'.format(rc_fusion))
        LOG.info('run_fusion_matlab() generated {}'.format(matlab_file))

        # Now that we've computed the Matlab file, convert to a NetCDF file...
        rc_fusion, fused_l1b_file = self.convert_matlab_to_netcdf(matlab_file,
                                                                  inputs['l1b'],
                                                                  **kwargs)

        LOG.debug('convert_matlab_to_netcdf() return value: {}'.format(rc_fusion))
        LOG.info('convert_matlab_to_netcdf() generated {}'.format(fused_l1b_file))

        # The staging routine assumes that the output file is located in the work directory
        # "tmp******", and that the output path is to be prepended, so return the basename.

        out_fn = os.path.basename(fused_l1b_file)

        return {'fused_l1b': out_fn}
