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
from os.path import basename, dirname, curdir, abspath, isdir, isfile, exists, join as pjoin
import sys
import string
import shutil
import logging
import traceback
from glob import glob
from subprocess import check_output, check_call, CalledProcessError
from netCDF4 import Dataset
from pyhdf.SD import SD, SDC

from flo.computation import Computation
from flo.builder import WorkflowNotReady
from timeutil import TimeInterval, datetime, timedelta
from flo.util import symlink_inputs_to_working_dir

from glutil.software import delivered_software, support_software, runscript
from glutil.catalogs import dawg_catalog

from utils import create_dir

# every module should have a LOG object
LOG = logging.getLogger(__file__)

class FUSION_MATLAB(Computation):

    parameters = ['granule', 'satellite', 'version']

    outputs = ['fused_l1b']

    def find_contexts(self, time_interval, satellite, delivery_id):
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

    def _add_modis_l1b_geo_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=5)

        myd03_interval = TimeInterval(granule, granule+granule_length-timedelta(seconds=1))
        myd03 = dawg_catalog.files('aqua','MOD03', myd03_interval, version='c6')
        if myd03 == []:
            raise WorkflowNotReady('Unable to find matching MYD03 granule for interval {}'.format(myd03_interval))
        myd03_file = myd03[0]
        LOG.debug('MYD03 granule path: {}'.format(myd03_file.path))
        #if not exists(myd03_file.path):
            #raise WorkflowNotReady('Nominally valid MYD03 path {} for granule {} or interval {} does not exist, possible DB corruption'.format(myd03_file.path, granule, myd03_interval))

        task.input('geo', myd03_file)

    def _add_modis_l1b_m_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=5)

        myd021km_interval = TimeInterval(granule, granule+granule_length-timedelta(seconds=1))
        myd021km = dawg_catalog.files('aqua','MOD021KM', myd021km_interval, version='c6')
        if myd021km == []:
            raise WorkflowNotReady('Unable to find matching MYD021KM granule for interval {}'.format(myd021km_interval))
        myd021km_file = myd021km[0]
        LOG.debug('MYD021KM granule path: {}'.format(myd021km_file.path))
        #if not exists(myd021km_file.path):
            #raise WorkflowNotReady('Nominally valid MYD021KM path {} for granule {} or interval {} does not exist, possible DB corruption'.format(myd021km_file.path, granule, myd021km_interval))

        task.input('l1b', myd021km_file)

    def _add_airs_l1b_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=5)

        myd03_interval = TimeInterval(granule, granule+granule_length-timedelta(seconds=1))
        myd03 = dawg_catalog.files('aqua','MOD03', myd03_interval, version='c6')
        if myd03 == []:
            raise WorkflowNotReady('Unable to find matching MYD03 granule for interval {}'.format(myd03_interval))
        myd03_file = myd03[0]
        LOG.debug('MYD03 granule path: {}'.format(myd03_file.path))
        #if not exists(myd03_file.path):
            #raise WorkflowNotReady('Nominally valid MYD03 path {} for granule {} or interval {} does not exist, possible DB corruption'.format(myd03_file.path, granule, myd03_interval))

        buf = 0 # seconds
        airs_begin = myd03_file.begin_time - timedelta(seconds=buf)
        airs_end = myd03_file.end_time + timedelta(seconds=buf)
        airs_interval = TimeInterval(airs_begin, airs_end)
        airs = dawg_catalog.files('aqua', 'AIRIBRAD', airs_interval)
        if airs == []:
            raise WorkflowNotReady('Unable to find matching AIRS granules for interval {}'.format(airs_interval))

        for idx, airs_file in enumerate(airs):
            LOG.debug('AIRS granule {}: {} -> {}'.format(idx, airs_file.begin_time, airs_file.end_time))
            task.input('sounder_{}'.format(idx),  airs_file)

    def build_task_aqua(self, context, task):
        '''
        Build up a set of inputs for a single context
        '''

        LOG.debug("context = {}".format(context))

        self._add_modis_l1b_geo_input(context, task)
        self._add_modis_l1b_m_input(context, task)
        self._add_airs_l1b_input(context, task)

    def _add_viirs_l1b_geo_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        vgeom = dawg_catalog.file(satellite, 'VGEOM', granule, version='2.0.2')
        task.input('geo', vgeom)

    def _add_viirs_l1b_m_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        vl1b = dawg_catalog.file(satellite, 'VL1BM', granule, version='2.0.2')
        task.input('l1b', vl1b)

    def _add_cris_l1b_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=6)
        cris_interval = TimeInterval(granule, granule+granule_length-timedelta(seconds=1))
        cris = dawg_catalog.files(satellite, 'CL1B', cris_interval, version='v1.0rc8')
        if cris == []:
            raise WorkflowNotReady('Unable to find matching CrIS granules for interval {}'.format(cris_interval))

        for idx, cris_file in enumerate(cris):
            LOG.debug('CrIS granule {}: {} -> {}'.format(idx, cris_file.begin_time, cris_file.end_time))
            task.input('sounder_{}'.format(idx),  cris_file)

    def build_task_snpp(self, context, task):
        '''
        Build up a set of inputs for a single context
        '''

        LOG.debug("context = {}".format(context))

        self._add_viirs_l1b_geo_input(context, task)
        self._add_viirs_l1b_m_input(context, task)
        self._add_cris_l1b_input(context, task)

    def build_task(self, context, task):

        satellite = context['satellite']

        if satellite == 'snpp':
            self.build_task_snpp(context, task)
        elif satellite == 'aqua':
            self.build_task_aqua(context, task)
        else:
            raise ValueError('Invalid satellite: {}'.format(satellite))

    def mend_viirs_l1b(self, geo, l1b):

        out_fn = basename(l1b).replace('.nc', '.bowtie_restored.nc')
        LOG.info('out_fn = {}'.format(out_fn))
        shutil.copy(l1b, out_fn)

        viirsmend = support_software.lookup('viirsmend', version='1.2.12')
        py_exe = pjoin(viirsmend.path,'bin/python')
        viirsmend_exe = pjoin(viirsmend.path,'bin/viirsl1mend')
        cmd = '{} {} {} {}'.format(py_exe, viirsmend_exe, out_fn, geo)

        LOG.info('cmd = {}'.format(cmd))
        runscript(cmd, requirements=[viirsmend])

        return out_fn

    def cris_viirs_collocation(self, inputs):

        LOG.info('inputs = {}'.format(inputs))
        input_keys = inputs.keys()
        input_keys.sort()

        sounder_keys = [key for key in input_keys if 'sounder' in key]
        cris_files = [inputs[key] for key in sounder_keys]

        vgeom_file = inputs['geo']

        crisviirs = support_software.lookup('collopak', version='0.1.65')
        crisviirs_exe = pjoin(crisviirs.path,'bin/crisviirs')
        for cris_file in cris_files:
            cmd = '{} {} {} > /dev/null'.format(crisviirs_exe, cris_file, vgeom_file)

            LOG.info('cmd = {}'.format(cmd))
            runscript(cmd, requirements=[])

        return glob('colloc.*.nc')

    def airs_modis_collocation(self, inputs):

        LOG.info('inputs = {}'.format(inputs))
        input_keys = inputs.keys()
        input_keys.sort()

        sounder_keys = [key for key in input_keys if 'sounder' in key]
        airs_files = [inputs[key] for key in sounder_keys]

        modis_file = inputs['geo']

        airsmodis = support_software.lookup('collopak', version='0.1.65')
        airsmodis_exe = pjoin(airsmodis.path,'bin/airsmod')
        for airs_file in airs_files:
            cmd = '{} {} {} > /dev/null'.format(airsmodis_exe, airs_file, modis_file)

            LOG.info('cmd = {}'.format(cmd))
            runscript(cmd, requirements=[])

        return glob('colloc.*.nc')

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
            support_software.lookup('matlab', '2015b').path,
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
            rc_fusion = 0
            runscript(cmd, requirements=[], env=env)
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
            matlab_file = glob(pjoin(out_dir, matlab_file))[0]
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
        satellite = kwargs['satellite']
        granule = kwargs['granule']

        rc_fusion = 0

        dt = datetime.strptime(basename(matlab_file), matlab_file_dt_filespec)
        dt_string = dt.strftime(matlab_file_dt_str)

        LOG.debug('dt_string = {}'.format(dt_string))
        LOG.debug('granule = {}'.format(granule))

        # Create the output directory
        current_dir = os.getcwd()

        # Copy the un-fused level-1b file to the work directory as a template...
        unfused_l1b_file = pjoin(out_dir, 'unfused', basename(l1b_file))
        unfused_l1b_dir = dirname(unfused_l1b_file)
        create_dir(unfused_l1b_dir)

        if exists(unfused_l1b_file):
            LOG.debug('{} exists, removing...'.format(unfused_l1b_file))
            os.remove(unfused_l1b_file)

        LOG.debug('Copying {} to {}'.format(l1b_file, unfused_l1b_file))
        shutil.copy(l1b_file, unfused_l1b_file)

        # Removing the fused file if it exists
        fused_l1b_file = pjoin(out_dir, basename(unfused_l1b_file))
        if exists(fused_l1b_file):
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
            rc_fusion = 0
            runscript(cmd, requirements=[], env=env)
        except CalledProcessError as err:
            rc_fusion = err.returncode
            LOG.error("CF converter {} returned a value of {}".format(conversion_bin, rc_fusion))
            return rc_fusion, []

        # Determine success...
        fused_l1b_file = glob(fused_l1b_file)
        if len(fused_l1b_file) != 0:
            fused_l1b_file = fused_l1b_file[0]
            LOG.debug('Found final fused output file "{}"'.format(fused_l1b_file))
        else:
            LOG.error('There is no fused file {}, aborting'.format(fused_l1b_file))
            rc_fusion = 1
            return rc_fusion, []

        # Remove the unfused dir...
        LOG.debug('Removing the unfused level-1b dir {} ...'.format(unfused_l1b_dir))
        shutil.rmtree(unfused_l1b_dir)

        # Move the final fused file to the work directory
        if satellite=='snpp':
            if 'VL1BM' in l1b_file:
                fused_l1b_file_new = dt.strftime('VNP02FSN.A%Y%j.%H%M.000.CTIME.nc')
            elif 'VNP02MOD' in l1b_file:
                fused_l1b_file_new = 'VNP02FSN.{}.CTIME.hdf'.format('.'.join(l1b_file.split('.')[1:4]))
            else:
                pass
        if satellite=='aqua':
            fused_l1b_file_new = 'MYD02FSN.{}.CTIME.hdf'.format('.'.join(l1b_file.split('.')[1:4]))

        dt_create = datetime.utcnow()
        fused_l1b_file_new = fused_l1b_file_new.replace('CTIME', dt_create.strftime('%Y%j%H%M%S'))

        LOG.debug('Moving "{}" to "{}" ...'.format(fused_l1b_file, fused_l1b_file_new))
        shutil.move(fused_l1b_file, pjoin(current_dir, fused_l1b_file_new))

        fused_l1b_file = glob(pjoin(current_dir, fused_l1b_file_new))[0]

        # Remove the fused_outputs directory
        LOG.debug('Removing the fused_outputs dir {} ...'.format(out_dir))
        shutil.rmtree(out_dir)

        return rc_fusion, fused_l1b_file

    def update_global_attrs(self, netcdf_file, readme_file, **kwargs):

        satellite = kwargs['satellite']

        # Get the git repo information
        repo_attrs = []
        try:
            LOG.debug('Opening {}...'.format(readme_file))
            readme_obj = open(readme_file, 'ro')
            line_obj = readme_obj.readlines()
            for idx, line in enumerate(line_obj):
                if '.git' in line:
                    repo_line = line.lstrip(' -*,*').rstrip(' -*,;')
                    commit_line = line_obj[idx+1].lstrip(' -*,;').rstrip(' -*,;')
                    git_line = '{}; {}'.format(repo_line, commit_line).replace('\n', '')
                    LOG.debug('{}'.format(git_line))
                    repo_attrs.append(git_line)
        except Exception:
            LOG.debug(traceback.format_exc())

        readme_obj.close()

        # Update the various file global attributes
        LOG.debug('Adding attributes to {} ...'.format(netcdf_file))
        if netcdf_file.split('.')[-1] == 'nc':
            args = (netcdf_file, "a")
            kwargs = {'format': "NETCDF4"}
            file_open = Dataset
        elif netcdf_file.split('.')[-1] == 'hdf':
            args = (netcdf_file, SDC.WRITE)
            kwargs = {}
            file_open = SD

        try:

            file_obj = file_open(*args, **kwargs)

            # Update the attributes, moving to the end
            for idx, attr in enumerate(repo_attrs):
                LOG.debug('{}'.format(attr))
                setattr(file_obj, 'source_git_repo_{}'.format(idx), attr)

            setattr(file_obj, 'SIPS_version', basename(dirname(readme_file)))

        except Exception:
            LOG.warning("\tProblem setting attributes in output file {}".format(netcdf_file))
            LOG.debug(traceback.format_exc())

        if netcdf_file.split('.')[-1] == 'nc':
            file_obj.close()
        elif netcdf_file.split('.')[-1] == 'hdf':
            file_obj.end()

        return

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
        work_dir = abspath(curdir)
        LOG.debug("working dir = {}".format(work_dir))

        # What are out inputs?
        for input in inputs.keys():
            inputs_dir = dirname(inputs[input])
            LOG.debug("inputs['{}'] = {}".format(input,inputs[input]))
        LOG.debug("Inputs dir = {}".format(inputs_dir))

        # Run viirsmend on the viirs level-1b, and generate the CrIS/VIIRS collocation
        if satellite == 'snpp':
            geo = inputs['geo']
            l1b = self.mend_viirs_l1b(inputs['geo'], inputs['l1b'])
            sounder_keys = [key for key in inputs.keys() if 'sounder' in key]
            sounder = [inputs[key] for key in sounder_keys]
            collo = self.cris_viirs_collocation(inputs)

        # Generate the AIRS/MODIS collocation
        if satellite == 'aqua':
            geo = inputs['geo']
            l1b = inputs['l1b']
            sounder_keys = [key for key in inputs.keys() if 'sounder' in key]
            sounder = [inputs[key] for key in sounder_keys]
            collo = self.airs_modis_collocation(inputs)

        LOG.info('geo = {}'.format(geo))
        LOG.info('l1b = {}'.format(l1b))
        LOG.info('sounder = {}'.format(sounder))
        LOG.info('collo = {}'.format(collo))

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
        kwargs['satellite'] = satellite
        kwargs['granule'] = granule

        if satellite=='snpp':
            kwargs['fusion_binary'] = 'run_imagersounderfusion_V.sh'
            kwargs['matlab_file_glob'] = 'fusion_viirs_*_t*.mat'
            kwargs['matlab_file_dt_str'] = '%Y%m%d_t%H%M%S'
            kwargs['matlab_file_dt_filespec'] = 'fusion_viirs_{}.mat'.format(kwargs['matlab_file_dt_str'])
            kwargs['conversion_bin'] = pjoin(envroot, 'bin', 'l1b-fusion-viirs-cris')
        elif satellite=='aqua':
            kwargs['fusion_binary'] = 'run_imagersounderfusion_M.sh'
            kwargs['matlab_file_glob'] = 'fusion_modis_*.*.mat'
            kwargs['matlab_file_dt_str'] = '%Y%j.%H%M'
            kwargs['matlab_file_dt_filespec'] = 'fusion_modis_{}.mat'.format(kwargs['matlab_file_dt_str'])
            kwargs['conversion_bin'] = pjoin(envroot, 'bin', 'l1b-fusion-modis-airs')
        else:
            return {}

        # Run the fusion_matlab package
        rc_fusion, matlab_file = self.run_fusion_matlab(
                                                   geo,
                                                   l1b,
                                                   sounder,
                                                   collo,
                                                   **kwargs
                                                  )

        LOG.debug('run_fusion_matlab() return value: {}'.format(rc_fusion))
        LOG.info('run_fusion_matlab() generated {}'.format(matlab_file))

        # Now that we've computed the Matlab file, convert to a NetCDF file...
        rc_fusion, fused_l1b_file = self.convert_matlab_to_netcdf(matlab_file,
                                                                  l1b,
                                                                  **kwargs)

        LOG.debug('convert_matlab_to_netcdf() return value: {}'.format(rc_fusion))
        LOG.info('convert_matlab_to_netcdf() generated {}'.format(fused_l1b_file))

        # Update some global attributes in the output file
        readme_file =  pjoin(delivery.path, 'README.txt')
        self.update_global_attrs(basename(fused_l1b_file), readme_file, **kwargs)

        # The staging routine assumes that the output file is located in the work directory
        # "tmp******", and that the output path is to be prepended, so return the basename.

        out_fn = basename(fused_l1b_file)

        return {'fused_l1b': out_fn}
