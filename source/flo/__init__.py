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
        VNP03MOD.A2018033.1836.001.2018033235737.uwssec.nc
        VNP02MOD.A2018033.1836.001.2018033235822.uwssec.nc

    * NASA CrIS L1B files
        SNDR.SNPP.CRIS.20180202T1836.m06.g187.L1B.std.v2_0_15.T.180203004403.nc

    * CrIS/VIIRS collocation files
        colloc.cris_snpp.viirs_m_snpp.20180202T183600_183600.nc



Copyright (c) 2017 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os
from os.path import basename, dirname, curdir, abspath, isdir, isfile, exists, splitext, join as pjoin
import sys
import string
import shutil
import logging
import traceback
from glob import glob
from itertools import chain
from subprocess import check_output, check_call, CalledProcessError
from netCDF4 import Dataset
from pyhdf.SD import SD, SDC

from flo.computation import Computation
from flo.builder import WorkflowNotReady
from timeutil import TimeInterval, datetime, timedelta
from flo.util import symlink_inputs_to_working_dir

import sipsprod
from glutil import (
    #check_call,
    dawg_catalog,
    delivered_software,
    support_software,
    runscript,
    get_viirs_l1_luts,
    #prepare_env,
    #nc_gen,
    nc_compress,
    hdf_compress,
    reraise_as,
    set_official_product_metadata,
    FileNotFound
)

from utils import create_dir
from detect_bad_fusion import SFX, detect

# every module should have a LOG object
LOG = logging.getLogger(__name__)

class FUSION_MATLAB(Computation):

    parameters = ['granule', 'satellite', 'version']

    outputs = ['fused_l1b']

    def find_contexts(self, time_interval, satellite, version):
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

        return [{'satellite': satellite, 'version': version, 'granule': g}
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
        LOG.debug("Ingesting inputs for M02FSN version {} ...".format(context['version']))

        # Get the product definition for 'V02FSN'
        product = sipsprod.lookup_product_recurse('V02FSN', version=context['version'])

        LOG.debug("context = {}".format(context))

        self._add_modis_l1b_geo_input(context, task)
        self._add_modis_l1b_m_input(context, task)
        self._add_airs_l1b_input(context, task)

        # Make the product definition available to build_task()
        task.option('product', product)

    def _add_viirs_l1b_geo_input(self, product, context, task):
        satellite = context['satellite']
        granule = context['granule']
        version = product.input('viirs_l1').version
        input_name = sipsprod.satellite_esdt('V03MOD', satellite)
        LOG.debug("Ingesting input {} ({}) for V02FSN version {}".format(input_name, version, product.version))
        vgeom = dawg_catalog.file(satellite, input_name, granule, version=version)
        task.input('geo', vgeom)

    def _add_viirs_l1b_m_input(self, product, context, task):
        satellite = context['satellite']
        granule = context['granule']
        version = product.input('viirs_l1').version
        input_name = sipsprod.satellite_esdt('V02MOD', satellite)
        LOG.debug("Ingesting input {} ({}) for V02FSN version {}".format(input_name, version, product.version))
        vl1b = dawg_catalog.file(satellite, input_name, granule, version=version)
        task.input('l1b', vl1b)

    def _add_cris_l1b_input(self, product, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=6)
        cris_interval = TimeInterval(granule, granule+granule_length-timedelta(seconds=1))
        version = product.input('cris_l1').version
        input_name = 'CL1B'
        LOG.debug("Ingesting input {} ({}) for V02FSN version {}".format(input_name, version, product.version))
        cris = dawg_catalog.files(satellite, input_name, cris_interval, version=version)
        if cris == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                input_name, version, cris_interval))

        for idx, cris_file in enumerate(cris):
            LOG.debug('CrIS granule {}: {} -> {}'.format(idx, cris_file.begin_time, cris_file.end_time))
            task.input('sounder_{}'.format(idx),  cris_file)

    def build_task_snpp(self, context, task):
        '''
        Build up a set of inputs for a single context
        '''
        LOG.debug("Ingesting inputs for V02FSN version {} ...".format(context['version']))

        # Get the product definition for 'V02FSN'
        product = sipsprod.lookup_product_recurse('V02FSN', version=context['version'])

        # Ingest the required inputs, defined in the VNP02 product definition for context['version']
        self._add_viirs_l1b_geo_input(product, context, task)
        self._add_viirs_l1b_m_input(product, context, task)
        self._add_cris_l1b_input(product, context, task)

        # Make the product definition available to build_task()
        task.option('product', product)

    @reraise_as(WorkflowNotReady, FileNotFound, prefix='V02FSN')
    def build_task(self, context, task):

        satellite = context['satellite']

        if satellite == 'snpp':
            self.build_task_snpp(context, task)
        elif satellite == 'aqua':
            self.build_task_aqua(context, task)
        else:
            raise ValueError('Invalid satellite: {}'.format(satellite))

    def mend_viirs_l1b(self, product, geo, l1b, dummy=False):

        out_fn = basename(l1b).replace('.nc', '.bowtie_restored.nc')
        LOG.info('out_fn = {}'.format(out_fn))
        shutil.copy(l1b, out_fn)

        version = product.input('viirsmend').version
        viirsmend = support_software.lookup('viirsmend', version=version)
        py_exe = pjoin(viirsmend.path,'bin/python')
        viirsmend_exe = pjoin(viirsmend.path,'bin/viirsl1mend')

        cmd = '{} {} {} {}'.format(py_exe, viirsmend_exe, out_fn, geo)
        LOG.debug('cmd = {}'.format(cmd))

        dummy_bowtie_rest_file = '/mnt/sdata/geoffc/fusion_matlab/work/snpp_temp_outputs/outputs/tmp5PwThb/VNP02MOD.A2018033.1836.001.2018033235822.uwssec.bowtie_restored.nc'

        if not dummy:
            runscript(cmd, requirements=[viirsmend])
        else:
            LOG.debug('dummy cmd = "cp {} {}"'.format(dummy_bowtie_rest_file, out_fn))
            shutil.copy(dummy_bowtie_rest_file, out_fn)

        return out_fn

    def cris_viirs_collocation(self, product, inputs, dummy=False):

        LOG.info('inputs = {}'.format(inputs))
        input_keys = inputs.keys()
        input_keys.sort()

        sounder_keys = [key for key in input_keys if 'sounder' in key]
        cris_files = [inputs[key] for key in sounder_keys]

        vgeom_file = inputs['geo']

        version = product.input('collopak').version
        crisviirs = support_software.lookup('collopak', version=version)
        crisviirs_exe = pjoin(crisviirs.path,'bin/crisviirs')

        dummy_collo_file = '/mnt/sdata/geoffc/fusion_matlab/work/snpp_temp_outputs/outputs/tmp5PwThb/colloc.cris_snpp.viirs_m_snpp.20180202T183600_183600.nc'

        if not dummy:
            for cris_file in cris_files:
                cmd = '{} {} {} > /dev/null'.format(crisviirs_exe, cris_file, vgeom_file)

                LOG.info('cmd = {}'.format(cmd))
                runscript(cmd, requirements=[])
        else:
            LOG.debug('dummy cmd = "cp {} {}"'.format(dummy_collo_file, abspath(curdir)))
            shutil.copy(dummy_collo_file, './')

        collo_files = glob('colloc.*.nc')
        collo_files.sort()

        return collo_files

    def airs_modis_collocation(self, product, inputs, dummy=False):

        LOG.info('inputs = {}'.format(inputs))
        input_keys = inputs.keys()
        input_keys.sort()

        sounder_keys = [key for key in input_keys if 'sounder' in key]
        airs_files = [inputs[key] for key in sounder_keys]

        modis_file = inputs['geo']

        version = product.input('collopak').version
        airsmodis = support_software.lookup('collopak', version=version)
        airsmodis_exe = pjoin(airsmodis.path,'bin/airsmod')

        dummy_collo_files = glob('/mnt/sdata/geoffc/fusion_matlab/work/aqua_temp_outputs/outputs/tmpMUoHF3/colloc.airs_aqua.modis_aqua.*.nc')

        if not dummy:
            for airs_file in airs_files:
                cmd = '{} {} {} > /dev/null'.format(airsmodis_exe, airs_file, modis_file)

                LOG.info('cmd = {}'.format(cmd))
                runscript(cmd, requirements=[])
        else:
            LOG.debug('dummy cmd = "cp {} {}"'.format(' '.join(dummy_collo_files), abspath(curdir)))
            for dummy_collo_file in dummy_collo_files:
                shutil.copy(dummy_collo_file, './')

        collo_files = glob('colloc.*.nc')
        collo_files.sort()

        return collo_files

    def run_fusion_matlab(self, product, geo_file, l1b_file, sounder_files, collo_files, **kwargs):
        '''
        Run the Matlab fusion binary on the input level-1b files to generate the *.mat output file.
        '''

        bin_dir = kwargs['bin_dir']
        anc_paths = kwargs['anc_paths']
        fusion_binary = kwargs['fusion_binary']
        fused_output_dir = kwargs['fused_output_dir']
        matlab_file_glob = kwargs['matlab_file_glob']
        env = kwargs['env']
        granule = kwargs['granule']
        satellite = kwargs['satellite']

        dummy = kwargs['dummy']
        if dummy:
            if satellite=='snpp':
                dummy_matlab_file = '/mnt/sdata/geoffc/fusion_matlab/work/snpp_temp_outputs/outputs/tmp5PwThb/fusion_output.mat'
            if satellite=='aqua':
                dummy_matlab_file = '/mnt/sdata/geoffc/fusion_matlab/work/aqua_temp_outputs/outputs/tmpMUoHF3/fusion_output.mat'

        rc_fusion = 0

        # Get the matlab runtim version that we require
        matlab_version = '2015b'
        #matlab_version = product.input('fusion_matlab').options['matlab_version']

        #run matlab
        cmd = '{}/{} {} {} {} {} {} {}'.format(
            bin_dir,
            fusion_binary,
            support_software.lookup('matlab', matlab_version).path,
            geo_file,
            l1b_file,
            ' '.join(sounder_files),
            ' '.join(collo_files),
            ' '.join(anc_paths)
            )

        # Create the fused outputs directory in the temp work dir (which should be the current dir)
        tmp_work_dir = os.getcwd()
        LOG.debug('The current directory is {}'.format(tmp_work_dir))
        create_dir(fused_output_dir)

        # Run the Matlab Fusion code
        try:
            LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))
            rc_fusion = 0
            if not dummy:
                runscript(cmd, requirements=[], env=env)
            else:
                LOG.debug('dummy cmd = "cp {} {}"'.format(dummy_matlab_file, tmp_work_dir))
                shutil.copy(dummy_matlab_file, tmp_work_dir) # DEBUG
        except CalledProcessError as err:
            rc_fusion = err.returncode
            LOG.error("Matlab binary {} returned a value of {}".format(fusion_binary, rc_fusion))
            return rc_fusion, []

        # Move matlab file to the fused outputs directory
        matlab_file = glob(matlab_file_glob)
        if len(matlab_file) != 0:
            matlab_file = matlab_file[0]
            LOG.debug('Found Matlab file "{}", moving to {}...'.format(matlab_file, fused_output_dir))
            if exists(pjoin(fused_output_dir, matlab_file)):
                LOG.info('{} exists, removing...'.format(pjoin(fused_output_dir, matlab_file)))
                os.remove(pjoin(fused_output_dir, matlab_file))
            shutil.move(matlab_file, fused_output_dir)
            matlab_file = glob(pjoin(fused_output_dir, matlab_file))[0]
        else:
            LOG.error('There are no Matlab files "{}" to convert, aborting'.format(matlab_file_glob))
            rc_fusion = 1
            return rc_fusion, []

        return rc_fusion, matlab_file


    def convert_matlab_to_netcdf(self, product, matlab_file, l1b_file, **kwargs):
        '''
        Transcode the Matlab *.mat file into a CF-compliant HDF4 (aqua) or NetCDF4 (snpp) file.
        '''

        py_interp = kwargs['py_interp']
        bin_dir = kwargs['bin_dir']
        fused_output_dir = kwargs['fused_output_dir']
        conversion_bin = kwargs['conversion_bin']
        env = kwargs['env']
        satellite = kwargs['satellite']
        granule = kwargs['granule']

        dummy = kwargs['dummy']

        if dummy:
            if satellite=='snpp':
                dummy_fusion_file = '/mnt/sdata/geoffc/fusion_matlab/work/snpp_temp_outputs/outputs/tmp5PwThb/VNP02FSN.A2018033.1836.001.2018058173216.nc'
            if satellite=='aqua':
                dummy_fusion_file = '/mnt/sdata/geoffc/fusion_matlab/work/aqua_temp_outputs/outputs/tmpMUoHF3/MYD02FSN.A2015107.1755.006.2018058170733.hdf'

        rc_fusion = 0
        dt_create = datetime.utcnow()

        # Create the output directory
        tmp_work_dir = os.getcwd()
        LOG.debug('tmp_work_dir (CWD): "{}"'.format(tmp_work_dir))
        LOG.debug('fused_output_dir: "{}"'.format(fused_output_dir))

        # Copy the un-fused level-1b file to the work directory as a template...
        unfused_l1b_file = pjoin(tmp_work_dir, 'unfused', basename(l1b_file))
        unfused_l1b_dir = dirname(unfused_l1b_file)
        create_dir(unfused_l1b_dir)

        if exists(unfused_l1b_file):
            LOG.debug('{} exists, removing...'.format(unfused_l1b_file))
            os.remove(unfused_l1b_file)

        LOG.debug('Copying {} to {}'.format(l1b_file, unfused_l1b_file))
        shutil.copy(l1b_file, unfused_l1b_file)

        # Removing the fused file if it exists
        fused_l1b_file = pjoin(fused_output_dir, basename(unfused_l1b_file))
        LOG.debug('Checking for existing fused output file "{}"'.format(fused_l1b_file))
        if exists(fused_l1b_file):
            LOG.debug('{} exists, removing...'.format(fused_l1b_file))
            os.remove(fused_l1b_file)

        # Convert the Matlab file to the desired format...
        cmd = '{} {}  {} {} {}'.format(
                py_interp,
                conversion_bin,
                unfused_l1b_file,
                matlab_file,
                fused_output_dir
                )
        try:
            LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))
            rc_fusion = 0
            if not dummy:
                runscript(cmd, requirements=[], env=env)
            else:
                LOG.debug('dummy cmd = "cp {} {}"'.format(dummy_fusion_file,
                                                          pjoin(fused_output_dir, basename(l1b_file))))
                shutil.copy(dummy_fusion_file, pjoin(fused_output_dir, basename(l1b_file))) # DEBUG
        except CalledProcessError as err:
            rc_fusion = err.returncode
            LOG.error("CF converter {} returned a value of {}".format(conversion_bin, rc_fusion))
            return rc_fusion, []

        # Determine success...
        LOG.debug('Looking for fused output file "{}"'.format(fused_l1b_file))
        fused_l1b_file = glob(fused_l1b_file)
        if len(fused_l1b_file) != 0:
            fused_l1b_file = fused_l1b_file[0]
            LOG.debug('Found final fused output file "{}"'.format(fused_l1b_file))
        else:
            LOG.error('There is no fused file {}, aborting'.format(fused_l1b_file))
            rc_fusion = 1
            return rc_fusion, []

        # Remove the unfused dir...
        #LOG.debug('Removing the unfused level-1b dir {} ...'.format(unfused_l1b_dir))
        #shutil.rmtree(unfused_l1b_dir)

        # Determine the name of the output fused file
        if satellite=='snpp' or satellite=='jpss1':
            esdt = sipsprod.satellite_esdt('V02FSN', satellite)
            product.options['collection'] = int(basename(l1b_file).split('.')[3])
            fused_l1b_file_new = sipsprod.product_filename(esdt, product.options['collection'],
                                                  granule, dt_create)
        if satellite=='aqua':
            esdt = sipsprod.satellite_esdt('M02FSN', satellite)
            product.options['collection'] = int(basename(l1b_file).split('.')[3])
            fused_l1b_file_new = sipsprod.product_filename(esdt, product.options['collection'],
                                                  granule, dt_create)
            fused_l1b_file_new = fused_l1b_file_new.replace('.nc', '.hdf')

        # Move the HDF4/NetCDF4 file to its new filename
        LOG.debug('Moving "{}" to "{}" ...'.format(fused_l1b_file, pjoin(tmp_work_dir, fused_l1b_file_new)))
        shutil.move(fused_l1b_file, pjoin(tmp_work_dir, fused_l1b_file_new))
        fused_l1b_file = glob(pjoin(tmp_work_dir, fused_l1b_file_new))[0]

        # Move the matlab file to its new filename
        if satellite=='snpp' or satellite=='jpss1':
            matlab_file_new = basename(fused_l1b_file_new).replace('.nc','.mat')
        if satellite=='aqua':
            matlab_file_new = basename(fused_l1b_file_new).replace('.hdf','.mat')

        LOG.debug('Moving "{}" to {}...'.format(matlab_file, pjoin(fused_output_dir, matlab_file_new)))
        shutil.move(matlab_file, pjoin(fused_output_dir, matlab_file_new))
        matlab_file = glob(pjoin(fused_output_dir, matlab_file_new))[0]

        # Remove the fused_outputs directory
        #LOG.debug('Removing the fused_outputs dir {} ...'.format(fused_output_dir))
        #shutil.rmtree(fused_output_dir)

        output_attrs = {'esdt': esdt, 'collection': product.options['collection'],
                        'created': dt_create}

        return rc_fusion, fused_l1b_file, output_attrs

    def output_QC(self, l1b_file, fused_l1b_file, band=None, input_rms=0.2, **kwargs):

        satellite = kwargs['satellite']
        LOG.debug('satellite = {}'.format(satellite))

        band_default = {'aqua':[31, 32], 'snpp':[15, 16]}
        if band is None:
            band = band_default[satellite]

        input_files = [fused_l1b_file] if satellite=='snpp' else [l1b_file, fused_l1b_file]

        generators = list(SFX[splitext(p)[-1].lower()](p, band) for p in input_files)
        stuff = list(chain(*generators))
        if len(stuff) != 2:
            raise AssertionError("requires original + fusion input")
        passfail, rms = detect(input_rms, *stuff)

        if passfail:
            LOG.info("Output QC PASS (rms {} < threshold {})".format(rms, input_rms))
        else:
            LOG.warn("Output QC FAIL (rms {} > threshold {})".format(rms, input_rms))

        return 0 if passfail else 1


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

    @reraise_as(WorkflowNotReady, FileNotFound, prefix='V02FSN')
    def run_task(self, inputs, context):

        LOG.debug("Running run_task()...")

        for key in context.keys():
            LOG.debug("run_task() context['{}'] = {}".format(key, context[key]))

        granule = context['granule']
        satellite = context['satellite']
        version = context['version']

        # Get the location of the binary package
        product = context['product']
        delivery = delivered_software.lookup(
                'fusion_matlab', delivery_id=product.input('fusion_matlab').version)
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

        # What are our inputs?
        for input in inputs.keys():
            inputs_dir = dirname(inputs[input])
            LOG.debug("inputs['{}'] = {}".format(input,inputs[input]))
        LOG.debug("Inputs dir = {}".format(inputs_dir))

        # Are we doing a dummy run?
        #dummy = True
        dummy = False

        # Run viirsmend on the viirs level-1b, and generate the CrIS/VIIRS collocation
        if satellite == 'snpp':
            geo = inputs['geo']
            l1b = self.mend_viirs_l1b(product, inputs['geo'], inputs['l1b'], dummy=dummy)
            sounder_keys = [key for key in inputs.keys() if 'sounder' in key]
            sounder = [inputs[key] for key in sounder_keys]
            collo = self.cris_viirs_collocation(product, inputs, dummy=dummy)

        # Generate the AIRS/MODIS collocation
        if satellite == 'aqua':
            geo = inputs['geo']
            l1b = inputs['l1b']
            sounder_keys = [key for key in inputs.keys() if 'sounder' in key]
            sounder = [inputs[key] for key in sounder_keys]
            collo = self.airs_modis_collocation(product, inputs, dummy=dummy)

        LOG.info('geo = {}'.format(geo))
        LOG.info('l1b = {}'.format(l1b))
        LOG.info('sounder = {}'.format(sounder))
        LOG.info('collo = {}'.format(collo))

        bin_dir = pjoin(dist_root, 'bin')
        anc_dir = pjoin(dist_root, 'luts')
        fused_output_dir = pjoin(work_dir, 'fused_outputs')

        # Setup the require keyword arguments for the fusion_matlab package
        kwargs = {}
        kwargs['py_interp'] = py_interp
        kwargs['bin_dir'] = bin_dir
        kwargs['env'] = env
        kwargs['fused_output_dir'] = fused_output_dir
        kwargs['satellite'] = satellite
        kwargs['granule'] = granule
        kwargs['dummy'] = dummy

        if satellite=='snpp':
            kwargs['anc_paths'] = [pjoin(anc_dir, 'modis_aqua.srf.nc'),
                                   pjoin(anc_dir, 'NG_VIIRS_NPP_RSR_filtered_Oct2011_BA/')]
            kwargs['fusion_binary'] = 'run_imagersounderfusion_V.sh'
            kwargs['matlab_file_glob'] = 'fusion_output.mat'
            kwargs['conversion_bin'] = pjoin(envroot, 'bin', 'l1b-fusion-viirs-cris')
        elif satellite=='aqua':
            kwargs['anc_paths'] = [pjoin(anc_dir, 'L2.chan_prop.2005.03.01.v9.5.1.txt'),
                                   pjoin(anc_dir, 'modis_aqua.srf.nc'),
                                   pjoin(anc_dir, 'modis_conv_error_2005.mat')]
            kwargs['fusion_binary'] = 'run_imagersounderfusion_M.sh'
            kwargs['matlab_file_glob'] = 'fusion_output.mat'
            kwargs['conversion_bin'] = pjoin(envroot, 'bin', 'l1b-fusion-modis-airs')
        else:
            return {}

        # Run the fusion_matlab package
        rc_fusion, matlab_file = self.run_fusion_matlab(
                                                   product,
                                                   geo,
                                                   l1b,
                                                   sounder,
                                                   collo,
                                                   **kwargs
                                                  )

        LOG.debug('run_fusion_matlab() return value: {}'.format(rc_fusion))
        LOG.info('run_fusion_matlab() generated {}'.format(matlab_file))

        # Now that we've computed the Matlab file, convert to a NetCDF file...
        rc_fusion, fused_l1b_file, output_attrs = self.convert_matlab_to_netcdf(
                                                                               product,
                                                                               matlab_file,
                                                                               l1b,
                                                                               **kwargs)

        LOG.debug('convert_matlab_to_netcdf() return value: {}'.format(rc_fusion))
        LOG.info('convert_matlab_to_netcdf() generated {}'.format(fused_l1b_file))

        # Update some global attributes in the output file
        #readme_file =  pjoin(delivery.path, 'README.txt')
        #self.update_global_attrs(basename(fused_l1b_file), readme_file, **kwargs)

        # Run a QC check on the output file
        #rc_qc = self.output_QC(l1b, fused_l1b_file, **kwargs)
        #LOG.debug('output_QC() return value: {}'.format(rc_qc))
        #if rc_qc != 0:
            #raise RuntimeError('Output fusion file {} failed RMS error QC check, output aborted.'.format(fused_l1b_file))

        # The staging routine assumes that the output file is located in the work directory
        # "tmp******", and that the output path is to be prepended, so return the basename.

        out_fn = basename(fused_l1b_file)

        # Set metadata to be put in the output file.
        if satellite == 'snpp':
            l1_version = product.input('viirs_l1').version
            input_fns, lut_version, lut_created = get_viirs_l1_luts(l1b, geo_fn=geo)
            ancillary_fns = []
            out_compress = nc_compress

            set_official_product_metadata(
                output_attrs['esdt'], product.version, output_attrs['collection'],
                product.input('fusion_matlab').version, context['satellite'],
                fused_l1b_file, geo,
                input_fns, ancillary_fns,
                l1_version, lut_version, lut_created,
                product.inputstr(), output_attrs['created'])

        if satellite == 'aqua':
            out_compress = hdf_compress
            #l1_version = 'ingest'
            #input_fns = [l1b, geo]
            #lut_version = ''
            #lut_created = granule
            #ancillary_fns = []


        LOG.debug('We are in {}'.format(os.getcwd()))
        LOG.debug('Compressing {}'.format(out_fn))
        return {'fused_l1b': out_compress(out_fn)}
        #return {'fused_l1b': out_fn}

        #return {
            #'fused_l1b': {
                #'file': out_fn,
                #'extra_attrs': extra_attrs,
            #},
        #}
