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
import netCDF4
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
    hdf_compress,
    reraise_as,
    set_official_product_metadata,
    FileNotFound
)
from glutil.nc import (nc_gen, nc_compress, nc_remove_unlimited_dims, nc_setattrs, nc_getattrs)

from utils import create_dir
from detect_bad_fusion import SFX, detect

# every module should have a LOG object
LOG = logging.getLogger(__name__)

class FusionProcessFailed(Exception):
    exit_code = 6001
class FusionFailedToProduceMat(Exception):
    exit_code = 6002
class CFConversionFailed(Exception):
    exit_code = 6003
class CFConversionFailedToProductNC(Exception):
    exit_code = 6004
class NCL2MetadataFailed(Exception):
    exit_code = 6005

class FUSION_MATLAB(Computation):

    parameters = ['granule', 'satellite', 'version']

    outputs = ['fused_l1b']

    def find_contexts(self, time_interval, satellite, version):
        '''
        Here we assume that the granule boundaries fall along 6-minute (snpp/noaa20) or
        5-minute (aqua) increments, starting at the top of the hour:

        SNPP:  [0.,   6.,  12.,  18.,  24.,  30.,  36.,  42.,  48.,  54.]
        NOAA20:[0.,   6.,  12.,  18.,  24.,  30.,  36.,  42.,  48.,  54.]
        AQUA:  [0.,   5.,  10.,  15.,  20.,  25.,  30.,  35.,  40.,  45.,  50., 55.]
        '''

        if satellite=='snpp' or satellite=='noaa20':
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
        wedge = timedelta(seconds=1)

        interval = TimeInterval(granule, granule+granule_length-wedge)

        version = 'c6'
        input_name = 'MOD03'

        LOG.debug("Ingesting input {} ({})...".format(input_name, version))
        mgeo = dawg_catalog.files(satellite, input_name, interval, version=version)
        if mgeo == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                input_name, version, interval))

        for idx, mgeo_file in enumerate(mgeo):
            LOG.debug('MODIS GEO granule {}: {} -> {}'.format(idx, mgeo_file.begin_time, mgeo_file.end_time))
            task.input('geo_{}'.format(idx),  mgeo_file)

    def _add_modis_l1b_m_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=5)
        wedge = timedelta(seconds=1)

        interval = TimeInterval(granule, granule+granule_length-wedge)

        version = 'c6'
        input_name = 'MOD021KM'

        LOG.debug("Ingesting input {} ({})...".format(input_name, version))
        myd021km = dawg_catalog.files(satellite, input_name, interval, version=version)
        if myd021km == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                input_name, version, interval))

        for idx, myd021km_file in enumerate(myd021km):
            LOG.debug('MODIS L1B granule {}: {} -> {}'.format(idx, myd021km_file.begin_time, myd021km_file.end_time))
            task.input('l1b_{}'.format(idx),  myd021km_file)

    def _add_airs_l1b_input(self, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=5)
        wedge = timedelta(seconds=1)

        myd03_interval = TimeInterval(granule, granule+granule_length-wedge)

        version = 'c6'
        input_name = 'MOD03'

        myd03 = dawg_catalog.files(satellite, input_name, myd03_interval, version=version)
        if myd03 == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                input_name, version, myd03_interval))

        myd03 = sorted(myd03, key=lambda x: basename(str(x.begin_time)))

        buf = 1 # seconds
        # airs_begin = myd03[0].begin_time - timedelta(seconds=buf)
        # airs_end = myd03[-1].end_time + timedelta(seconds=buf)
        airs_begin = myd03[0].begin_time - granule_length #+ timedelta(seconds=buf)
        airs_end = myd03[-1].end_time + granule_length #- timedelta(seconds=buf)
        airs_interval = TimeInterval(airs_begin, airs_end)

        input_name = 'AIRIBRAD'

        LOG.debug("Ingesting input {}...".format(input_name))
        airs = dawg_catalog.files(satellite, 'AIRIBRAD', airs_interval)
        if airs == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for interval {}'.format(
                input_name, airs_interval))

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
        granule_length = timedelta(minutes=6)

        viirs_l1 = product.input('viirs_l1')
        v02mod_bt_sc = product.input('V02MOD-bt-sc')

        interval = TimeInterval(granule-granule_length, granule+granule_length)

        if viirs_l1:
            version = viirs_l1.version
        elif v02mod_bt_sc:
            version = v02mod_bt_sc.input('viirs_l1').version

        input_name = sipsprod.satellite_esdt('V03MOD', satellite)

        LOG.debug("Ingesting input {} ({}) for FSNRAD_L2_VIIRS_CRIS version {}".format(input_name, version, product.version))
        vgeom = dawg_catalog.files(satellite, input_name, interval, version=version)

        if vgeom == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                input_name, version, interval))
        if len(vgeom) != 3:
            LOG.warn('Insufficient VIIRS GEO granules (should be 3 granules) for {}...'.format(granule))
            for idx, vgeom_file in enumerate(vgeom):
                LOG.debug('\tVIIRS GEO granule {}: {} -> {}'.format(idx, vgeom_file.begin_time, vgeom_file.end_time))
            raise WorkflowNotReady('FUSION_MATLAB: Insufficient {} inputs for version {} and interval {}'.format(
                input_name, version, interval))

        for idx, vgeom_file in enumerate(vgeom):
            LOG.debug('VIIRS GEO granule {}: {} -> {}'.format(idx, vgeom_file.begin_time, vgeom_file.end_time))
            task.input('geo_{}'.format(idx),  vgeom_file)

    def _add_viirs_l1b_m_input(self, product, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=6)

        viirs_l1 = product.input('viirs_l1')
        v02mod_bt_sc = product.input('V02MOD-bt-sc')

        interval = TimeInterval(granule-granule_length, granule+granule_length)

        # Standard NASA VIIRS level-1 radiances and BT
        if viirs_l1:
            version = viirs_l1.version
            input_name = sipsprod.satellite_esdt('V02MOD', satellite)
            LOG.debug("Ingesting input {} ({}) for FSNRAD_L2_VIIRS_CRIS version {}".format(input_name, version, product.version))
            vl1b = dawg_catalog.files(satellite, input_name, interval, version=version)
            if vl1b == []:
                raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                    input_name, version, interval))

        # Standard NASA VIIRS level-1 radiances and BT, bias corrected to match MODIS (or something...)
        elif v02mod_bt_sc:
            version = v02mod_bt_sc.version
            input_name = sipsprod.satellite_esdt('V02MOD-bt-sc', satellite)
            LOG.debug("Ingesting input {} ({}) for FSNRAD_L2_VIIRS_CRIS version {}".format(input_name, version, product.version))
            from flo.sw.v02mod_bt_sc import v02mod_bt_sc
            vl1b = v02mod_bt_sc().dataset('out').product(
                {'granule': context['granule'], 'satellite': context['satellite'],
                 'nrt': False, 'version': version})
        else:
            raise ValueError('V02MOD -- missing one of (V02MOD, V02MOD-bc) input for VIIRS')

        if vl1b == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                input_name, version, interval))
        if len(vl1b) != 3:
            LOG.warn('Insufficient VIIRS L1B granules (should be 3 granules) for  {}...'.format(granule))
            for idx, vl1b_file in enumerate(vl1b):
                LOG.debug('\tVIIRS L1B granule {}: {} -> {}'.format(idx, vl1b_file.begin_time, vl1b_file.end_time))
            raise WorkflowNotReady('FUSION_MATLAB: Insufficient {} inputs for version {} and interval {}'.format(
                input_name, version, interval))

        for idx, vl1b_file in enumerate(vl1b):
            LOG.debug('VIIRS L1B granule {}: {} -> {}'.format(idx, vl1b_file.begin_time, vl1b_file.end_time))
            task.input('l1b_{}'.format(idx),  vl1b_file)

    def get_cris_l1b_version(self, product, context):
        '''
        Return the correct CrIS version and input name depending on the target granule, as the
        required version and input name may differ from that given in the product definition for a
        given granule.
        '''

        satellite = context['satellite']
        granule = context['granule']

        input_name = 'CL1B'
        version = product.input('cris_l1').version

        if satellite == 'snpp':

            while True:
                # Use CL1B_nsr (2.0.15) before dt_switch_nsr
                dt_switch_nsr = datetime(2015, 12, 1, 0)
                if (granule - dt_switch_nsr).total_seconds() <= 0.:
                    input_name = 'CL1B_nsr'
                    version = '2.0.15'
                    break

                # Use CL1B (2.0.15) before dt_switch_side2
                dt_switch_side2 = datetime(2019, 6, 24, 18, 48)
                if (granule - dt_switch_side2).total_seconds() <= 0.:
                    input_name = 'CL1B'
                    version = '2.0.15'
                    break

                # We're after the special cases, use what's in the product definition
                break

        return input_name, version

    def _add_cris_l1b_input(self, product, context, task):
        satellite = context['satellite']
        granule = context['granule']
        granule_length = timedelta(minutes=6)

        interval = TimeInterval(granule-granule_length, granule+granule_length)

        input_name, version = self.get_cris_l1b_version(product, context)

        LOG.debug("Ingesting input {} ({}) for FSNRAD_L2_VIIRS_CRIS version {}".format(input_name, version, product.version))
        cris = dawg_catalog.files(satellite, input_name, interval, version=version)

        if cris == []:
            raise WorkflowNotReady('FUSION_MATLAB: Missing {} inputs for version {} and interval {}'.format(
                input_name, version, interval))
        if len(cris) != 3:
            LOG.warn('Insufficient CrIS L1B granules (should be 3 granules) for {}...'.format(granule))
            for idx, cris_file in enumerate(cris):
                LOG.debug('\tCrIS L1B granule {}: {} -> {}'.format(idx, cris_file.begin_time, cris_file.end_time))
            raise WorkflowNotReady('FUSION_MATLAB: Insufficient {} inputs for version {} and interval {}'.format(
                input_name, version, interval))

        for idx, cris_file in enumerate(cris):
            LOG.debug('CrIS L1B granule {}: {} -> {}'.format(idx, cris_file.begin_time, cris_file.end_time))
            task.input('sounder_{}'.format(idx),  cris_file)

    def build_task_snpp(self, context, task):
        '''
        Build up a set of inputs for a single context
        '''
        LOG.debug("Ingesting inputs for FSNRAD_L2_VIIRS_CRIS version {} ...".format(context['version']))

        # Get the product definition for 'FSNRAD_L2_VIIRS_CRIS'. Different versions may use regular or bias
        # corrected VIIRS level-1b files.
        product_name = 'FSNRAD_L2_VIIRS_CRIS'
        product = sipsprod.lookup_product_recurse(product_name, version=context['version'])

        # Ingest the required inputs, defined in the VNP02 product definition for context['version']
        self._add_viirs_l1b_geo_input(product, context, task)
        self._add_viirs_l1b_m_input(product, context, task)
        self._add_cris_l1b_input(product, context, task)

        # Make the product definition available to build_task()
        task.option('product', product)

    @reraise_as(WorkflowNotReady, FileNotFound, prefix='FSNRAD_L2_VIIRS_CRIS')
    def build_task(self, context, task):

        satellite = context['satellite']

        if satellite=='snpp' or satellite=='noaa20':
            self.build_task_snpp(context, task)
        elif satellite == 'aqua':
            self.build_task_aqua(context, task)
        else:
            raise ValueError('Invalid satellite: {}'.format(satellite))

    def mend_viirs_l1b(self, product, geo, l1b, dummy=False):

        out_fn = basename(l1b).replace('.nc', '.bowtie_restored.nc')
        LOG.debug('out_fn = {}'.format(out_fn))
        shutil.copy(l1b, out_fn)

        version = product.input('viirsmend').version
        viirsmend = support_software.lookup('viirsmend', version=version)
        py_exe = pjoin(viirsmend.path,'bin/python')
        viirsmend_exe = pjoin(viirsmend.path,'bin/viirsl1mend')

        cmd = '{} {} {} {}'.format(py_exe, viirsmend_exe, out_fn, geo)

        dummy_bowtie_rest_file = sorted(glob(pjoin('/data/geoffc/fusion_matlab/work/local_processing',
                                                   'snpp_fusion_output/outputs',
                                                   'VNP02MOD*bowtie*.nc')))

        if not dummy:
            runscript(cmd, requirements=[viirsmend])
        else:
            LOG.debug('dummy cmd = "cp {} {}"'.format(dummy_bowtie_rest_file, out_fn))
            rc_copy = [shutil.copy(bowtie_file, out_fn) for bowtie_file in dummy_bowtie_rest_file]
            LOG.debug('rc_copy =  {}'.format(rc_copy))

        return out_fn

    def cris_viirs_collocation(self, product, context, geo, sounder, dummy=False):

        satellite = context['satellite']

        LOG.debug("geo = \n\t{}".format('\n\t'.join(geo)))
        LOG.debug("sounder = \n\t{}".format('\n\t'.join(sounder)))

        version = product.input('collopak').version
        crisviirs = support_software.lookup('collopak', version=version)
        crisviirs_exe = pjoin(crisviirs.path,'bin/crisviirs')

        dummy_collo_file = sorted(glob(pjoin('/data/geoffc/fusion_matlab/work/local_processing',
                                                   'snpp_fusion_output/outputs',
                                                   'colloc.cris_snpp.viirs_m_snpp.*.nc')))

        if not dummy:
            for cris, vgeo in zip(sounder, geo):
                collo_dt = datetime.strptime(basename(cris).split('.')[3],'%Y%m%dT%H%M')
                collo_log = 'colloc.cris_{0:}.viirs_m_{0:}.{1:%Y%m%dT%H%M%S}_{1:%H%M%S}'.format(
                        satellite, collo_dt)
                #cmd = '{} {} {} > /dev/null'.format(crisviirs_exe, cris, vgeo)
                cmd = '{} {} {} > {}.log'.format(crisviirs_exe, cris, vgeo, collo_log)

                runscript(cmd, requirements=[])
        else:
            LOG.debug('dummy cmd = "cp {} {}"'.format(dummy_collo_file, abspath(curdir)))
            rc_copy = [shutil.copy(collo_file, './') for collo_file in dummy_collo_file]
            LOG.debug('rc_copy =  {}'.format(rc_copy))

        collo_files = sorted(glob('colloc.*.nc'))

        return collo_files

    def airs_modis_collocation(self, product, inputs, dummy=False):

        LOG.debug('inputs = {}'.format(inputs))
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

                LOG.debug('cmd = {}'.format(cmd))
                runscript(cmd, requirements=[])
        else:
            LOG.debug('dummy cmd = "cp {} {}"'.format(' '.join(dummy_collo_files), abspath(curdir)))
            for dummy_collo_file in dummy_collo_files:
                shutil.copy(dummy_collo_file, './')

        collo_files = glob('colloc.*.nc')
        collo_files.sort()

        return collo_files

    def run_fusion_matlab(self, product, geo_files, l1b_files, sounder_files, collo_files, **kwargs):
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
            dummy_matlab_file = pjoin('/data/geoffc/fusion_matlab/work/local_processing',
                                                    'snpp_fusion_output/outputs/fusion_output.mat')

        rc_fusion = 0

        # Get the matlab runtim version that we require
        #matlab_version = '2015b'
        matlab_version = '2018b'
        #matlab_version = product.input('fusion_matlab').options['matlab_version']

        #run matlab
        cmd = '{}/{} {}/ {} {} {} {} {}'.format(
            bin_dir,
            fusion_binary,
            support_software.lookup('matlab', matlab_version).path,
            ' '.join(geo_files),
            ' '.join(l1b_files),
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
                shutil.copy(dummy_matlab_file, tmp_work_dir)
        except CalledProcessError as err:
            rc_fusion = err.returncode
            LOG.error("Matlab binary {} returned a value of {}".format(fusion_binary, rc_fusion))
            raise FusionProcessFailed('Matlab binary {} failed or was killed, returning a value of {}'.format(fusion_binary, rc_fusion))

        # Move matlab file to the fused outputs directory
        matlab_file = glob(matlab_file_glob)
        if len(matlab_file) != 0:
            matlab_file = matlab_file[0]
            LOG.debug('Found Matlab file "{}", moving to {}...'.format(matlab_file, fused_output_dir))
            if exists(pjoin(fused_output_dir, matlab_file)):
                LOG.debug('{} exists, removing...'.format(pjoin(fused_output_dir, matlab_file)))
                os.remove(pjoin(fused_output_dir, matlab_file))
            shutil.move(matlab_file, fused_output_dir)
            matlab_file = glob(pjoin(fused_output_dir, matlab_file))[0]
        else:
            LOG.error('There are no Matlab files "{}" to convert, aborting'.format(matlab_file_glob))
            raise FusionFailedToProduceMat('Matlab binary {} failed to produce the file "{}"'.format(fusion_binary, matlab_file_glob))

        return rc_fusion, matlab_file

    def convert_matlab_to_viirs_netcdf(self, product, viirs_l1b_file, cris_l1b_file, matlab_file,  **kwargs):
        '''
        Transcode the Matlab *.mat file into a CF-compliant NetCDF4 file.
        '''

        py_interp = kwargs['py_interp']
        bin_dir = kwargs['bin_dir']
        cdl_dir = kwargs['cdl_dir']
        fused_output_dir = kwargs['fused_output_dir']
        conversion_bin = kwargs['conversion_bin']
        env = kwargs['env']
        satellite = kwargs['satellite']
        granule = kwargs['granule']

        dummy = kwargs['dummy']

        rc_fusion = 0
        dt_create = datetime.utcnow()

        # Create the output directory
        tmp_work_dir = os.getcwd()
        LOG.debug('tmp_work_dir (CWD): "{}"'.format(tmp_work_dir))
        LOG.debug('fused_output_dir: "{}"'.format(fused_output_dir))

        # Construct the candidate filename for the new level-2 file
        esdt = 'FSNRAD_L2_VIIRS_CRIS' + ('_SNPP' if satellite=='snpp' else '_NOAA20')
        collection = product.options['collection']
        viirs_fused_l1b_file = sipsprod.product_filename(esdt, collection, granule, created=dt_create)
        viirs_fused_l1b_file = pjoin(fused_output_dir, basename(viirs_fused_l1b_file))

        # Remove the output level-2 template file if it exists
        LOG.debug('Checking for existing fused output file "{}"'.format(viirs_fused_l1b_file))
        if exists(viirs_fused_l1b_file):
            LOG.debug('{} exists, removing...'.format(viirs_fused_l1b_file))
            os.remove(viirs_fused_l1b_file)

        if not dummy:
            # Create the new level2 template file and add a dimension from the level1 file
            cdl_file = pjoin(cdl_dir, 'fusion_cris_viirs.cdl')
            LOG.info('Creating template file {} from CDL file {}'.format(viirs_fused_l1b_file, cdl_file))

            try:
               nc_gen(cdl_file, viirs_fused_l1b_file)
            except CalledProcessError as err:
               rc = err.returncode
               LOG.error("ncgen returned a value of {}".format(rc))
               return rc, []

            # Open the l1b file and get the 'number_of_LUT_values' dimension
            nc_l1b = netCDF4.Dataset(viirs_l1b_file, 'r')
            dim = nc_l1b.dimensions['number_of_LUT_values']

            # Copy the 'number_of_LUT_values' dim to the level-2 template file
            nc_l2 = netCDF4.Dataset(viirs_fused_l1b_file, 'r+')
            nc_l2.createDimension(dim.name, size=dim.size)

            # Close the level1 and level2 files
            nc_l2.close()
            nc_l1b.close()

            # Copy the Matlab file data to the level-2 tamplate file...
            cmd = '{} {} {} {} {} {}'.format(
                                          py_interp,
                                          conversion_bin,
                                          viirs_l1b_file,
                                          cris_l1b_file,
                                          matlab_file,
                                          viirs_fused_l1b_file)

            try:
                LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))
                rc_fusion = 0
                runscript(cmd, requirements=[], env=env)
            except CalledProcessError as err:
                rc_fusion = err.returncode
                LOG.error("CF converter {} returned a value of {}".format(conversion_bin, rc_fusion))
                raise CFConversionFailed('CF converter {} failed or was killed, returning a value of {}'.format(conversion_bin, rc_fusion))
                #return rc_fusion, None, {}

        else:

            if satellite=='snpp':
                dummy_fusion_file = '/data/geoffc/fusion_matlab/work/local_processing/snpp_fusion_output/FSNRAD_L2_VIIRS_CRIS_SNPP.A2019218.1824.001.2019220192142.nc'
            if satellite=='noaa20':
                dummy_fusion_file = '/data/geoffc/fusion_matlab/work/local_processing/snpp_fusion_output/VNP02FSN.A2018033.1836.001.2018058173216.nc'

            LOG.debug('dummy cmd = "cp {} {}"'.format(dummy_fusion_file, fused_output_dir))
            shutil.copy(dummy_fusion_file, fused_output_dir)
            viirs_fused_l1b_file = pjoin(fused_output_dir, basename(dummy_fusion_file))


        # Determine success...
        LOG.debug('Looking for fused output file "{}"'.format(viirs_fused_l1b_file))
        viirs_fused_l1b_file = glob(viirs_fused_l1b_file)
        if len(viirs_fused_l1b_file) != 0:
            viirs_fused_l1b_file = viirs_fused_l1b_file[0]
            LOG.debug('Found final fused output file "{}"'.format(viirs_fused_l1b_file))
        else:
            LOG.error('There is no fused file {}, aborting'.format(viirs_fused_l1b_file))
            raise CFConversionFailedToProductNC('CF converter {} failed to produce the file "{}"'.format(conversion_bin, viirs_fused_l1b_file))
            #rc_fusion = 1
            #return rc_fusion, None, {}


        # Move the NetCDF4 file to its new filename
        LOG.debug('Moving "{}" to "{}" ...'.format(viirs_fused_l1b_file, tmp_work_dir))
        shutil.move(viirs_fused_l1b_file, tmp_work_dir)
        viirs_fused_l1b_file = glob(pjoin(tmp_work_dir, basename(viirs_fused_l1b_file)))[0]
        LOG.debug('Final fused output file "{}"'.format(viirs_fused_l1b_file))

        # Move the matlab file to its new filename
        matlab_file_new = viirs_fused_l1b_file.replace('.nc','.mat')

        LOG.debug('Moving "{}" to {}...'.format(matlab_file, matlab_file_new))
        shutil.move(matlab_file, matlab_file_new)
        matlab_file = glob(matlab_file_new)[0]
        LOG.debug('Final fused output matlab file "{}"'.format(matlab_file))

        # Remove the fused_outputs directory
        LOG.debug('Removing the fused_outputs dir {} ...'.format(fused_output_dir))
        shutil.rmtree(fused_output_dir)

        output_attrs = {'esdt': esdt, 'collection': product.options['collection'],
                        'created': dt_create}

        return rc_fusion, viirs_fused_l1b_file, output_attrs

    def output_QC(self, l1b_file, fused_l1b_file, band=None, input_rms=0.2, **kwargs):

        satellite = kwargs['satellite']
        LOG.debug('satellite = {}'.format(satellite))

        band_default = {'aqua':[31, 32], 'snpp':[15, 16], 'noaa20':[15, 16]}
        if band is None:
            band = band_default[satellite]

        input_files = [fused_l1b_file] if (satellite=='snpp' or satellite=='noaa20') else [l1b_file, fused_l1b_file]

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

    def set_l2_metadata(self, l2_file, viirs_l1b_file, viirs_geo_file, cris_l1b_file, product, context):
        '''
        Construct various metadata strings and write them to the level2 file
        '''
        LOG.debug('Updating the metadata in {}'.format(l2_file))

        context['nrt'] = True

        satname = 'SNPP' if context['satellite']=='snpp' else 'NOAA20'
        esdt = product.name + ('_SNPP' if context['satellite']=='snpp' else '_NOAA20')
        dt_create = datetime.strptime(splitext(basename(l2_file))[0].split('.')[-1], '%Y%j%H%M%S')
        print("Creation date is {}".format(dt_create))

        viirs_input_fns, viirs_lut_version, viirs_lut_created = get_viirs_l1_luts(viirs_l1b_file, geo_fn=viirs_geo_file)
        ancillary_fns = []
        viirs_l1_version = product.input('viirs_l1').version

        # Make any input version mods that aren't in the product definition...
        cris_input_name, cris_l1_version = self.get_cris_l1b_version(product, context)

        set_official_product_metadata(
            esdt,
            product.version,
            product.options['collection'],
            product.input('fusion_matlab').version,
            context['satellite'],
            l2_file,
            viirs_geo_file,
            viirs_input_fns,
            ancillary_fns,
            viirs_l1_version,
            viirs_lut_version,
            viirs_lut_created,
            product.inputstr(),
            dt_create,
            context['nrt'])

        # Remove any troublesome global attributes
        try:
            nc_l2 = netCDF4.Dataset(l2_file, 'a')
            xmlmetadata = nc_l2.getncattr('xmlmetadata')
            nc_l2.delncattr('xmlmetadata')
            nc_l2.delncattr('l1_version')
            nc_l2.delncattr('l1_lut_version')
            nc_l2.delncattr('l1_lut_created')
            nc_l2.delncattr('input_files')
            source = nc_l2.getncattr('source')
            nc_l2.delncattr('source')
        except Exception as err:
            LOG.warning('There was a problem removing nc_l2 attrs in {}'.format(l2_file))
            LOG.debug(traceback.format_exc())
            nc_l2.close()
            raise NCL2MetadataFailed('NetCDF4 global metadata update failed')

        nc_l2.close()

        nc_remove_unlimited_dims(basename(l2_file))

        cris_l1_version += ' (NSR)' if 'nsr' in cris_input_name else ''
        source = ', '.join(
                ['cris_l1 {}'.format(cris_l1_version) if 'cris_l1' in x else x for x in source.split(',')]
                ).replace('  ',' ')

        FIX_ATTRS = [
            ('title', '{0:} VIIRS+CrIS Fusion ({1:})'.format(satname, esdt)),
            ('platform', {'snpp':'Suomi-NPP', 'noaa20':'NOAA-20'}[context['satellite']]),
            ('instrument', 'VIIRS+CrIS'),
            ('conventions', 'CF-1.6, ACDD-1.3'),
            ('AlgorithmType', 'OPS'),
            ('long_name', '{} VIIRS+CrIS Fusion 6-Min L2 Swath 750m'.format(satname)),
            ('project', 'NASA Atmosphere Discipline'),
            ('creator_name', 'NASA Atmosphere SIPS'),
            ('processing_version', product.input('fusion_matlab').version),
            ('viirs_l1_version', viirs_l1_version),
            ('viirs_lut_version', viirs_lut_version),
            ('viirs_lut_created', str(viirs_lut_created.isoformat())),
            ('source', source),
            ('cris_l1_version', cris_l1_version),
            ('input_files', ', '.join([basename(x) for x in [viirs_geo_file, viirs_l1b_file, cris_l1b_file]])),
            ('xmlmetadata', xmlmetadata),
            ]
        nc_setattrs(l2_file, FIX_ATTRS)

        extra_attrs = {'esdt': esdt,
                       'collection': product.options['collection'],
                       'ecstype': 'SCIENCE',
                       'viirs_lut_version': viirs_lut_version,
                       'viirs_lut_created': viirs_lut_created,
                       'created': dt_create}


        return {
            'fused_l1b': {
                'file': basename(l2_file),
                'extra_attrs': extra_attrs,
            }
    }

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

    @reraise_as(WorkflowNotReady, FileNotFound, prefix='FSNRAD_L2_VIIRS_CRIS')
    def run_task(self, inputs, context):

        LOG.debug("Running run_task()...")

        for key in context.keys():
            LOG.debug("run_task() context['{}'] = {}".format(key, context[key]))

        granule = context['granule']
        satellite = context['satellite']

        # Get the location of the binary package
        product = context['product']
        delivery = delivered_software.lookup(
                'fusion_matlab', delivery_id=product.input('fusion_matlab').version)
        dist_root = pjoin(delivery.path, 'dist')
        envroot = pjoin(dist_root, 'env')

        LOG.debug("dist_root = '{}'".format(dist_root))

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

        geo_keys = sorted([key for key in inputs.keys() if 'geo' in key])
        l1b_keys = sorted([key for key in inputs.keys() if 'l1b' in key])
        sounder_keys = sorted([key for key in inputs.keys() if 'sounder' in key])

        LOG.debug("geo_keys = {}".format(geo_keys))
        LOG.debug("l1b_keys = {}".format(l1b_keys))
        LOG.debug("sounder_keys = {}".format(sounder_keys))

        geo = [inputs[key] for key in geo_keys]
        l1b = [inputs[key] for key in l1b_keys]
        sounder = [inputs[key] for key in sounder_keys]

        LOG.debug("geo = \n\t{}".format('\n\t'.join(geo)))
        LOG.debug("l1b = \n\t{}".format('\n\t'.join(l1b)))
        LOG.debug("sounder = \n\t{}".format('\n\t'.join(sounder)))

        # Run viirsmend on the viirs level-1b, and generate the CrIS/VIIRS collocation
        if satellite=='snpp' or satellite=='noaa20':

            for idx, l1b_file in enumerate(l1b):
                if 'bowtie' in basename(l1b_file):
                    # Bias corrected l1b has already been mended, so just copy to working dir
                    l1b[idx] = basename(l1b_file)
                    shutil.copy(l1b_file, l1b[idx])
                else:
                    # Mend bowtie pixels of NASA VIIRS l1b file...
                    l1b[idx] = self.mend_viirs_l1b(product, geo[idx], l1b[idx], dummy=dummy)

            collo = self.cris_viirs_collocation(product, context, geo, sounder, dummy=dummy)

        # Generate the AIRS/MODIS collocation
        if satellite == 'aqua':
            geo = inputs['geo']
            l1b = inputs['l1b']
            sounder_keys = [key for key in inputs.keys() if 'sounder' in key]
            sounder = [inputs[key] for key in sounder_keys]
            collo = self.airs_modis_collocation(product, inputs, dummy=dummy)

        LOG.debug("geo = \n\t{}".format('\n\t'.join(geo)))
        LOG.debug("l1b = \n\t{}".format('\n\t'.join(l1b)))
        LOG.debug("sounder = \n\t{}".format('\n\t'.join(sounder)))
        LOG.debug("collo = \n\t{}".format('\n\t'.join(collo)))

        bin_dir = pjoin(dist_root, 'bin')
        cdl_dir = pjoin(dist_root, 'cdl')
        anc_dir = pjoin(dist_root, 'luts')
        fused_output_dir = pjoin(work_dir, 'fused_outputs')

        # Setup the require keyword arguments for the fusion_matlab package
        kwargs = {}
        kwargs['py_interp'] = py_interp
        kwargs['bin_dir'] = bin_dir
        kwargs['cdl_dir'] = cdl_dir
        kwargs['env'] = env
        kwargs['fused_output_dir'] = fused_output_dir
        kwargs['satellite'] = satellite
        kwargs['granule'] = granule
        kwargs['dummy'] = dummy

        if satellite=='snpp' or satellite=='noaa20':
            geo_file = geo[1]
            l1b_file = l1b[1]
            sounder_file = sounder[1]
            kwargs['anc_paths'] = [pjoin(anc_dir, 'modis_aqua.srf.nc'),
                                   pjoin(anc_dir, 'NG_VIIRS_NPP_RSR_filtered_Oct2011_BA/')]
            kwargs['fusion_binary'] = 'run_imagersounderfusion_V.sh'
            kwargs['matlab_file_glob'] = 'fusion_output.mat'
            kwargs['conversion_bin'] = pjoin(envroot, 'bin', 'l1b-fusion-viirs-cris')
        elif satellite=='aqua':
            geo_file = geo[0]
            l1b_file = l1b[0]
            sounder_file = sounder[0]
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

        if matlab_file is None:
            raise RuntimeError('Output fusion file fusion_output.mat not created.')

        #kwargs['dummy'] = True # dummy
        #kwargs['dummy'] = False # dummy

        # Now that we've computed the Matlab file, convert to a NetCDF file...
        if satellite=='snpp' or satellite=='noaa20':
            rc_fusion, fused_l1b_file, output_attrs = self.convert_matlab_to_viirs_netcdf(
                                                                                   product,
                                                                                   l1b_file,
                                                                                   sounder_file,
                                                                                   matlab_file,
                                                                                   **kwargs)
        #elif satellite=='aqua':
            #rc_fusion, fused_l1b_file, output_attrs = self.convert_matlab_to_netcdf(
                                                                                   #product,
                                                                                   #matlab_file,
                                                                                   #l1b_file,
                                                                                   #**kwargs)

        LOG.debug('convert_matlab_to_netcdf() return value: {}'.format(rc_fusion))
        LOG.info('convert_matlab_to_netcdf() generated {}'.format(fused_l1b_file))

        if fused_l1b_file is None:
            raise RuntimeError('NetCDF output fusion file was not created.')

        # The staging routine assumes that the output file is located in the work directory
        # "tmp******", and that the output path is to be prepended, so return the basename.
        out_fn = basename(fused_l1b_file)

        # Update the output file metadata to conform to LAADS requirements.
        out_dict = self.set_l2_metadata(fused_l1b_file, l1b_file, geo_file, sounder_file, product, context)

        if satellite=='snpp' or satellite=='noaa20':
            out_compress = nc_compress
        if satellite == 'aqua':
            out_compress = hdf_compress

        LOG.debug('We are in {}'.format(os.getcwd()))

        #LOG.debug('Compressing {}'.format(out_fn))
        #return {'fused_l1b': out_compress(out_fn)}

        return out_dict


class FUSION_MATLAB_QL(Computation):

    parameters = ['granule', 'satellite', 'version']

    outputs = ['fused_l1b_ql_band27_asc', 'fused_l1b_ql_band27_desc', 'fused_l1b_ql_band33_asc', 'fused_l1b_ql_band33_desc']

    def find_contexts(self, time_interval, satellite, version):
        '''
        Takes the input time interval, and returns whatever 0Z datetimes fall within the interval.
        '''

        if satellite=='snpp' or satellite=='noaa20':
            granule_length = timedelta(minutes=6)
        elif satellite=='aqua':
            granule_length = timedelta(minutes=5)
        else:
            return []

        LOG.debug('interval = {}'.format(time_interval))
        LOG.debug('satellite = {}'.format(satellite))
        LOG.debug('version = {}'.format(version))

        days = [day.left for day in time_interval.overlapping_interval_series(timedelta(days=1))]

        return [{'satellite': satellite, 'version': version, 'granule': day} for day in days]

    def _add_viirs_l1b_geo_input(self, product, context, task):
        satellite = context['satellite']
        granule = context['granule']
        version = product.input('viirs_l1').version
        input_name = sipsprod.satellite_esdt('V03MOD', satellite)
        interval = TimeInterval(granule, granule+timedelta(days=1.00)-timedelta(seconds=1))
        LOG.debug("Ingesting input {} ({}) for V02FSN_DailyQL version {}".format(input_name, version, product.version))
        vgeom = dawg_catalog.files(satellite, input_name, interval, version=version)
        if vgeom == []:
            raise WorkflowNotReady('Missing {} inputs for version {} and interval {}'.format(
                input_name, version, interval))
        elapsed_time = ((datetime.utcnow()-context['granule']).total_seconds())/86400.
        LOG.info("After {:4.2f} days we have {} {} files ({}).".format(elapsed_time, len(vgeom), input_name, version))
        if len(vgeom) < 228 and elapsed_time < 2.:
            raise WorkflowNotReady('Number of available {} inputs is < 228, for version {} and interval {}, aborting...'.format(
                input_name, version, interval))
        for idx, geo_file in enumerate(vgeom):
            LOG.debug('V03MOD granule {}: {} -> {}'.format(idx, geo_file.begin_time, geo_file.end_time))
            task.input('geo_{}'.format(idx),  geo_file)

    def _add_cris_viirs_fusion_l1b_input(self, product, context, task):
        satellite = context['satellite']
        granule = context['granule']
        version = product.input('FSNRAD_L2_VIIRS_CRIS').version
        input_name = sipsprod.satellite_esdt('FSNRAD_L2_VIIRS_CRIS', satellite)
        interval = TimeInterval(granule, granule+timedelta(days=1.00)-timedelta(seconds=1))
        LOG.debug("Ingesting input {} ({}) for V02FSN_DailyQL version {}".format(input_name, version, product.version))
        vl1b = dawg_catalog.files(satellite, input_name, interval, version=version)
        if vl1b == []:
            raise WorkflowNotReady('Missing {} inputs for version {} and interval {}'.format(
                input_name, version, interval))
        elapsed_time = ((datetime.utcnow()-context['granule']).total_seconds())/86400.
        LOG.info("After {:4.2f} days we have {} {} files ({}).".format(elapsed_time, len(vl1b), input_name, version))
        if len(vl1b) < 228 and elapsed_time < 2.:
            raise WorkflowNotReady('Number of available {} inputs is < 228, for version {} and interval {}, aborting...'.format(
                input_name, version, interval))
        for idx, l1b_file in enumerate(vl1b):
            LOG.debug('FSNRAD_L2_VIIRS_CRIS granule {}: {} -> {}'.format(idx, l1b_file.begin_time, l1b_file.end_time))
            task.input('l1b_{}'.format(idx),  l1b_file)

    def build_task_snpp(self, context, task):
        '''
        Build up a set of inputs for a single context
        '''
        LOG.debug("Ingesting inputs for V02FSN_DailyQL version {} ...".format(context['version']))

        # Get the product definition for 'V02FSN_DailyQL'
        product = sipsprod.lookup_product_recurse('V02FSN_DailyQL', version=context['version'])

        # Ingest the required inputs, defined in the VNP02 product definition for context['version']
        self._add_viirs_l1b_geo_input(product, context, task)
        self._add_cris_viirs_fusion_l1b_input(product, context, task)

        # Make the product definition available to build_task()
        task.option('product', product)

    @reraise_as(WorkflowNotReady, FileNotFound, prefix='FUSION_MATLAB_QL')
    def build_task(self, context, task):

        satellite = context['satellite']

        if satellite=='snpp' or satellite=='noaa20':
            self.build_task_snpp(context, task)
        elif satellite == 'aqua':
            raise ValueError('Satellite option "{}" not yet implemented.'.format(satellite))
        else:
            raise ValueError('Invalid satellite: {}'.format(satellite))

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

    def run_fusion_quicklooks(self, fsn_dir, geo_dir, **kwargs):

        bin_dir = kwargs['bin_dir']
        fusion_ql_binary = kwargs['fusion_ql_binary']
        granule = kwargs['granule']
        satellite = kwargs['satellite']
        env = kwargs['env']

        dt = granule
        year = dt.utctimetuple().tm_year
        jday = dt.utctimetuple().tm_yday

        rc_fusion_ql = 0

        # Get the matlab runtim version that we require
        matlab_version = '2015b'
        #matlab_version = product.input('fusion_matlab').options['matlab_version']

        #run matlab
        #cmd = '{}/{} {} {} {} {}/ {}/  >> fusion_quicklooks.log'.format(
        cmd = '{}/{} {} {} {} {}/ {}/'.format(
            bin_dir,
            fusion_ql_binary,
            support_software.lookup('matlab', matlab_version).path,
            year,
            jday,
            fsn_dir,
            geo_dir
            )

        # Run the Matlab Fusion Quicklook code
        try:
            LOG.debug("cmd = \\\n\t{}".format(cmd.replace(' ',' \\\n\t')))
            rc_fusion_ql = 0
            runscript(cmd, requirements=[], env=env)
        except CalledProcessError as err:
            rc_fusion_ql = err.returncode
            LOG.error("Matlab binary {} returned a value of {}".format(fusion_ql_binary, rc_fusion_ql))
            return rc_fusion_ql, []

        # Move matlab file to the output directory
        orig_fusion_ql_files = glob('*.png')
        if len(orig_fusion_ql_files) != 0:
            LOG.info('Found Fusion quicklook files {}.'.format(
                ', '.join([basename(x) for x in orig_fusion_ql_files])
                ))
        else:
            LOG.error('There are no Fusion quicklook files "*.png", aborting')
            rc_fusion_ql = 1
            return rc_fusion_ql

        fusion_ql_files = []
        for orig_fusion_ql_file in orig_fusion_ql_files:
            filename_chunks = splitext(orig_fusion_ql_file)[0].split('_')
            fusion_ql_files.append(
                    '{}.A{}.{}.{}.png'.format(
                        sipsprod.satellite_esdt('FSNRAD_L2_VIIRS_CRIS', satellite),
                        granule.strftime('%Y%j'),
                        filename_chunks[1],
                        filename_chunks[3]
                        )
                    )

        for oldfile, newfile in zip(orig_fusion_ql_files, fusion_ql_files):
            shutil.move(oldfile, newfile)

        return rc_fusion_ql, fusion_ql_files

    @reraise_as(WorkflowNotReady, FileNotFound, prefix='FUSION_MATLAB_QL')
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

        LOG.debug("dist_root = '{}'".format(dist_root))

        # Get the required  environment variables
        env = self.prepare_env(dist_root, inputs, context)

        # What is the path of the python interpreter
        py_interp = "{}/bin/python".format(envroot)
        LOG.debug("py_interp = '{}'".format(py_interp))

        bin_dir = pjoin(dist_root, 'bin')

        # Where are we running the package
        work_dir = abspath(curdir)
        LOG.debug("working dir = {}".format(work_dir))

        # What are our inputs?
        for input in inputs.keys():
            inputs_dir = dirname(inputs[input])
            LOG.debug("inputs['{}'] = {}".format(input,inputs[input]))
        LOG.debug("Inputs dir = {}".format(inputs_dir))

        current_dir = os.getcwd()
        geo_dir = pjoin(current_dir,'GEO')
        fsn_dir = pjoin(current_dir,'FSN')
        create_dir(geo_dir)
        create_dir(fsn_dir)

        os.chdir(geo_dir)
        geo_keys = [key for key in inputs.keys() if 'geo' in key]
        geo_inputs = {key: inputs[key] for key in geo_keys}
        symlink_inputs_to_working_dir(geo_inputs)
        for geo_file in glob('VJ103*.nc'):
            shutil.move(geo_file, geo_file.replace('VJ103', 'VNP03'))
        os.chdir(current_dir)

        os.chdir(fsn_dir)
        fsn_keys = [key for key in inputs.keys() if 'l1b' in key]
        fsn_inputs = {key: inputs[key] for key in fsn_keys}
        symlink_inputs_to_working_dir(fsn_inputs)
        for geo_file in glob('VJ102*.nc'):
            shutil.move(geo_file, geo_file.replace('VJ102', 'VNP02'))
        os.chdir(current_dir)

        # Setup the require keyword arguments for the fusion_matlab package
        kwargs = {}
        kwargs['py_interp'] = py_interp
        kwargs['bin_dir'] = bin_dir
        kwargs['env'] = env
        kwargs['satellite'] = satellite
        kwargs['granule'] = granule
        kwargs['fusion_ql_binary'] = 'run_plot_globalVIIRSfusion_fct.sh'

        # Run the fusion_matlab package
        rc_fusion_ql, fusion_ql_files = self.run_fusion_quicklooks(
                                                   fsn_dir,
                                                   geo_dir,
                                                   **kwargs
                                                  )

        extra_attrs = {
                       'begin_time': context['granule'],
                       'end_time': context['granule']+timedelta(days=1)-timedelta(seconds=1)
                      }

        LOG.debug('extra_attrs = {}'.format(extra_attrs))

        LOG.debug('run_fusion_quicklooks() return value: {}'.format(rc_fusion_ql))
        LOG.info('run_fusion_quicklooks() generated {}'.format(fusion_ql_files))

        return {
                'fused_l1b_ql_band27_asc': {'file': [x for x in fusion_ql_files if 'Band27.asc'  in x][0], 'extra_attrs': extra_attrs},
                'fused_l1b_ql_band27_desc': {'file': [x for x in fusion_ql_files if 'Band27.desc' in x][0], 'extra_attrs': extra_attrs},
                'fused_l1b_ql_band33_asc': {'file': [x for x in fusion_ql_files if 'Band33.asc' in x][0], 'extra_attrs': extra_attrs},
                'fused_l1b_ql_band33_desc': {'file': [x for x in fusion_ql_files if 'Band33.desc' in x][0], 'extra_attrs': extra_attrs},
                }
