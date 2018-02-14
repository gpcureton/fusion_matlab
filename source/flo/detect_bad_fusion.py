#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Detect obviously-bad fusion data
"""
import os
import sys
import numpy as np
from pyhdf.SD import SD, SDC
import netCDF4 as nc4
#from scipy.io import loadmat
import logging
from itertools import chain

LOG = logging.getLogger(__name__)

DEFAULT_BAND = 32
DEFAULT_THRESHOLD = 0.2


def detect(threshold, tru, tst):
    """Compute RMS error and return pass/fail boolean"""
    yisq = (tst.ravel() - tru.ravel()) ** 2.0
    rms = np.sqrt(np.sum(yisq) / float(len(yisq)))
    LOG.debug("RMS error is {} (threshold {})".format(rms, threshold))
    return bool(rms <= threshold), rms


#def mat_band(filename, band_number=DEFAULT_BAND):
    #"""Yield a single F from a MAT file output from fusion science code"""
    #LOG.info("reading band {} from {} Fusion MAT".format(band_number, filename))
    #mat = loadmat(filename)
    #bandex = np.argwhere(mat['imager_bands'].squeeze() == band_number)[0][0]
    #rad = mat['fusion_rad'][bandex].squeeze()
    #yield rad


def modis_band_nc(filename, band_number=DEFAULT_BAND):
    """Yield M or F from a single MODIS or MODIS-AIRS fusion file, NetCDF4"""
    LOG.info("reading band {} from {} MODIS HDF".format(band_number, filename))
    nc = nc4.Dataset(filename)
    em = nc.variables['EV_1KM_Emissive']
    bn,sc,of = dict(((b,i) for (i,b) in enumerate(map(int, em.band_names.split(','))))), em.radiance_scales, em.radiance_offsets
    band = lambda v, b: sc[bn[b]] * (v[bn[b]] - of[bn[b]])
    yield band(em, band_number)

def modis_band(filename, band_number=DEFAULT_BAND):
    """Yield M or F from a single MODIS or MODIS-AIRS fusion file, using HDF4"""
    LOG.info("reading band {} from {} MODIS HDF".format(band_number, filename))
    print(filename)
    file_obj = SD(str(filename))
    em = file_obj.select('EV_1KM_Emissive')
    bn,sc,of = dict(((b,i) for (i,b) in enumerate(map(int, em.band_names.split(','))))), em.radiance_scales, em.radiance_offsets
    band = lambda v, b: sc[bn[b]] * (v[bn[b]] - of[bn[b]])
    yield band(em, band_number)

def modis_band_new(filename, band_number=DEFAULT_BAND):
    """Yield M or F from a single MODIS or MODIS-AIRS fusion file, using HDF4"""
    LOG.info("reading band {} from {} MODIS HDF".format(band_number, filename))
    print(filename)
    file_obj = SD(str(filename))
    em = file_obj.select('EV_1KM_Emissive')
    bn,sc,of = dict(((b,i) for (i,b) in enumerate(map(int, em.band_names.split(','))))), em.radiance_scales, em.radiance_offsets
    band = lambda v, b: sc[bn[b]] * (v[bn[b]] - of[bn[b]])
    yield band(em, band_number)


def viirs_band(filename, band_number=DEFAULT_BAND):
    """Yield both Mband and Fband from a single file"""
    LOG.info("reading band {} (M and Fusion) from {} VIIRS NetCDF".format(band_number, filename))
    nc = nc4.Dataset(filename)
    ncob = nc['observation_data']
    m = ncob['M%02d' % band_number][:]
    yield m
    f = ncob['Fusion%02d' % band_number][:]
    yield f

def viirs_band_new(filename, band_number=DEFAULT_BAND):
    """Yield both Mband and Fband from a single file"""
    LOG.info("reading band {} (M and Fusion) from {} VIIRS NetCDF".format(band_number, filename))
    nc = nc4.Dataset(filename)
    ncob = nc['observation_data']
    m = ncob['M%02d' % band_number][:]
    f = ncob['Fusion%02d' % band_number][:]
    yield {'M{}'.format(band_number):[m,f]}

def _debug(type_, value, tb):
    """enable with sys.excepthook = debug"""
    if not sys.stdin.isatty():
        sys.__excepthook__(type_, value, tb)
    else:
        import traceback, pdb
        traceback.print_exception(type_, value, tb)
        # …then start the debugger in post-mortem mode.
        pdb.post_mortem(tb)  # more “modern”


#SFX = {'.mat': mat_band, '.hdf': modis_band, '.nc': viirs_band}
SFX = {'.hdf': modis_band, '.nc': viirs_band}
SFX_new = {'.hdf': modis_band, '.nc': viirs_band_new}

def file_QC(input_rms, band, input_files):
    generators = list(SFX[os.path.splitext(p)[-1].lower()](p, band) for p in input_files)
    stuff = list(chain(*generators))
    if len(stuff) != 2:
        raise AssertionError("requires original + fusion input")
    passfail, rms = detect(input_rms, *stuff)
    print("> pass (rms {} < threshold {})".format(rms, input_rms) if passfail else "> FAIL (rms {} > threshold {})".format(
        rms, input_rms))
    return 0 if passfail else 1


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="",
        epilog="",
        fromfile_prefix_chars='@')
    parser.add_argument('-v', '--verbose', dest='verbosity', action="count", default=0,
                        help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true',
                        help="enable interactive PDB debugger on exception")
    parser.add_argument('-B', '--band', dest='band', type=int, default=DEFAULT_BAND,
                        help="choose band number")
    parser.add_argument('-R', '--rms', dest='rms', type=float, default=DEFAULT_THRESHOLD,
                        help="maximum tolerable RMS error")
    parser.add_argument('inputs', nargs='*',
                        help="input files to process (MAT or HDF), requires 2 files")
    args = parser.parse_args()
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=levels[min(3, args.verbosity)])
    if args.debug:
        sys.excepthook = _debug
    # if not args.inputs:
    #     unittest.main()
    #     return 0

    # Creates file reading generators for each file type in the inputs.
    generators = list(SFX[os.path.splitext(p)[-1].lower()](p, args.band) for p in args.inputs)
    stuff = list(chain(*generators))
    if len(stuff) != 2:
        raise AssertionError("requires original + fusion input")
    passfail, rms = detect(args.rms, *stuff)
    print("> pass (rms {} < threshold {})".format(rms, args.rms) if passfail else "> FAIL (rms {} > threshold {})".format(
        rms, args.rms))
    return 0 if passfail else 1


if __name__ == '__main__':
    sys.exit(main())
