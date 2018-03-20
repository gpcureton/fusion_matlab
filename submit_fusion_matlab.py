#!/usr/bin/env python
# encoding: utf-8

import sys
import traceback
import calendar
import logging
from calendar import monthrange
from time import sleep

from flo.ui import safe_submit_order
from timeutil import TimeInterval, datetime, timedelta

from flo.sw.fusion_matlab import FUSION_MATLAB
from flo.sw.fusion_matlab.utils import setup_logging

# every module should have a LOG object
LOG = logging.getLogger(__name__)

setup_logging(2)

#
# General information
#
delivery_id = '20180225-1'
version = '1.0dev2'
wedge = timedelta(seconds=1.)
day = timedelta(days=1.)

#
# Satellite specific information
#
#satellite = 'aqua'
#granule_length = timedelta(minutes=5)
satellite = 'snpp'
granule_length = timedelta(minutes=6)

#
# Specify the intervals
#
#granule = datetime(2015, 4, 17, 17, 55) # Aqua
granule = datetime(2018, 2, 2, 18, 36) # SNPP
#intervals = [
    #TimeInterval(granule, granule + granule_length - wedge)
    #TimeInterval(granule, granule + timedelta(hours=1) - wedge)
    #TimeInterval(granule, granule + timedelta(days=1) - wedge)
    #TimeInterval(granule - timedelta(hours=1), granule + timedelta(hours=1))
    #TimeInterval(datetime(2014, 7, 1, 0, 5), datetime(2014, 8, 1) - wedge)
    #TimeInterval(datetime(2014, 7, 1, 0, 0), datetime(2014, 7, 1, 0, 10) - wedge)
    #TimeInterval(datetime(2015, 4, 1, 0, 0), datetime(2015, 5, 1, 0, 0) - wedge)
    #TimeInterval(datetime(2015, 4, 25, 13,  0), datetime(2015, 4, 25,  14,  0) - wedge)
#]
intervals = []
years = 2018
intervals += [TimeInterval(datetime(years,month,1), datetime(years,month,calendar.monthrange(years,month)[1])+day-wedge) for month in range(1,3) ]

#
# Initialize the computation
#
comp = FUSION_MATLAB()

#
# Submit the jobs
#
LOG.info("Submitting intervals...")

dt = datetime.utcnow()
log_name = '/home/flo/geoffc/fusion_matlab_logs/fusion_matlab_{}_s{}_e{}_c{}.log'.format(
    satellite,
    intervals[0].left.strftime('%Y%m%d%H%M'),
    intervals[-1].right.strftime('%Y%m%d%H%M'),
    dt.strftime('%Y%m%d%H%M%S'))

try:
    for interval in intervals:
        LOG.info("Submitting interval {} -> {}".format(interval.left, interval.right))

        contexts =  comp.find_contexts(interval, satellite, version)

        LOG.info("Opening log file {}".format(log_name))
        file_obj = open(log_name,'a')

        LOG.info("\tThere are {} contexts in this interval".format(len(contexts)))
        contexts.sort()

        if contexts != []:
            #for context in contexts:
                #LOG.info(context)

            LOG.info("\tFirst context: {}".format(contexts[0]))
            LOG.info("\tLast context:  {}".format(contexts[-1]))

            try:
                job_nums = []
                job_nums = safe_submit_order(comp, [comp.dataset('fused_l1b')], contexts)

                if job_nums != []:
                    #job_nums = range(len(contexts))
                    #LOG.info("\t{}".format(job_nums))

                    file_obj.write("contexts: [{}, {}]; job numbers: [{}..{}]\n".format(contexts[0], contexts[-1], job_nums[0],job_nums[-1]))
                    LOG.info("contexts: [{}, {}]; job numbers: [{},{}]".format(contexts[0], contexts[-1], job_nums[0],job_nums[-1]))
                    LOG.info("job numbers: [{}..{}]\n".format(job_nums[0],job_nums[-1]))
                else:
                    LOG.info("contexts: [{}, {}]; --> no jobs\n".format(contexts[0], contexts[-1]))
                    file_obj.write("contexts: [{}, {}]; --> no jobs\n".format(contexts[0], contexts[-1]))
            except Exception:
                LOG.warning(traceback.format_exc())

            #sleep(30.)

        LOG.info("Closing log file {}".format(log_name))
        file_obj.close()

except Exception:
    LOG.warning(traceback.format_exc())
