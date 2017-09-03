import sys
from datetime import datetime, timedelta
import traceback
import logging

from flo.time import TimeInterval
from flo.ui import local_prepare, local_execute

from flo.sw.fusion_matlab import FUSION_MATLAB
from flo.sw.fusion_matlab.utils import setup_logging

# every module should have a LOG object
LOG = logging.getLogger(__file__)

setup_logging(2)

comp = FUSION_MATLAB()

satellite = 'snpp'
#satellite = 'aqua'
delivery_id = '20170831-1'

# Specify the intervals
wedge = timedelta(seconds=1.)
intervals = [
    # TimeInterval(datetime(2014, 1, 1, 0, 0), datetime(2014, 2, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 2, 1, 0, 0), datetime(2014, 3, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 3, 1, 0, 0), datetime(2014, 4, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 4, 1, 0, 0), datetime(2014, 5, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 5, 1, 0, 0), datetime(2014, 6, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 6, 1, 0, 0), datetime(2014, 7, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 7, 1, 0, 0), datetime(2014, 8, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 8, 1, 0, 0), datetime(2014, 9, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 9, 1, 0, 0), datetime(2014, 10, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 10, 1, 0, 0), datetime(2014, 11, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 11, 1, 0, 0), datetime(2014, 12, 1, 0, 0) - wedge),
    #TimeInterval(datetime(2014, 12, 1, 0, 0), datetime(2015, 1, 1, 0, 0) - wedge)
    # TimeInterval(datetime(2014, 1, 1, 0, 0), datetime(2014, 2, 1, 0, 0) - wedge)
    #TimeInterval(datetime(2016, 5, 17, 5, 30), datetime(2016, 5, 17, 6, 0) - wedge),
    TimeInterval(datetime(2016, 5, 18, 3, 36), datetime(2016, 5, 18, 3, 49) - wedge)
]

LOG.info("Submitting intervals...")
for interval in intervals:
    LOG.info("Submitting interval {} -> {}".format(interval.left, interval.right))
    
    contexts =  comp.find_contexts(interval, satellite, delivery_id)

    LOG.info("\tThere are {} contexts in this interval".format(len(contexts)))
    contexts.sort()
    
    for context in contexts:
        print context

    LOG.info("\tFirst context: {}".format(contexts[0]))
    LOG.info("\tLast context:  {}".format(contexts[-1]))
    #LOG.info("\t{}".format(safe_submit_order(comp, [comp.dataset('fused_l1b')], contexts))
    #time.sleep(30.)
