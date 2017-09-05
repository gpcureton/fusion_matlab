import sys
from datetime import datetime, timedelta
import traceback
import logging

from flo.time import TimeInterval
from flo.ui import safe_submit_order

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
granule = datetime(2015, 4, 17, 0, 0)
wedge = timedelta(seconds=1.)
intervals = [
    #TimeInterval(granule, granule + timedelta(hours=1) - wedge)
    TimeInterval(granule, granule + timedelta(days=1) - wedge)
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
    LOG.info("\t{}".format(safe_submit_order(comp, [comp.dataset('fused_l1b')], contexts)))
    #time.sleep(30.)
