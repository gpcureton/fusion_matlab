import sys
import traceback
import calendar
import logging

from flo.ui import safe_submit_order
from timeutil import TimeInterval, datetime, timedelta

from flo.sw.fusion_matlab import FUSION_MATLAB
from flo.sw.fusion_matlab.utils import setup_logging

# every module should have a LOG object
LOG = logging.getLogger(__name__)

setup_logging(2)

comp = FUSION_MATLAB()

#satellite = 'aqua'
#granule_length = timedelta(minutes=5)
satellite = 'snpp'
granule_length = timedelta(minutes=6)

version = '1.0dev0'
delivery_id = '20170920-1'
version = '1.0dev0'

# Specify the intervals
#granule = datetime(2014, 7, 6, 2, 0) # Gala Wind
granule = datetime(2015, 4, 17) # Bryan Baum
#granule = datetime(2015, 4, 17, 17, 55) # Bryan Baum
wedge = timedelta(seconds=1.)
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

intervals = [TimeInterval(datetime(2015, 4, day), datetime(2015, 4, day, 23, 59)) for day in range(1, calendar.monthrange(2015, 4)[1]+1)]


LOG.info("Submitting intervals...")
for interval in intervals:
    LOG.info("Submitting interval {} -> {}".format(interval.left, interval.right))

    contexts =  comp.find_contexts(interval, satellite, version)

    LOG.info("\tThere are {} contexts in this interval".format(len(contexts)))
    contexts.sort()

    #for context in contexts:
        #print context

    LOG.info("\tFirst context: {}".format(contexts[0]))
    LOG.info("\tLast context:  {}".format(contexts[-1]))
    LOG.info("\t{}".format(safe_submit_order(comp, [comp.dataset('fused_l1b')], contexts)))
    #time.sleep(30.)
