import sys
import traceback
import logging

from timeutil import TimeInterval, datetime, timedelta
from flo.ui import local_prepare, local_execute

from flo.sw.fusion_matlab import FUSION_MATLAB
from flo.sw.fusion_matlab.utils import setup_logging

# every module should have a LOG object
LOG = logging.getLogger(__name__)

comp = FUSION_MATLAB()

#
# Local execution
#

# General information
#granule = datetime(2015, 4, 17, 14)
#interval = TimeInterval(granule, granule+timedelta(minutes=0))

satellite = 'snpp'
#satellite = 'aqua'
version = '1.0dev4' # base VIIRS l1vel-1b
version = '1.0dev5' # bias corrected VIIRS level-1b


def local_execute_example(interval, satellite, version, skip_prepare=False, skip_execute=False, verbosity=2):

    setup_logging(verbosity)

    if satellite == 'snpp':
        LOG.info("We are doing NPP...")
    elif satellite == 'aqua':
        LOG.info("We are doing AQUA...")
    else:
        LOG.error("Invalid satellite.")

    # Get the required context...
    contexts =  comp.find_contexts(interval, satellite, version)

    if len(contexts) != 0:
        LOG.info("Candidate contexts in interval...")
        for context in contexts:
            print("\t{}".format(context))

        try:
            if not skip_prepare:
                LOG.info("Running fusion_matlab local_prepare()...")
                LOG.info("Preparing context... {}".format(contexts[0]))
                local_prepare(comp, contexts[0])
            if not skip_execute:
                LOG.info("Running local_execute()...")
                LOG.info("Running context... {}".format(contexts[0]))
                local_execute(comp, contexts[0])
        except Exception, err:
            LOG.error("{}".format(err))
            LOG.debug(traceback.format_exc())
    else:
        LOG.error("There are no valid {} contexts for the interval {}.".format(satellite, interval))

def print_contexts(interval, satellite, type, version):
    contexts = comp.find_contexts(interval, satellite, version)
    for context in contexts:
        print context
