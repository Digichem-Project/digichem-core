# Classes for parsing calculation result data.
#
# Most of the heavy lifting is done by cclib, we just extract additional data not currently handed by cclib.
from .base import from_log_files, parse_calculation, parse_calculations