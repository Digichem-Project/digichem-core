"""
Classes for parsing calculation result data.

Most of the heavy lifting is done by cclib, we just extract additional data not currently handed by cclib.
"""

# These alignment classes are needed to parse correctly.
from digichem.result.alignment.AAA import Adjusted_average_angle
from digichem.result.alignment.AA import Average_angle
from digichem.result.alignment.FAP import Furthest_atom_pair
from digichem.result.alignment import Minimal

from digichem.parse.util import from_log_files, parse_calculation, parse_and_merge_calculations, parse_multiple_calculations, parse_and_merge_multiple_calculations