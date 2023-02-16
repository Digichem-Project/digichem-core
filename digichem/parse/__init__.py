"""
Classes for parsing calculation result data.

Most of the heavy lifting is done by cclib, we just extract additional data not currently handed by cclib.
"""

from silico.parse.util import from_log_files, parse_calculation, parse_and_merge_calculations, parse_multiple_calculations, parse_and_merge_multiple_calculations