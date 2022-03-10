# Various input file formats that are supported by silico.
# Most input files are actually coordinate files (.xyz, .com etc) and these are all converted to the same internal coordinate representation (the Input_coord class).
# However, sometimes calculations use other types of files as input, for example entire directories (often containing a previous calculation),
# or program specific binary files (.chk is the most obious example).
# For these more complicated scenarios, specific input file classes are provided.
#
# These files might properly belong in silico.submit,
# but because the Silico_coord class in particular is used in various places (it is used for all file interconversions by silico convert for example),
# they are placed more centrally here (for now at least).

from .base import Input_file
from .coord import Silico_coords