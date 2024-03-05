# Various input file formats that are supported by digichem.
# Most input files are actually coordinate files (.xyz, .com etc) and these are all converted to the same internal coordinate representation (the Input_coord class).
# However, sometimes calculations use other types of files as input, for example entire directories (often containing a previous calculation),
# or program specific binary files (.chk is the most obious example).
# For these more complicated scenarios, specific input file classes are provided.
#

from .base import Input_file
from digichem.input.digichem import Digichem_coords, si_from_file
from .directory import Calculation_directory_input
from .chk import Chk_input