from silico.parser.base import Sub_parser
import scipy.constants
from silico.misc.base import is_int
import numpy


class Turbomole_orbitals(Sub_parser):
    """
    Sub parser for turbomole orbitals.
    """
       
    HOMO_LUMO_HEADER = "HOMO-LUMO Separation"
    
    def __init__(self, *args, **kwargs):
        """
        """
        super().__init__(*args, **kwargs)
        
        self.homo = None
        self.alpha_homo = None
        self.beta_homo = None
        
    def parse(self, current_line, file):
        """
        """
        # This strategy is difficult to implement because of rounding errors...
#         if self.HOMO_LUMO_HEADER in current_line:
#             # The HOMO LUMO section.
#             # The next line contains the HOMO energy.
#             # Both Hartrees and Ev are printed, but rounding errors make eV unusable.
#             # Take Hartree and convert ourselves.
#             # HOMO         :   -0.25753158 H =     -7.00779 eV
#             homo_line = next(file)
#             homo = float(homo_line.split()[2]) / scipy.constants.physical_constants['electron volt-hartree relationship'][0]
#             
#             # Next line is the LUMO energy.
#             lumo_line = next(file)
#             lumo = float(lumo_line.split()[-2])

        # Parse number of occupied orbitals
        if current_line == "$closed shells\n":
            # Next line.
            homo_line = next(file)
            # This line has the following format:
            # The first symbol is the symmetry.
            # The second identifies orbitals
            # The third is the occupancy (this seems redundant as the occupancy can be inferred from the section name; $closed shells or $alpha/$beta shells).
            #
            # We are interested in the second line which tells us which orbitals are occupied.
            #  a       1-21                                   ( 2 )
            self.homo = self.homo_index_from_line(homo_line)
            
        elif current_line == "$alpha shells\n":
            # Next line.
            homo_line = next(file)
            self.alpha_homo = self.homo_index_from_line(homo_line)
        
        elif current_line == "$beta shells\n":
            # Next line.
            homo_line = next(file)
            self.beta_homo = self.homo_index_from_line(homo_line)
            
    def finalize(self):
        """
        """
        if self.homo is not None:
            # Restricted.
            self.parser.data.homos = numpy.array([self.homo])
        elif self.alpha_homo is not None and self.beta_homo is not None:
            # Unrestricted.
            self.parser.data.homos = numpy.array([self.alpha_homo, self.beta_homo])
        
        
    @classmethod
    def homo_index_from_line(self, homo_line):
        """
        Parse the highest occupied molecular orbital from the relevant line in the control file. 
        """
        # First, get the second part which tells us about occupation.
        orbital_nums = homo_line.split()[1]
        
        # We want the last number from this list.
        homo_num = ""
        
        for char in reversed(orbital_nums):
            # Give up if the character is not an int.
            if not is_int(char):
                break
            
            # Add our character (but in reverse).
            homo_num = char + homo_num
            
        # All done. Return the index (n -1).
        return int(homo_num) -1
        