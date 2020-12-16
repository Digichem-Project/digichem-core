from silico.parser.base import Sub_parser
import scipy.constants


class Turbomole_orbitals(Sub_parser):
    """
    Sub parser for turbomole orbitals.
    """
       
    HOMO_LUMO_HEADER = "HOMO-LUMO Separation"
        
    def parse(self, current_line, file):
        """
        """
        if self.HOMO_LUMO_HEADER in current_line:
            # The HOMO LUMO section.
            # The next line contains the HOMO energy.
            # Both Hartrees and Ev are printed, but rounding errors make eV unnusuable.
            # Take Hartree and convert ourselves.
            # HOMO         :   -0.25753158 H =     -7.00779 eV
            homo_line = next(file)
            homo = float(homo_line.split()[2]) / scipy.constants.physical_constants['electron volt-hartree relationship'][0]
            
            # Next line is the LUMO energy.
            lumo_line = next(file)
            lumo = float(lumo_line.split()[-2])