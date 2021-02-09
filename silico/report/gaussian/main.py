# General imports.
from pathlib import Path

# Silico imports.
from silico.report.base.pdf import PDF_report
from silico.file.fchk import Chk_to_fchk
from silico.file.cube import Fchk_to_spin_cube, Fchk_to_cube
import silico.file.types as file_types

class Gaussian_report(PDF_report):
    """
    A specialised report object for processing Gaussian results.
    """
    
    # A dictionary of file types accepted by this report.
    INPUT_FILES = {
            file_types.gaussian_chk_file: 'chk_file_path',
            file_types.gaussian_fchk_file: 'fchk_file_path'
        }
    
    def __init__(self, result, *, chk_file_path = None, fchk_file_path = None, options):
        """
        Constructor for Gaussian reports.
        
        One (or both) of either chk_file_path or fchk_file_path must be given in order for orbital/molecular images to be rendered.
        
        :param result: A result_set object.
        :param chk_file_path: Optional path to a chk_file.
        :param fchk_file_path: Optional path to an fchk_file.
        """
        self.chk_file_path = chk_file_path
        self.fchk_file_path = fchk_file_path
        
        # Continue in parent.
        super().__init__(result, options = options)
        
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, get our fchk file (from which cubes are made in Gaussian.
        self.fchk_file = Chk_to_fchk(
            Path(output_dir, output_name + ".fchk"),
            chk_file = self.chk_file_path,
            fchk_file = self.fchk_file_path,
        )
        
        ################
        # Spin density #
        ################
        self.cubes['spin_density'] = Fchk_to_spin_cube.from_options(
            Path(output_dir, "Spin Density", output_name + ".spin.cube"),
            fchk_file = self.fchk_file,
            spin_density = "SCF",
            options = self.options
        )
        
        
        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals):
            for orbital in orbital_list:
                # First, decide what type of orbital we need.
                if orbital.spin_type == "alpha":
                    cubegen_type = "AMO"
                elif orbital.spin_type == "beta":
                    cubegen_type = "BMO"
                else:
                    cubegen_type = "MO"
                
                # Save cube.
                self.cubes[orbital.label] = Fchk_to_cube.from_options(
                    Path(output_dir, orbital.label, output_name + ".{}.cube".format(orbital.label)),
                    fchk_file = self.fchk_file,
                    cubegen_type = cubegen_type,
                    orbital = orbital.level,
                    options = self.options)
        
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO']
        elif "HOMO (alpha)" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO (alpha)']
        else:
            # No MO cubes available, create one for structure.
            # We'll just use the HOMO to get our cube, as it almost certainly should exist.
            self.cubes['structure'] = Fchk_to_cube.from_options(
                Path(output_dir, "Structure", output_name + ".structure.cube"),
                fchk_file = self.fchk_file,
                cubegen_type = "MO",
                orbital = "HOMO",
                options = self.options)
    
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        # Delete our fchk file.
        self.fchk_file.delete(lazy = True)
        
        # Continue.
        super().cleanup()
      
        