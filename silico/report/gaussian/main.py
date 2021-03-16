# General imports.
from pathlib import Path

# Silico imports.
from silico.report.main import PDF_report
from silico.file.fchk import Chk_to_fchk
from silico.file.cube import Fchk_to_spin_cube, Fchk_to_cube,\
    Fchk_to_density_cube

class Gaussian_report(PDF_report):
    """
    A specialised report object for processing Gaussian results.
    """
    
    def __init__(self, result, *, options):
        """
        Constructor for Gaussian reports.
        """
        # Get the chk and fchk files we'll be using.
        self.chk_file_paths = {"structure": None, "spin": None}
        self.fchk_file_paths = {"structure": None, "spin": None}
        
        # We want the first available file of each type, so go backwards.
        for metadata in reversed(result.metadatas):
            # First look for general 'structure' files.
            if "chk_file" in metadata.auxiliary_files:
                self.chk_file_paths['structure'] = metadata.auxiliary_files['chk_file']
                
                # If this metadata also has spin data, this will do for our spin chk too.
                if metadata.multiplicity != 1:
                    self.chk_file_paths['spin'] = metadata.auxiliary_files['chk_file']
                
            if "fchk_file" in metadata.auxiliary_files:
                self.fchk_file_paths['structure'] = metadata.auxiliary_files['fchk_file']
                
                # If this metadata also has spin data, this will do for our spin fchk too.
                if metadata.multiplicity != 1:
                    self.fchk_file_paths['spin'] = metadata.auxiliary_files['fchk_file']
        
        # Actual fchk file maker objects.
        # These cannot be set here as they depend on our output dir.
        self.fchk_files = {}
        
        # Continue in parent.
        super().__init__(result, options = options)
        
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, get our fchk files (from which cubes are made in Gaussian.
        self.fchk_files = {
            "structure": Chk_to_fchk(
                Path(output_dir, output_name + ".fchk"),
                chk_file = self.chk_file_paths['structure'],
                fchk_file = self.fchk_file_paths['structure'],
            )
        }
        
        self.fchk_files = {
            "spin": Chk_to_fchk(
                Path(output_dir, output_name + ".spin.fchk"),
                chk_file = self.chk_file_paths['spin'],
                fchk_file = self.fchk_file_paths['spin'],
            )
        }

        
        ################
        # Spin density #
        ################
        self.cubes['spin_density'] = Fchk_to_spin_cube.from_options(
            Path(output_dir, "Spin Density", output_name + ".spin.cube"),
            fchk_file = self.fchk_files['spin'],
            spin_density = "SCF",
            options = self.options
        )
        
        #################
        # Total density #
        #################
        self.cubes['SCF'] = Fchk_to_density_cube.from_options(
            Path(output_dir, "Density", output_name + ".SCF.cube"),
            fchk_file = self.fchk_files['structure'],
            density_type = "SCF",
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
                    fchk_file = self.fchk_files['structure'],
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
                fchk_file = self.fchk_files['structure'],
                cubegen_type = "MO",
                orbital = "HOMO",
                options = self.options)
    
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        # Delete our fchk files.
        for fchk_file in self.fchk_files:
            self.fchk_files[fchk_file].delete(lazy = True)
        
        # Continue.
        super().cleanup()
      
        