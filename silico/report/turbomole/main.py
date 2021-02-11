# General imports
from pathlib import Path

# Silico imports.
from silico.report.base.pdf import PDF_report
from silico.file.cube import Turbomole_to_cube, Turbomole_to_orbital_cube


class Turbomole_report(PDF_report):
    """
    A specialised report object for processing Turbomole results.
    """
    
    def __init__(self, *args, log_file_path, **kwargs):
        super().__init__(*args, **kwargs)
        self.cube_maker = None
        log_file_path = Path(log_file_path)
        self.calculation_directory = log_file_path.parent
        
    
    
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, we need to decide which cubes we're making.
        required_orbitals = [orbital for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals) for orbital in orbital_list]
                
        # Now get our main cube maker.
        self.cube_maker = Turbomole_to_cube.from_options(
            Path(output_dir, "Cubes"),
            calculation_directory = self.calculation_directory,
            orbitals = required_orbitals,
            restricted = self.result.molecular_orbitals.spin_type == "none",
            options = self.options
        )
        
        ################
        # Spin density #
        ################
        # Spin density is not yet supported for Turbomole.
        self.cubes['spin_density'] = None
        
        
        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals):
            for orbital in orbital_list:                
                # Save cube.
                self.cubes[orbital.label] = Turbomole_to_orbital_cube(turbomole_to_cube = self.cube_maker, irrep = orbital.irrep, spin = orbital.spin_type)        
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO']
        elif "HOMO (alpha)" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO (alpha)']
            
