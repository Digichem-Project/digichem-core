# General imports
from pathlib import Path

# Silico imports.
from silico.report.main import PDF_report
from silico.file.cube import Turbomole_to_cube, Turbomole_to_orbital_cube,\
    Turbomole_to_spin_cube, Turbomole_to_density_cube


class Turbomole_report(PDF_report):
    """
    A specialised report object for processing Turbomole results.
    """
    
    def __init__(self, *args, log_file_path, turbomole_calculation = None, **kwargs):
        """
        Constructor for Turbomole report generator.
        
        :param: turbomole_program: An optional turbomole calculation which will be used as a  template to generate cubes.
        """
        super().__init__(*args, **kwargs)
        self.cube_maker = None
        log_file_path = Path(log_file_path)
        self.calculation_directory = log_file_path.parent
        self.turbomole_calculation = turbomole_calculation
    
    
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, we need to decide which cubes we're making.
        #required_orbitals = [orbital for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals) for orbital in orbital_list]
        required_orbitals = self.orbitals_to_render
                
        # Now get our main cube maker.
        if self.turbomole_calculation is not None:
            self.cube_maker = Turbomole_to_cube.from_calculation(
                Path(output_dir, "Cubes"),
                turbomole_calculation = self.turbomole_calculation,
                orbitals = required_orbitals,
                density = True,
                spin = self.result.metadata.system_multiplicity != 1,
                options = self.options
            )
        else:
            self.cube_maker = Turbomole_to_cube.from_options(
                Path(output_dir, "Cubes"),
                calculation_directory = self.calculation_directory,
                orbitals = required_orbitals,
                density = True,
                spin = self.result.metadata.system_multiplicity != 1,
                options = self.options
            )
        
        ################
        # Spin density #
        ################
        self.cubes['spin_density'] = Turbomole_to_spin_cube(turbomole_to_cube = self.cube_maker)
        
        #################
        # Total density #
        #################
        self.cubes['SCF'] = Turbomole_to_density_cube(turbomole_to_cube = self.cube_maker)
        
        
        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals):
            for orbital in orbital_list:                
                # Save cube.
                self.cubes[orbital.label] = Turbomole_to_orbital_cube(turbomole_to_cube = self.cube_maker, orbital = orbital)    
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO']
        elif "HOMO (alpha)" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO (alpha)']
            
