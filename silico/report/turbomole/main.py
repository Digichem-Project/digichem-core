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
    
    def __init__(self, result, *, turbomole_calculation = None, options):
        """
        Constructor for Turbomole report generator.
        
        :param: turbomole_program: An optional turbomole calculation which will be used as a  template to generate cubes.
        """
        super().__init__(result, options = options)
        self.turbomole_calculation = turbomole_calculation
        
        # Find which turbomole directories we can use to generate cubes.
        self.calculation_directories = {'structure': None, 'spin': None}
        
        for metadata in reversed(result.metadatas):
            try:
                self.calculation_directories['structure'] = metadata.log_files[0].parent
                
                # See if this will also do for our spin calc.
                if metadata.multiplicity != 1:
                    self.calculation_directories['spin'] = metadata.log_files[0].parent
                
            except (IndexError, KeyError):
                # Either no log files or empty list of log files.
                pass
            
        self.cube_makers = {}
    
    
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
            # Normal structure maker.
            self.cube_makers['structure'] = Turbomole_to_cube.from_calculation(
                Path(output_dir, "Cubes"),
                turbomole_calculation = self.turbomole_calculation,
                calculation_directory = self.calculation_directories['structure'],
                orbitals = required_orbitals,
                density = True,
                spin = self.result.metadata.multiplicity != 1,
                options = self.options
            )
            
            # And spin.
            self.cube_makers['spin'] = Turbomole_to_cube.from_calculation(
                Path(output_dir, "Cubes"),
                turbomole_calculation = self.turbomole_calculation,
                calculation_directory = self.calculation_directories['spin'],
                orbitals = required_orbitals,
                density = True,
                spin = self.result.metadata.multiplicity != 1,
                options = self.options
            )
        else:
            # Normal structure maker.
            self.cube_makers['structure'] = Turbomole_to_cube.from_options(
                Path(output_dir, "Cubes"),
                calculation_directory = self.calculation_directories['structure'],
                orbitals = required_orbitals,
                density = True,
                spin = self.result.metadata.multiplicity != 1,
                options = self.options
            )
            
            # And spin.
            self.cube_makers['spin'] = Turbomole_to_cube.from_options(
                Path(output_dir, "Cubes"),
                calculation_directory = self.calculation_directories['spin'],
                orbitals = required_orbitals,
                density = True,
                spin = self.result.metadata.multiplicity != 1,
                options = self.options
            )
        
        ################
        # Spin density #
        ################
        self.cubes['spin_density'] = Turbomole_to_spin_cube(turbomole_to_cube = self.cube_makers['spin'])
        
        #################
        # Total density #
        #################
        self.cubes['SCF'] = Turbomole_to_density_cube(turbomole_to_cube = self.cube_makers['structure'])
        
        
        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals):
            for orbital in orbital_list:                
                # Save cube.
                self.cubes[orbital.label] = Turbomole_to_orbital_cube(turbomole_to_cube = self.cube_makers['structure'], orbital = orbital)    
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO']
        elif "HOMO (alpha)" in self.cubes:
            self.cubes['structure'] = self.cubes['HOMO (alpha)']
            
