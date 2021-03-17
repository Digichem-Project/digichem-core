# General imports.
from pathlib import Path

# Silico imports.
from silico.report.image.main import Image_setup
from silico.file.cube import Turbomole_to_cube, Turbomole_to_spin_cube,\
    Turbomole_to_density_cube, Turbomole_to_orbital_cube


class Turbomole_setup(Image_setup):
    """
    Class for setting up Turbomole images.
    """
    
    def __init__(self, report, metadata, options, orbitals, spin, calculation = None):
        """
        """
        super().__init__(report, metadata, options, calculation = calculation)
        
        # Find which turbomole directories we can use to generate cubes.
        self.calculation_directories = {'structure': None, 'spin': None}
        
        try:
            self.calculation_directories['structure'] = metadata.log_files[0].parent
            
            # See if this will also do for our spin calc.
            if metadata.multiplicity != 1 and metadata.orbital_spin_type == "unrestricted":
                self.calculation_directories['spin'] = metadata.log_files[0].parent
            
        except (IndexError, KeyError):
            # Either no log files or empty list of log files.
            pass
            
        self.cube_makers = {}
        
    def setup(self, output_dir, output_name):
        """
        Perform setup.
        
        Calling this method will set cube objects in the parent report.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        self.setup_cubes(output_dir, output_name)
    
    
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, we need to decide which cubes we're making.
        required_orbitals = self.report.orbitals_to_render
                
        # Now get our main cube maker.
        if self.calculation is not None:
            # Normal structure maker.
            self.cube_makers['structure'] = Turbomole_to_cube.from_calculation(
                Path(output_dir, "Cubes"),
                turbomole_calculation = self.calculation,
                calculation_directory = self.calculation_directories['structure'],
                orbitals = required_orbitals,
                density = True,
                # If we can use the same directory for both structure and spin, set spin to yes please.
                spin = self.calculation_directories['spin'] is not None and self.calculation_directories['spin'] == self.calculation_directories['structure'],
                options = self.options
            )
            
            if self.calculation_directories["spin"] != None:
                # If we can use the same directory for both structure and spin, use the same cube maker.
                if self.calculation_directories['spin'] == self.calculation_directories['structure']:
                    self.cube_makers['spin'] = self.cube_makers['structure']
                else:
                    # We need to use a different directory for spin.
                    self.cube_makers['spin'] = Turbomole_to_cube.from_calculation(
                        Path(output_dir, "Cubes"),
                        turbomole_calculation = self.calculation,
                        calculation_directory = self.calculation_directories['spin'],
                        orbitals = [],
                        density = True,
                        spin = True,
                        options = self.options
                    )
        else:
            # Normal structure maker.
            self.cube_makers['structure'] = Turbomole_to_cube.from_options(
                Path(output_dir, "Cubes"),
                calculation_directory = self.calculation_directories['structure'],
                orbitals = required_orbitals,
                density = True,
                # If we can use the same directory for both structure and spin, set spin to yes please.
                spin = self.calculation_directories['spin'] is not None and self.calculation_directories['spin'] == self.calculation_directories['structure'],
                options = self.options
            )
            
            if self.calculation_directories["spin"] is not None:
                # If we can use the same directory for both structure and spin, use the same cube maker.
                if self.calculation_directories['spin'] == self.calculation_directories['structure']:
                    self.cube_makers['spin'] = self.cube_makers['structure']
                    print("Using same")
                else:
                    # We need to use a different directory for spin.
                    self.cube_makers['spin'] = Turbomole_to_cube.from_options(
                        Path(output_dir, "Cubes"),
                        calculation_directory = self.calculation_directories['spin'],
                        orbitals = [],
                        density = True,
                        spin = True,
                        options = self.options
                    )
        
        ################
        # Spin density #
        ################
        if "spin" in self.cube_makers:
            self.report.cubes['spin_density'] = Turbomole_to_spin_cube(turbomole_to_cube = self.cube_makers['spin'])
        
        #################
        # Total density #
        #################
        self.report.cubes['SCF'] = Turbomole_to_density_cube(turbomole_to_cube = self.cube_makers['structure'])
        
        
        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.report.result.molecular_orbitals, self.report.result.beta_orbitals):
            for orbital in orbital_list:                
                # Save cube.
                self.report.cubes[orbital.label] = Turbomole_to_orbital_cube(turbomole_to_cube = self.cube_makers['structure'], orbital = orbital)    
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.report.cubes:
            self.report.cubes['structure'] = self.report.cubes['HOMO']
        elif "HOMO (alpha)" in self.report.cubes:
            self.report.cubes['structure'] = self.report.cubes['HOMO (alpha)']