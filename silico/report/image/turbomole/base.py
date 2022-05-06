# General imports.
from pathlib import Path

# Silico imports.
from silico.report.image.main import Cube_setup, Partial_cube_setup
from silico.file.cube import Turbomole_to_cube, Turbomole_to_spin_cube,\
    Turbomole_to_density_cube, Turbomole_to_orbital_cube,\
    Turbomole_to_anadens_cube
from silico.file.base import Dummy_file_maker


class Turbomole_setup(Partial_cube_setup):
    """
    Class for setting up Turbomole images.
    """
    
    def __init__(self, report, *, result, do_orbitals, do_spin, options, calculation = None):
        """
        :param report: The report object we will setup images for.
        :param result: The result corresponding to the (sub) calculation we will make images from.
        :param do_orbitals: Whether to generate orbitals.
        :param do_spin: Whether to generate spin density plots.
        :param options: Dictionary of config options.
        :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
        """
        super().__init__(report, result = result, options = options, calculation = calculation)
        
        # Find which turbomole directories we can use to generate cubes.
        self.calculation_directories = {'structure': None, 'spin': None}
        
        try:
            if do_orbitals:
                self.calculation_directories['structure'] = result.metadata.log_files[0].parent
            
            if do_spin:
                self.calculation_directories['spin'] = result.metadata.log_files[0].parent
            
        except (IndexError, AttributeError):
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
                # Cubes can only be created from a real image maker, the dummy class lacks the methods required by the Turbomole_to_orbital_cube init() function.
                #TODO: Don't really like this, ideally the constructor would seamlesly use a Dummy_file_maker
                if not isinstance(self.cube_makers['structure'], Dummy_file_maker):
                    self.report.cubes[orbital.label] = Turbomole_to_orbital_cube(turbomole_to_cube = self.cube_makers['structure'], orbital = orbital)
                
                else:
                    self.report.cubes[orbital.label] = Dummy_file_maker(orbital.label+ " cube file", self.cube_makers['structure'].message)
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.report.cubes:
            self.report.cubes['structure'] = self.report.cubes['HOMO']
        elif "HOMO (alpha)" in self.report.cubes:
            self.report.cubes['structure'] = self.report.cubes['HOMO (alpha)']


class Turbomole_anadens_setup(Cube_setup):
    """
    Class for setting up anadens cubes.
    """
    
    def __init__(self, report, options, calculation = None):
        """
        :param report: The report object we will setup images for.
        :param options: Dictionary of config options.
        :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
        """
        self.report = report
        self.options = options
        self.calculation = calculation

    def setup(self, output_dir, output_name):
        """
        Perform setup.
        
        Calling this method will set cube objects in the parent report.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # We can only do setup if we have a ground state density and at least some excited state densities.
        if 'ground_state_cao_file' in self.report.result.metadata.auxiliary_files and 'excited_state_cao_files' in self.report.result.metadata.auxiliary_files:
            ground_state_cao_file = self.report.result.metadata.auxiliary_files['ground_state_cao_file']
            excited_state_cao_files = self.report.result.metadata.auxiliary_files['excited_state_cao_files']
            
            # Handle each available excited state.
            for excited_state in self.report.result.excited_states:
                if excited_state.state_symbol in excited_state_cao_files:
                    # Got a relevant density file.
                    excited_state_cao_file = excited_state_cao_files[excited_state.state_symbol]
                    
                    # We use a different constructor depending on whether we have a calculation to use as a template.
                    if self.calculation is not None:
                        cube_maker = Turbomole_to_anadens_cube.from_calculation(
                            Path(output_dir, "Cubes", excited_state.state_symbol + "_differential_density.cub"),
                            turbomole_calculation = self.calculation,
                            # TODO: This relies on the density files still being located in the correct turbomole calculation directory, which may not always be the case...
                            calculation_directory = excited_state_cao_file.parent,
                            first_density = ground_state_cao_file,
                            second_density = excited_state_cao_file,
                            options = self.options
                        )
                        
                    else:
                        cube_maker = Turbomole_to_anadens_cube.from_options(
                            Path(output_dir, "Cubes", excited_state.state_symbol + "_differential_density.cub"),
                            # TODO: This relies on the density files still being located in the correct turbomole calculation directory, which may not always be the case...
                            calculation_directory = excited_state_cao_file.parent,
                            first_density = ground_state_cao_file,
                            second_density = excited_state_cao_file,
                            options = self.options
                        )
                        
                    self.report.cubes[excited_state.state_symbol + "_differential_density"] = cube_maker
