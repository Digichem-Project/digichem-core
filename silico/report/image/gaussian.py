# General imports.
from pathlib import Path

# Silico imports
from silico.report.image import Partial_cube_setup
from silico.file.fchk import Chk_to_fchk
from silico.file.cube import Fchk_to_spin_cube, Fchk_to_density_cube,\
    Fchk_to_cube
from silico.file.NTO import Chk_to_NTO_chk


class Gaussian_setup(Partial_cube_setup):
    """
    Class for setting up Gaussian images.
    """
    
    def __init__(self, report, *, result, do_orbitals, do_spin, options, calculation = None):
        """
        :param report: The report object we will setup images for.
        :param result: The result corresponding to the (sub) calculation we will make images from.
        :param do_orbitals: Whether to generate orbitals.
        :param do_spin: Whether to generate spin density plots.
        :param do_NTOs: Whether to generate natural transition orbitals.
        :param options: Dictionary of config options.
        :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
        """
        super().__init__(report, result = result, options = options, calculation = calculation)
        
        # The memory we'll use for formchk and cubegen.
        self.memory = calculation.memory if calculation is not None else None
        
        # Get the chk and fchk files we'll be using.        
        self.chk_file_paths = {"structure": None, "spin": None}
        self.fchk_file_paths = {"structure": None, "spin": None}
        
        # First look for general 'structure' files.
        if "chk_file" in result.metadata.auxiliary_files:
            if do_orbitals:
                self.chk_file_paths['structure'] = result.metadata.auxiliary_files['chk_file']
            
            if do_spin:
                self.chk_file_paths['spin'] = result.metadata.auxiliary_files['chk_file']
                
            if True:
                # NTOs can only be generated from virgin .chk files (chk files converted back from fchk files will not work; the necessary information is apparently lost).
                self.chk_file_paths['NTO'] = result.metadata.auxiliary_files['chk_file']
            
        if "fchk_file" in result.metadata.auxiliary_files:
            if do_orbitals:
                self.fchk_file_paths['structure'] = result.metadata.auxiliary_files['fchk_file']
            
            if do_spin:
                self.fchk_file_paths['spin'] = result.metadata.auxiliary_files['fchk_file']
                
        # Actual chk file makes objects.
        # Currently only used for making NTOs.
        self.chk_files = {}
        
        # Actual fchk file maker objects.
        # These cannot be set here as they depend on our output dir.
        self.fchk_files = {}
        
    def setup(self, output_dir, output_name):
        """
        Perform setup.
        
        Calling this method will set cube objects in the parent report.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        self.setup_fchk(output_dir, output_name)
        self.setup_cubes(output_dir, output_name)
        
    def setup_fchk(self, output_dir, output_name):
        """
        Setup the fchk files which will be used to create cube files.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        # First, get our fchk files (from which cubes are made in Gaussian.
        self.fchk_files = {
            "structure": Chk_to_fchk(
                Path(output_dir, output_name + ".fchk"),
                chk_file = self.chk_file_paths['structure'],
                fchk_file = self.fchk_file_paths['structure'],
                memory = self.memory,
                formchk_executable = self.options['external']['formchk']
            )
        }
        
        # Only get our spin fchk if we have a calc with spin density available.
        if self.chk_file_paths['spin'] != None or self.fchk_file_paths['spin'] != None:
            self.fchk_files['spin'] = Chk_to_fchk(
                Path(output_dir, output_name + ".spin.fchk"),
                chk_file = self.chk_file_paths['spin'],
                fchk_file = self.fchk_file_paths['spin'],
                memory = self.memory
            )
            
        ########
        # NTOs #
        ########
        # NTOs are slightly more complicated because each requires a separate chk file, so there's an additional layer of file conversion.
        if 'NTO' in self.chk_file_paths:
            # CAREFUL: This is very fragile.
            # If this report contains multiple excited state calculations merged together, then the excited states
            # of each individual calc (which here is self.results.excited_states) will have been modified in place.
            # This means that the 'total_level' attribute of each is now incorrect with respect to the order in the
            # chk file. However, the excited_states list itself is not modified, so if we use the index of that list
            # the states should still match. This works but relies on some unusual implementation details so could
            # break in the future.
            for index, excited_state in enumerate(self.result.excited_states):
                if self.calculation is not None:
                    self.chk_files[excited_state.state_symbol + "_NTO"] = Chk_to_NTO_chk.from_calculation(
                        Path(output_dir, "NTOs", "{}.{}.chk".format(output_name, excited_state.state_symbol)),
                        calculation = self.calculation,
                        chk_file = self.chk_file_paths['NTO'],
                        transition = index+1,
                        options = self.options
                    )
                
                else:
                    self.chk_files[excited_state.state_symbol + "_NTO"] = Chk_to_NTO_chk.from_options(
                        Path(output_dir, "NTOs", "{}.{}.chk".format(output_name, excited_state.state_symbol)),
                        chk_file = self.chk_file_paths['NTO'],
                        transition = index+1,
                        options = self.options
                    )
                
                # Also create the relevant fchk file.
                self.fchk_files[excited_state.state_symbol + "_NTO"] = Chk_to_fchk(
                    Path(output_dir, "NTOs", "{}.{}.fchk".format(output_name, excited_state.state_symbol)),
                    chk_file = self.chk_files[excited_state.state_symbol + "_NTO"],
                    memory = self.memory,
                    formchk_executable = self.options['external']['formchk']
                )
            
        
        
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """        
        ################
        # Spin density #
        ################
        if "spin" in self.fchk_files:
            self.report.cubes['spin_density'] = Fchk_to_spin_cube.from_options(
                Path(output_dir, "Spin Density", output_name + ".spin.cube"),
                fchk_file = self.fchk_files['spin'],
                spin_density = "SCF",
                options = self.options,
                memory = self.memory
            )
        
        #################
        # Total density #
        #################
        if "structure" in self.fchk_files:
            self.report.cubes['SCF'] = Fchk_to_density_cube.from_options(
                Path(output_dir, "Density", output_name + ".SCF.cube"),
                fchk_file = self.fchk_files['structure'],
                density_type = "SCF",
                options = self.options,
                memory = self.memory
            )
        
        
        ############
        # Orbitals #
        ############
        if "structure" in self.fchk_files:
            # We need to set images for both alpha and beta orbitals (if we have them).
            for orbital_list in (self.report.result.orbitals, self.report.result.beta_orbitals):
                for orbital in orbital_list:
                    # First, decide what type of orbital we need.
                    if orbital.spin_type == "alpha":
                        cubegen_type = "AMO"
                    elif orbital.spin_type == "beta":
                        cubegen_type = "BMO"
                    else:
                        cubegen_type = "MO"
                    
                    # Save cube.
                    self.report.cubes[orbital.label] = Fchk_to_cube.from_options(
                        Path(output_dir, orbital.label, output_name + ".{}.cube".format(orbital.label)),
                        fchk_file = self.fchk_files['structure'],
                        cubegen_type = cubegen_type,
                        orbital = orbital.level,
                        options = self.options,
                        memory = self.memory
                    )
        
        
        #############
        # Structure #
        #############
        # If we have an orbital cube, we can just reuse this for our structure.
        if "HOMO" in self.report.cubes:
            self.report.cubes['structure'] = self.report.cubes['HOMO']
        elif "HOMO (alpha)" in self.report.cubes:
            self.report.cubes['structure'] = self.report.cubes['HOMO (alpha)']
        elif "structure" in self.fchk_files:
            # No MO cubes available, create one for structure.
            # We'll just use the HOMO to get our cube, as it almost certainly should exist.
            self.report.cubes['structure'] = Fchk_to_cube.from_options(
                Path(output_dir, "Structure", output_name + ".structure.cube"),
                fchk_file = self.fchk_files['structure'],
                cubegen_type = "MO",
                orbital = "HOMO",
                options = self.options,
                memory = self.memory
            )
            
        ########
        # NTOs #
        ########
        for excited_state in self.result.excited_states:
            if excited_state.state_symbol + "_NTO" in self.fchk_files:
                # There is an fchk_file we can use.
                # For NTOs there are two cubes we need, one for the occupied orbital and one for unoccupied.
                self.report.cubes[excited_state.state_symbol + "_NTO_occ"] = Fchk_to_cube.from_options(
                    Path(output_dir, excited_state.state_symbol, output_name + ".NTO.occupied.cube"),
                    fchk_file = self.fchk_files[excited_state.state_symbol + "_NTO"],
                    cubegen_type = "MO",
                    orbital = "HOMO",
                    options = self.options,
                    memory = self.memory
                )
                self.report.cubes[excited_state.state_symbol + "_NTO_unocc"] = Fchk_to_cube.from_options(
                    Path(output_dir, excited_state.state_symbol, output_name + ".NTO.unoccupied.cube"),
                    fchk_file = self.fchk_files[excited_state.state_symbol + "_NTO"],
                    cubegen_type = "MO",
                    orbital = "LUMO",
                    options = self.options,
                    memory = self.memory
                )
        
    
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        # Delete our chk files.
        for chk_file in self.chk_files:
            self.chk_files[chk_file].delete(lazy = True)
        
        # Delete our fchk files.
        for fchk_file in self.fchk_files:
            self.fchk_files[fchk_file].delete(lazy = True)
        
        # Continue.
        super().cleanup()