# General imports
from itertools import chain
from pathlib import Path

# Silico imports.
from silico.image.vmd import Spin_density_image_maker, Orbital_image_maker,\
    Combined_orbital_image_maker, Structure_image_maker, Dipole_image_maker,\
    Density_image_maker, Differential_density_image_maker
from silico.image.orbitals import Orbital_diagram_maker
from silico.exception.base import Result_unavailable_error
from silico.image.excited_states import Excited_states_diagram_maker
from silico.image.spectroscopy import Absorption_graph_maker,\
    Emission_graph_maker, Frequency_graph_maker
from silico.result.excited_states import Excited_state_list
from silico.image.graph import Convergence_graph_maker
from silico.result.molecular_orbitals import Molecular_orbital_list
import silico.report.image
from silico.image.structure import Skeletal_image_maker
from silico.file.input import Silico_coords
from silico.report.image.turbomole.base import Turbomole_anadens_setup


class Report():
    """
    An enhanced report set object that contains graphics and other objects for building graphical reports.
    """
    
    def __init__(self, result, *, options, calculation = None):
        """
        Constructor for Report objects.
        
        This Report class is calculation program agnostic, which is important because the results from multiple different calculation programs can be merged together and supplied as a merged result.
        However, image generation is deeply integrated with the calculation program that produced the individual results, so that is handled separately by Cube_setup objects which do depend on each calculation program type.
        
        :param result: A parsed result set object to render a report of.
        :param options: A silico Config dictionary which contains various options that control the appearance of this report.
        """    
        # Save our result set object.
        self.result = result
        
        # Save our image maker options.
        self.options = options
        
        # Two dictionaries of VMD image and cube maker objects which we use to render 3D images.
        self.images = {}
        self.cubes = {}
        
        # Path to the directory in which the current report is being written.
        self.report_directory = None
        
        # Decide which extra orbitals we want.
        self._init_get_orbitals_to_render(
            orbital_distances = options['report']['orbital_image']['orbital_distances'],
            beta_distances = options['report']['orbital_image']['beta_distances'],
            orbital_levels = options['report']['orbital_image']['orbital_levels'],
            beta_levels = options['report']['orbital_image']['beta_levels'],
            excited_state_transition_threshold = options['report']['orbital_image']['et_transition_threshold']
        )            
        
        # Get cube makers for each metadata.
        # We need to be careful with what we request from each image maker, especially for Turbomole
        # eg, We can't just ask for orbitals from all Turbomole makers because this is wasteful if 
        # we only need spin images from one of the makers; Turbomole will still render the orbital images
        # which we won't use. Worse, if this image maker is unrestricted but the main result is restricted,
        # the image maker will fail because it will expect restricted file names for orbital cubes from turbomole
        # (we can't control what turbomole names cube files), but these won't be there.
        #
        # We also move backwards through our results, this ensures that earlier calculations have precedence.
        self.image_setters = []
        for index, result in enumerate(reversed(self.result.results)):
            # Only request orbitals from the main result (first, but remember we are travelling backwards thro metadatas).
            do_orbitals = index == len(self.result.metadatas) -1
            # Only request spin images if available.
            do_spin = result.metadata.multiplicity != 1 and result.metadata.orbital_spin_type == "unrestricted"
            
            # Setup.
            self.image_setters.append(
                silico.report.image.from_result(self, result = result, do_orbitals = do_orbitals, do_spin = do_spin, options = options, calculation = calculation)
            )
            
        # Also create a cube setter for anadens cubes.
        # Anadens cubes are specific to calculations performed with Turbomole, but they are handled differently because they require multiple calculations (a ground state and excited state) to work.
        # Even if no Turbomole calculations are available here, this class will function fine.
        self.image_setters.append(Turbomole_anadens_setup(self, options = options, calculation = calculation))
            
        
    @property
    def rotations(self):
        """
        """
        return self.result.alignment.rotations
        
    def _init_get_orbitals_to_render(self,
            orbital_distances = None,
            beta_distances = None,
            orbital_levels = None,
            beta_levels = None,
            excited_state_transition_threshold = None
        ):
        """
        Init helper function which decides which orbitals to render later.
        """
        orbital_distances = orbital_distances if orbital_distances is not None else []
        beta_distances = beta_distances if beta_distances is not None else []
        orbital_levels = orbital_levels if orbital_levels is not None else []
        beta_levels = beta_levels if beta_levels is not None else []
        excited_state_transition_threshold = excited_state_transition_threshold if excited_state_transition_threshold is not None else 2
        
        try:
            # Init our list.
            self.orbitals_to_render = Molecular_orbital_list()
            
            # Now add those orbitals we've been asked to find.
            self.orbitals_to_render.extend(
                chain(
                    # First alpha.
                    chain.from_iterable(self.result.molecular_orbitals.search(HOMO_difference = HOMO_difference) for HOMO_difference in orbital_distances),
                    chain.from_iterable(self.result.molecular_orbitals.search(level = level) for level in orbital_levels),
                    
                    # Now beta.
                    chain.from_iterable(self.result.beta_orbitals.search(HOMO_difference = HOMO_difference) for HOMO_difference in beta_distances),
                    chain.from_iterable(self.result.beta_orbitals.search(level = level) for level in beta_levels),
                    
                    # Also add orbitals that are involved in excited state transitions.
                    chain.from_iterable(
                        (transition.starting_mo, transition.ending_mo) for excited_state in self.result.excited_states for transition in excited_state.transitions if transition.probability >= excited_state_transition_threshold
                    )
                )
            )
            
            # Now remove duplicates and reorder.
            self.orbitals_to_render = self.orbitals_to_render.ordered()
            
        except Exception:
            self.orbitals_to_render = []
            raise
          
    def setup_cubes(self, output_dir, output_name):
        """
        Setup the cube files which will be used to render images.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        for image_setup in self.image_setters:
            image_setup.setup(output_dir, output_name)
    
    def relative_image(self, image, sub_image = 'file'):
        """
        Get a path to a particular image relative to the report file.
        
        :param image: Name of the image to get (should correspond to an item in self.images)
        :param sub_image: Optional name of the sub image to get if image represents more than one real image.
        """
        return self.images[image].relative_path(sub_image, output_base = self.report_directory)
    
    def setup_images(self, output_dir, output_name):
        """
        Set the options that will be used to create images from this object.
        
        :param output_dir: A pathlib Path object to the directory within which our files should be created.
        :param output_name: A string that will be used as the start of the file name of the files we create.
        """
        ################
        # Spin density #
        ################
        # Only set spin image maker objects if we have a cube we can use.
        if "spin_density" in self.cubes:
            self.images['positive_spin_density'] = Spin_density_image_maker.from_options(
                Path(output_dir, "Spin Density", output_name + ".spin_pos.jpg"),
                cube_file = self.cubes['spin_density'],
                rotations = self.rotations,
                spin = "positive",
                options = self.options)
            
            self.images['negative_spin_density'] = Spin_density_image_maker.from_options(
                Path(output_dir, "Spin Density", output_name + ".spin_neg.jpg"),
                cube_file = self.cubes['spin_density'],
                rotations = self.rotations,
                spin = "negative",
                options = self.options)
            
            self.images['spin_density'] = Spin_density_image_maker.from_options(
                Path(output_dir, "Spin Density", output_name + ".spin_both.jpg"),
                cube_file = self.cubes['spin_density'],
                rotations = self.rotations,
                spin = "both",
                options = self.options)
        
        #################
        # Total density #
        #################
        self.images['SCF'] = Density_image_maker.from_options(
            Path(output_dir, "Density", output_name + ".SCF.jpg"),
            cube_file = self.cubes['SCF'],
            rotations = self.rotations,
            options = self.options)

        ############
        # Orbitals #
        ############
        # We need to set images for both alpha and beta orbitals (if we have them).
        for orbital_list in (self.result.molecular_orbitals, self.result.beta_orbitals):
            # First set individual orbital images.
            for orbital in orbital_list:
                # Save our orbital image.
                self.images[orbital.label] = Orbital_image_maker.from_options(
                    Path(output_dir, orbital.label, output_name + ".{}.jpg".format(orbital.label)),
                    cube_file = self.cubes[orbital.label],
                    rotations = self.rotations,
                    options = self.options)
                
            # Next set combined images and graphs.
            try:
                # Now get our spin type.
                if orbital_list.spin_type == "alpha":
                    spin_type = "alpha_"
                elif orbital_list.spin_type == "beta":
                    spin_type = "beta_"
                else:
                    spin_type = ""
                    
                # Save our orbital diagram.
                self.images[spin_type + 'orbital_energies'] = Orbital_diagram_maker.from_options(
                    Path(output_dir, "Orbital Diagram", output_name + ".{}orbitals.png".format(spin_type)),
                    molecular_orbitals = orbital_list,
                    options = self.options)
                
                # A version of the diagram with only the HOMO/LUMO
                self.images[spin_type + 'HOMO_LUMO_energies'] = Orbital_diagram_maker.from_options(
                    Path(output_dir, "Orbital Diagram", output_name + ".{}HOMO_LUMO.png".format(spin_type)),
                    molecular_orbitals = type(orbital_list)(orbital for orbital in orbital_list if orbital.HOMO_difference == 0 or orbital.HOMO_difference == 1),
                    options = self.options)
                
                # Also get our HOMO/LUMO combined image.
                # Get our FMOs.
                HOMO = orbital_list.get_orbital(HOMO_difference = 0)
                LUMO = orbital_list.get_orbital(HOMO_difference = 1)
                
                self.images[spin_type + 'HOMO_LUMO'] = Combined_orbital_image_maker.from_options(
                    Path(output_dir, "HOMO LUMO", output_name + ".{}HOMO_LUMO.jpg".format(spin_type)),
                    HOMO_cube_file = self.cubes[HOMO.label],
                    LUMO_cube_file = self.cubes[LUMO.label],
                    rotations = self.rotations,
                    options = self.options)
                
            except Result_unavailable_error:
                # We couldn't find our HOMO/LUMO.
                pass
        
        
        if 'structure' in self.cubes:
            #############
            # Structure #
            #############
            # We'll save both the aligned and unaligned structure.
            self.images['structure'] = Structure_image_maker.from_options(
                Path(output_dir, "Structure", output_name + ".structure.jpg"),
                cube_file = self.cubes['structure'],
                options = self.options)
            
            self.images['aligned_structure'] = Structure_image_maker.from_options(
                Path(output_dir, "Structure", output_name + ".structure.jpg"),
                cube_file = self.cubes['structure'],
                rotations = self.rotations,
                options = self.options)
        
        
            #######################
            # Dipole moment (PDM) #
            #######################
            self.images['dipole_moment'] = Dipole_image_maker.from_options(
                Path(output_dir, "Dipole Moment", output_name + ".dipole.jpg"),
                cube_file = self.cubes['structure'],
                dipole_moment = self.result.dipole_moment,
                rotations = self.rotations,
                options = self.options)
            
        ################
        # 2D Structure #
        ################
        self.images['skeletal'] = Skeletal_image_maker.from_options(
            Path(output_dir, "Structure", output_name + ".skeletal.png"),
            coords = Silico_coords.from_xyz(self.result.atoms.to_xyz()),
            options = self.options
            )
        # NOTE: We ask for the skeletal image here because if report: front_page == "rendered" (the default), then this image is never created.
        # This could be a shame for the user as the image might be useful even if it's not used in the report (and it takes no time to make so).
        self.images['skeletal'].get_image()
        
        ####################################################
        # Excited states & transition dipole moments (TDM) #
        ####################################################
        for excited_state in self.result.excited_states:
            # Set PDM.
            # Work out what we'll name our files.
            file_name = "{}_dipole".format(excited_state.state_symbol)
            sub_dir_name = "{} Transition Dipole Moment".format(excited_state.state_symbol)
                
            # Get our image.
            if "structure" in self.cubes:
                self.images[file_name] = Dipole_image_maker.from_options(
                    Path(output_dir, sub_dir_name, output_name + ".{}.jpg".format(file_name)),
                    cube_file = self.cubes['structure'],
                    dipole_moment = excited_state.transition_dipole_moment,
                    rotations = self.rotations,
                    options = self.options)
                
            # If we have differential density cubes, create images that can use them.
            if excited_state.state_symbol + "_differential_density" in self.cubes:
                self.images[excited_state.state_symbol + "_differential_density"] = Differential_density_image_maker.from_options(
                    Path(output_dir, excited_state.state_symbol, output_name + ".{}_differential_density.jpg".format(excited_state.state_symbol)),
                    cube_file = self.cubes[excited_state.state_symbol + "_differential_density"],
                    rotations = self.rotations,
                    options = self.options
                )
                
            # If we have NTO cubes, create images for those too.
            if excited_state.state_symbol + "_NTO_occ" in self.cubes and excited_state.state_symbol + "_NTO_unocc" in self.cubes:
                self.images[excited_state.state_symbol + "_NTO"] = Combined_orbital_image_maker.from_options(
                    Path(output_dir, excited_state.state_symbol, output_name + ".{}_NTO.jpg".format(excited_state.state_symbol)),
                    # For NTOs, we swap the apparent HOMO and LUMO.
                    # This is because we're really interested in how the excited state compares to the ground (rather than the reverse),
                    # So we'll set the HOMO to the orbital containing the excited electron.
                    HOMO_cube_file = self.cubes[excited_state.state_symbol + "_NTO_unocc"],
                    LUMO_cube_file = self.cubes[excited_state.state_symbol + "_NTO_occ"],
                    rotations = self.rotations,
                    options = self.options
                )
                
            
        # Also set our states diagram.
        self.images['excited_state_energies'] = Excited_states_diagram_maker.from_options(
            Path(output_dir, output_name + ".excited_states.png"),
            excited_states = self.result.excited_states,
            ground_state = self.result.ground_state,
            options = self.options
        )
        
        # Then our simulated absorption graph.
        self.images['simulated_absorption_graph'] = Absorption_graph_maker.from_options(
            Path(output_dir, output_name + ".simulated_absorption_spectrum.png"),
            excited_states = self.result.excited_states,
            options = self.options
        )
        
        
        #################################
        # Vertical & adiabatic emission #
        #################################
        for emission_type in (self.result.vertical_emission, self.result.adiabatic_emission):
            for emission_multiplicity in emission_type:
                emission = emission_type[emission_multiplicity]
                
                # First our states diagram.
                self.images['{}_{}_emission_energies'.format(emission.transition_type, emission.state_symbol)] = Excited_states_diagram_maker.from_options(
                    Path(output_dir, output_name + ".{}_{}_emission_states.png".format(emission.transition_type, emission.state_symbol)),
                    excited_states = Excited_state_list([emission]),
                    ground_state = self.result.ground_state,
                    options = self.options)
                
                # Now emission spectrum.
                self.images['simulated_{}_{}_emission_graph'.format(emission.transition_type, emission.state_symbol)] = Emission_graph_maker.from_options(
                    Path(output_dir, output_name + ".simulated_{}_{}_emission_spectrum.png".format(emission.transition_type, emission.state_symbol)),
                    excited_states = Excited_state_list([emission]),
                    options = self.options)
                    
        
        
        #############################
        # Energy convergence graphs #
        #############################
        for energy in (self.result.SCF_energies, self.result.MP_energies, self.result.CC_energies):
            self.images['{}_convergence_graph'.format(energy.energy_type)] = Convergence_graph_maker.from_options(
                Path(output_dir, output_name + ".{}_graph.png".format(energy.energy_type)),
                energies = energy,
                options = self.options
            )
        
        
        ##############
        # Vibrations #
        ##############
        # First our states diagram.
        self.images['simulated_IR_graph'] = Frequency_graph_maker.from_options(
            Path(output_dir, output_name + ".simulated_frequencies.png"),
            vibrations = self.result.vibrations,
            options = self.options
        )
                
    def cleanup(self):
        """
        Remove any intermediate files that may have been created by this object.
        """
        # Delete all our cubes.
        for name, cube in self.cubes.items():
            if cube is not None:
                cube.delete(lazy = True)
                
        # Cleanup image setup.
        for image_setup in self.image_setters:
            image_setup.cleanup()
    
    def _write(self, output, **kwargs):
        """
        Write the various elements of this report to file.
        
        :param output: Path to a directory in which the report will be written.
        """
        self.report_directory = output
        # Base directory for our images.
        image_dir = Path(self.report_directory, "image")
        # The base name for our images.
        #image_base_name = Path(self.result.metadata.name).with_suffix("").name
        image_base_name = self.result.metadata.name
        
        # Make our output directory.
        try:
            self.report_directory.mkdir(parents = True)
        except FileExistsError:
            # This will happen when the dir already exists, which is fine.
            pass
        
        # And our image dir.
        try:
            image_dir.mkdir(parents = True)
        except FileExistsError:
            # This will happen when the dir already exists, which is fine.
            pass
        
        # Setup cubes and images.
        self.setup_cubes(image_dir, image_base_name)
        self.setup_images(image_dir, image_base_name)
        
        
        
    def write(self, output, **kwargs):
        """
        Write this HTML_report to file.
        """
        self._write(output, **kwargs)
        
        # If we've been asked to remove intermediate files, do so.
        if self.options['report']['cleanup']:
            self.cleanup()
