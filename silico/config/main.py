# General imports.
from copy import deepcopy
import yaml

# Silico imports.
from silico.config.base import Auto_type
from silico.config.configurable.loader import Configurable_list
from silico.config.configurable.base import Configurable
from silico.config.configurable.option import Option, Method_target_option
from silico.config.configurable.options import Options
from silico.config.file.locations import user_config_location
from silico.misc.io import atomic_write
from silico.submit.calculation.turbomole import Turbomole_memory


class Silico_options(Configurable):
    """
    Class for holding main Silico options from various sources.
    """
    
    alignment = Option(help = "The default alignment method to use, MIN: minimal, FAP: furthest atom pair, AA: average angle, AAA: advanced average angle.", choices = ["MIN", "FAP", "AA", "AAA"], default = "MIN")
    angle_units = Option(help = "The default angle units to use, deg: degrees, rad: radians.", choices = ["deg", "rad"], default = "deg")
    
    external = Options(
        help = "Options specifying paths to various external programs that silico may use. If no path is given, then these programs will simply be executed by name (so relying on OS path resolution to find the necessary executables, which is normally fine.)",
        formchk = Option(help = "Gaussian's formchk utility https://gaussian.com/formchk/", default = "formchk"),
        cubegen = Option(help = "Gaussian's cubegen utility https://gaussian.com/cubegen/", default = "cubegen"),
        vmd = Option(help = "VMD (Visual Molecular Dynamics), used for rendering molecules https://www.ks.uiuc.edu/Research/vmd/", default = "vmd"),
        tachyon = Option(help = "The tachyon ray-tracing library, performs the actual rendering. Tachyon is typically packaged with VMD, but often isn't added to the path automatically", default = "tachyon")
    )
    
    skeletal_image = Options(
        help = "Options for controlling the rendering of 2D skeletal images.",
        render_backend = Option(help = "The library to use for rendering.", choices = ["rdkit", "obabel"], default = "rdkit"),
        resolution = Option("The resolution for rendered images.", type = int, default = 1024)
    )
    
    logging = Options(
        help = "Options relating to output of error messages. Note that the final logging level is determined by combining both 'log_level' and 'verbose', so a 'log_level' of 'OFF' and 'verbose' of '2' is equal to 'ERROR'.",
        log_level = Option(help = "The level of messages to output, one of OFF (no logging at all), CRITICAL (fewest messages), ERROR, WARNING, INFO or DEBUG (most messages)", choices = ["OFF", "CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], default = "INFO"),
        verbose = Option(help = "Increase the verbosity of the program by this amount. Each integer increase of verbosity will increase 'log_level' by 1 degree.", type = int, default = 0),
        vmd_logging = Option(help = "Whether to print output from vmd, the rendering program for density plots.", type = bool, default = False)
    )
    
    rendered_image = Options(
        help = "Options for controlling the appearance of 3D molecule images.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        use_existing = Option(help =\
"""If True, previously created files will be reused. If False, new images will rendered, replacing the old.
This is on by default for 3D rendered images because they are expensive (time-consuming) to render.""", type = bool, default = True
        ),
        rendering_style = Option(help =\
"""The render/display mode, changes the appearance of rendered molecules/orbitals.
Possible options are:
 pastel: The default style, uses light, pastel-ish colours for orbitals with low transparency. Normal atom colours.
 light-pastel: Similar to pastel, but with lighter orbital colours.
 dark-pastel: Similar to pastel, but with darker orbital colours, reminiscant of the 'sharp' style.
 sharp: Orbitals are darkly coloured around their edge but have high transparency in their center. Normal atom colours.
 gaussian: A style attempting to mimic that of GaussView.
 vesta: A style attempting to mimic that of VESTA.""", choices = ["pastel", "light-pastel", "dark-pastel", "sharp", "gaussian", "vesta"], default = "pastel"
        ),
        auto_crop = Option(help = "Whether to enable  automatic cropping of excess whitespace around the border of generated images.  If False, overall image rendering is likely to take less time, but molecules may only occupy a small portion of the true image.", type = bool, default = True),
        resolution = Option(help = "The target resolution for rendered images. Higher values will increase image quality, at the cost of increased render time and file size.", type = int, default = 1024),
        
        orbital = Options(help = "Specific options for orbital density plots.",
            cube_grid_size = Option(help =\
"""The size of grid used to generate Gaussian cube files.
This option is passed to cubegen as the 'npts' argument (https://gaussian.com/cubegen/) and controls how detailed MO surfaces are.
Note that generating cube files can take up the majority of total rendering time.
There are a number of valid options here (please see the cubegen manual), the most common options are listed here:
   - 0: The default grid size
   - -2: The 'Coarse' grid, this is the lowest quality grid (and is in fact similar to the default option). It is probably adequate for most uses.
   - -3: The 'Medium' grid, twice the size of the 'Coarse' grid and taking roughly 10x as long. This option is probably adequate for most publication purposes.
   - -4: The 'Fine' grid, twice the size again of 'Medium' and taking a further, rough 10x as long. This option can also consume large amounts of file space, and is probably unnecessary for most purposes.""", type = int, default = 0
            ),
            isovalue = Option(help = "The isovalue to use for rendering orbital density.", type = float, default = 0.02)
        ),
        spin = Options(help = "Specific options for spin density plots.",
            cube_grid_size = Option(help = "The size of the grid use to plot cube data. As cubes of this type are rendered with a smaller isovalue, it is often necessary to use a larger grid size than normal to maintain quality.", type = int, default = -3),
            isovalue = Option(help = "The isovalue to use for plotting spin density.", type = float, default = 0.0004)
        ),
        density = Options(help = "Specific options for total density plots.",
            cube_grid_size = Option(help = "The size of the grid use to plot cube data.", type = int, default = -3),
            isovalue = Option(help = "The isovalue to use for plotting total density.", type = float, default = 0.0004)
        ),
        difference_density = Options(help = "Specific options for excited states difference density plots.",
            isovalue = Option(help = "The isovalue to use for difference density plots.", type = float, default = 0.001)
        ),
        dipole_moment = Options(help = "Specific options for permanent dipole moment plots.",
            scaling = Option(help = "The value (x) to scale the TDM by, where 1 D = x Å.", type = float, default = 1.0)
            
        ),
        transition_dipole_moment = Options(help = "Specific options for transition dipole moment plots.",
            electric_scaling = Option(help = "The value (x) to scale the TEDM by, where 1 D = x Å.", type = float, default = 5.0),
            magnetic_scaling = Option(help = "The value (x) to scale the TMDM by, where 1 au = x Å.", type = float, default = 10.0),
        ),
            
    )
    
    orbital_diagram = Options(help = "Options that control the appearance of orbital energy diagrams.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'all': Limits are set so all orbitals are visible.
  - 'auto': Limits are set automatically so the 0 point (ie, 'X' axis), HOMO and LUMO are clearly visible.
  - 'center': Like 'auto' except the lower limit is set more negative so the HOMO/LUMO are closer to the middle of the diagram.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in eV).""", type = Auto_type, default = "auto" 
        ),
        full_axis_lines = Option(help = "If True, black lines are drawn around the boarder of the diagram. If False, a line is drawn for the Y axis but not for the other 3 sides of the diagram.", type = bool, default = False)
    )
    
    excited_states_diagram = Options(help = "Options that control the appearance of excited states energy diagrams.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'all': Limits are set so all excited states are visible.
  - 'auto': Limits are set automatically so the lowest excited state of each multiplicity is clearly visible (S1, D1, T1, Q1... N1 etc).
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in eV).""", type = Auto_type, default = "all"
        ),
        show_dest = Option(help = "Whether or not to show the dE(ST) label.", type = bool, default = True)
    )
    
    absorption_spectrum = Options(help = "Options for controlling the appearance of simulated UV-Vis like absorption spectra.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in oscillator strength).""", type = Auto_type, default = "auto"
        ),
        x_limits = Option(help =\
"""The method used to set the min/max limits of the X axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [x_min, x_max], where x_min is the most negative value on the X axis, and x_max is the most positive value on the X axis (both in nm).""", type = Auto_type, default = "auto"
        ),
        max_width = Option(help =\
"""The maximum image width in pixels.
Absorption graphs will grow/shrink their width to fit available data, keeping a constant scale (constant pixels to nm ratio) but only up to this maximum.
To disable the maximum width, set to null.""", type = int, default = 1200
        ),
        peak_cutoff = Option(help =\
"""The minimum oscillator strength that a peak must have to be shown in the graph, as a fraction ofthe highest peak.
Set to 0 for no cutoff (all peaks shown), which may results in the graph being extended well beyond the drawn peaks (because many peaks are too small to see).
This option has no effect when using manual x limits.""", type = float, default = 0.01
        ),
        x_padding = Option(help = "The amount (in nm) to extend the x axis past the highest/lowest energy peak.", type = int, default = 40),
        fwhm = Option(help = "The full-width at half-maximum; changes how wide the drawn peaks are. Note that the choice of peak width is essentially arbitrary; only the peak height is given by calculation. Units are eV.", type = float, default = 0.4),
        gaussian_cutoff = Option(help = "The minimum y value to plot using the Gaussian function (controls how close to the x axis we draw the gaussian) as a fraction of the max peak height.", type = float, default = 0.001),
        gaussian_resolution = Option(help = "The spacing between x values to plot using the Gaussian function, in eV. Values that are too large will result in 'curves' made up of a series of straight edges.", type = float, default = 0.01),
        use_jacobian = Option(help = "Whether or not to use the jacobian transformation to correctly scale the y-axis (see J. Phys. Chem. Lett. 2014, 5, 20, 3497)", type = bool, default = True),
        plot_bars = Option(help = "Whether to plot vertical bars for each excited state.", type = bool, default = True),
        plot_peaks = Option(help = "Whether to plot individual Gaussian functions for each excited state.", type = bool, default = False),
        plot_cumulative_peak = Option(help = "Whether to plot the sum of all Gaussian functions (most closely simulates a real spectrum).", type = bool, default = True)
    )
    
    emission_spectrum = Options(help = "Options for controlling the appearance of simulated emission spectra. 'emission_spectrum' and 'absorption_spectrum 'take the same options.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in oscillator strength).""", type = Auto_type, default = "auto"
        ),
        x_limits = Option(help =\
"""The method used to set the min/max limits of the X axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [x_min, x_max], where x_min is the most negative value on the X axis, and x_max is the most positive value on the X axis (both in nm).""", type = Auto_type, default = "auto"
        ),
        max_width = Option(help =\
"""The maximum image width in pixels.
Absorption graphs will grow/shrink their width to fit available data, keeping a constant scale (constant pixels to nm ratio) but only up to this maximum.
To disable the maximum width, set to null.""", type = int, default = 1200
        ),
        peak_cutoff = Option(help =\
"""The minimum oscillator strength that a peak must have to be shown in the graph, as a fraction ofthe highest peak.
Set to 0 for no cutoff (all peaks shown), which may results in the graph being extended well beyond the drawn peaks (because many peaks are too small to see).
This option has no effect when using manual x limits.""", type = float, default = 0.01
        ),
        x_padding = Option(help = "The amount (in nm) to extend the x axis past the highest/lowest energy peak.", type = int, default = 40),
        fwhm = Option(help = "The full-width at half-maximum; changes how wide the drawn peaks are. Note that the choice of peak width is essentially arbitrary; only the peak height is given by calculation. Units are eV.", type = float, default = 0.4),
        gaussian_cutoff = Option(help = "The minimum y value to plot using the Gaussian function (controls how close to the x axis we draw the gaussian) as a fraction of the max peak height.", type = float, default = 0.001),
        gaussian_resolution = Option(help = "The spacing between x values to plot using the Gaussian function, in eV. Values that are too large will result in 'curves' made up of a series of straight edges.", type = float, default = 0.01),
        use_jacobian = Option(help = "Whether or not to use the jacobian transformation to correctly scale the y-axis (see J. Phys. Chem. Lett. 2014, 5, 20, 3497)", type = bool, default = True),
        plot_bars = Option(help = "Whether to plot vertical bars for each excited state.", type = bool, default = True),
        plot_peaks = Option(help = "Whether to plot individual Gaussian functions for each excited state.", type = bool, default = False),
        plot_cumulative_peak = Option(help = "Whether to plot the sum of all Gaussian functions (most closely simulates a real spectrum).", type = bool, default = True)
    )
    
    IR_spectrum = Options(help = "Options for controlling the appearance of simulated IR like vibrational frequency spectra.",
        enable_rendering = Option(help = "Set to False to disable image rendering.", type = bool, default = True),
        y_limits = Option(help =\
"""The method used to set the min/max limits of the Y axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [y_min, y_max], where y_min is the most negative value on the Y axis, and y_max is the most positive value on the Y axis (both in km/mol).""", type = Auto_type, default = "auto"
        ),
        x_limits = Option(help =\
"""The method used to set the min/max limits of the X axis.
Possible options are:
  - 'auto': Limits are set automatically so all peaks are clearly visible.
  - a list of [x_min, x_max], where x_min is the most negative value on the X axis, and x_max is the most positive value on the X axis (both in cm-1).""", type = Auto_type, default = "auto"
        ),
        fwhm = Option(help = "The full-width at half-maximum; changes how wide the drawn peaks are. Note that the choice of peak width is essentially arbitrary; only the peak height is given by calculation. Units are cm-1.", type = int, default = 80),
        max_width = Option(help =\
"""The maximum image width in pixels.
IR spectra will grow/shrink their width to fit available data, keeping a constant scale (constant pixels to nm ratio) but only up to this maximum.
To disable the maximum width, set to null.""", type = int, default = 1500),
        gaussian_cutoff = Option(help = "The minimum y value to plot using the Gaussian function (controls how close to the x axis we draw the gaussian) as a fraction of the max peak height.", type = float, default = 0.001),
        gaussian_resolution = Option(help = "The spacing between x values to plot using the Gaussian function, in eV. Values that are too large will result in 'curves' made up of a series of straight edges.", type = float, default = 1.0)
    )
    
    report = Options(help = "Options for controlling the generation of calculation reports.",
        front_page_image = Option(help = "The image to use for the front page of the report.", choices = ["skeletal", "rendered"], default = "rendered"),
        turbomole = Options(help = "Options that control the running of Turbomole calculations to generate cube files from completed calculations. Note that when reports are created automatically following calculation completion these options will be overridden with the specifics of that calculation.",
            num_CPUs = Option(help = "The number of CPUs with which to run.", type = int, default = 1),
            memory = Option(help = "The amount of memory with which to run.", type = Turbomole_memory, default = Turbomole_memory("1GB")),
            program = Method_target_option(lambda option, config: config.programs, help = "A program definition from the internal library to run.", default = None)
        ),
        gaussian = Options(help = "Options that control the running of Gaussian calculations to generate NTO cube files from completed calculations. Note that when reports are created automatically following calculation completion these options will be overridden with the specifics of that calculation.",
            num_CPUs = Option(help = "The number of CPUs with which to run.", type = int, default = 1),
            memory = Option(help = "The amount of memory with which to run.", type = Turbomole_memory, default = Turbomole_memory("1GB")),
            program = Method_target_option(lambda option, config: config.programs, help = "A program definition from the internal library to run.", default = None),
            # TODO: This needs expanding.
            scratch_path = Option(help = "Path to the top of the scratch directory.", default = "/scratch")
        ),
        cleanup = Option(help =\
"""Whether to delete intermediate files that are written during the report generation process.
Intermediate files include:
   - .cube files.
   - .fchk files.
   - .html files.
   - .css files.
Note that this option will only delete new files written by the program; existing files given by the user are never deleted.""", type = bool, default = True
        ),
        orbital_table = Options(help =\
"""Options which control how many orbitals to print in the MO table.
These numbers are relative to the HOMO for both the min and max.
null can be specified for either/both the min/max to print all available orbitals in that direction.
If both alpha and beta orbitals are available (for unrestricted calculations, for example), then additional orbitals may be printed outside of the given min/max to ensure the given value is met for both sets of orbitals. This is common for triplet calculations, where the alpha and beta frontier MOs are at different levels.
Examples:
  min: -10, max: 16: From HOMO-10 to LUMO+15 inclusive (total of 27 orbitals).
  min: null, max: 11: All orbitals with energy less than or equal to LUMO+10.
  min: 0, max: 1: HOMO and LUMO only.""",
            max = Option(help = "The highest orbital to show in the molecular orbital table.", type = int, default = 16),
            min = Option(help = "The lowest orbital to show in the molecular orbital table.", type = int, default = -15)
        ),
        orbital_image = Options(help = \
"""Options which specify which orbitals to render 3D images of. The default is the HOMO and LUMO.
Orbitals can be specified by level (an index starting at 1 for the most negative orbital, useful if you want images of a particular orbital) and/or by distance from the HOMO (useful for more day-to-day operation.
In addition, excited_state_transition_threshold: can be used to add orbitals that are involved in an excited state transition with a probability above a certain threshold.
Duplicates will be automatically ignored, as will orbitals specified here that are not actually available in the calculation result file.
Alpha and beta are set separately. Beta will be ignored if there are no beta orbitals available. 
Example:
   orbital_levels:
     - 5
   orbital_distance:
     - -1
     - 0
     - 1
   excited_state_transition_threshold: 0.1
   Would render orbital 5, HOMO-1, HOMO, LUMO and any orbitals involved in an excited state transition with a probability 0.1 or greater.""",
            et_transition_threshold = Option(help = "Include orbitals involved in an excited state transition with a probability of this value or greater.", type = float, default = 0.1),
            orbital_levels = Option(help = "Include normal/alpha orbitals by level.", list_type = list, type = int, default = None),
            beta_levels = Option(help = "Include beta orbitals by level.", list_type = list, type = int, default = None),
            orbital_distances = Option(help = "Include normal/alpha orbitals by distance from the HOMO.", list_type = list, type = int, default = [0, 1]),
            beta_distances = Option(help = "Include beta orbitals by distance from the HOMO.", list_type = list, type = int, default = [0, 1])
        ),
        frequency_table = Options(help = "Options that control how many vibrational frequencies to list in the frequencies table. ",
            min_frequency = Option(help = "The most negative frequency to show in the table (remember that frequencies can go negative). 'null' is for no limit. Units are cm-1.", type = float, default = None),
            max_frequency = Option(help = "The most positive frequency to show in the table. 'null' is for no limit. Units are cm-1.", type = float, default = None),
            max_num = Option(help = "The maximum number of frequencies to show in the table.", type = int, default = None)
        )
    )
    
    def __init__(self, validate_now = True, palette = None, destinations = None, programs = None, calculations = None, **kwargs):
        """
        Constructor for Silico_options objects.
        """
        # Set Configurable lists.
        self.destinations = destinations
        self.programs = programs
        self.calculations = calculations
        
        # The palette to use for urwid.
        # A palette is a list of tuples, where the first item of each tuple identifies the name of an attribute, and the remaining specify how that attribute should appear.
        # However, we instead store the palette as a dictionary, as this would appear to fit the format better.
        # The palette can be accessed in the urwid format at 'urwid_palette'.
        self.palette = palette if palette is not None else []
        
        super().__init__(validate_now=validate_now, **kwargs)
        
    # These methods allow the main silico options object to be accessed as a dict.
    # This is largely to provide compatibility with legacy code, where the silico options object was, in fact, a dict.
    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
        
    def __delitem__(self, key):
        # Deleting a configurable option doesn't actually delete it (because that wouldn't make much sense),
        # Instead it reverts that option to its default.
        delattr(self, key)
    
    @property
    def methods(self):
        """
        A list of known methods which have been loaded from config files.
        
        This is a configurable list of destinations, which can be traversed to also find programs and then calculations.
        """
        return self.destinations
    
    @property
    def urwid_palette(self):
        """
        A shortcut for accessing the urwid palette (in a format that can be passed directly to urwid).
        """
        return self.yaml_to_palette(self.palette)
    
    def save(self):
        """
        Save the current value of these options to file, so that they will be reloaded on next program start.
        
        Changed settings are always saved to the user's config file.
        """
        data = yaml.dump(self.dump())
        
        try:
            user_config_location.parent.mkdir(exist_ok = True)
            atomic_write(user_config_location, data)
        
        except FileNotFoundError as e:
            # We lost the race, give up.
            raise Exception("Failed to write settings to file '{}'; one of the parent directories does not exist".format(user_config_location)) from e
            
    @classmethod
    def yaml_to_palette(self, yaml_palette):
        """
        Convert a palette loaded from a YAML config file to a format understood by urwid.
        """
        # The palette (a list of tuples)
        palette = []
        
        # Go through each item in the yaml_palette (a dict).
        for name, values in yaml_palette.items():
            # Add the name to the start of values.
            values.insert(0, name)
            
            # Add as a tuple.
            palette.append(tuple(values))
        
        # Done.
        return palette

    def __deepcopy__(self, memo):
        """
        Overriding deep copy.
        
        Because Silico_options can contain large lists of Configurables (complex objects), real deepcopy can be very slow.
        We overcome this by excluding these lists.
        """
        # TODO: Need to review whether this is necessary and how it could be improved.
        destinations = self.destinations
        self.destinations = []
        programs = self.programs
        self.programs = []
        calculations = self.calculations
        self.calculations = []
        # TMP remove __deepcopy__ so we don't recurse.
        copyfunc = self.__deepcopy__
        self.__deepcopy__ = None
        
        new = deepcopy(self, memo)

        
        # Restore.
        self.destinations = destinations
        new.destinations = destinations
        self.programs = programs
        new.programs = programs
        self.calculations = calculations
        new.calculations = calculations
        self.__deepcopy__ = copyfunc
        
        # And return.
        return new
            