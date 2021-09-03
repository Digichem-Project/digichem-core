from mako.lookup import TemplateLookup

import silico
from silico.exception import Configurable_exception
from silico.config.configurable.option import Option
from silico.submit.calculation import Concrete_calculation
from silico.config.configurable.options import Options

class Keyword():
    """
    A class that represents a Gaussian keyword.
    
    Keywords appear in the route section of a Gaussian input file and instruct Gaussian which calculation to perform. They have any of the following general syntax:
        keyword
        keyword=option
        keyword=(option=value)
        keyword=(option1, option2=value)
    """
    
    def __init__(self, keyword, *options):
        """
        Constructor for Gaussian Keyword objects.
        
        :param keyword: The Gaussian keyword.
        :param options: An optional number of dicts which contain options for the keyword. 
        """
        self.keyword = keyword
        self.options = options
        
    def __str__(self):
        """
        Stringify this Keyword.
        
        The returned string is suitable for inclusion in a Gaussian input file.
        """
        options_strings = []
        for option in self.options:
            for option_name, option_value in option:
                if option_value == "":
                    # 'Blank' option value, just add the option name.
                    options_strings.append(option_name)
                
                elif option_value != None:
                    # Add the name and the value, unless the value is None (in which case we add neither).
                    options_strings.append("{}={}".format(option_name, option_value))
        
        if len(options_strings) == 0:
            return self.keyword
        
        else:
            return "{}=({})".format(self.keyword, ", ".join(options_strings))
        

class Gaussian(Concrete_calculation):
    """
    DFT (density functional theory) calculations with Gaussian.
    """
    # Identifying handle.
    CLASS_HANDLE = ("Gaussian",)
    
    # A list of strings describing the expected input file types (file extensions) for calculations of this class. The first item of this list will be passed to obabel via the -o flag. 
    INPUT_FILE_TYPES = ["gau", "com", "gjf", "gjc"]
    
    # The format of the output file containing coordinates.
    OUTPUT_COORD_TYPE = "log"

    # Configurable options.
    CPU_list = Option(help = "A list of integers specifying specific CPUs to use for the calculation, starting at 0. CPU_list and num_CPUs are mutually exclusive", exclude = ("num_CPUs",), default = (), type = tuple)
    num_CPUs = Option(help = "An integer specifying the number of CPUs to use for this calculation. CPU_list and num_CPUs are mutually exclusive", exclude = ("CPU_list",), default = 1, type = int)
    functional = Option(help = "The DFT functional to use", type = str)
    basis_set = Option(help = "The basis set to use. 'Gen' or 'GenECP' should be given here if an external basis set is to be used", type = str)
    _external_basis_sets = Option(
        "external_basis_sets",
        help = "A list of external basis sets to use. The order given here is the order the basis sets will be appended to the input file",
        type = tuple,
        default = ()
    )
    _external_ECPs = Option(
        "external_ECPs",
        help = "A list of external ECPs (effective core potentials) to use",
        type = tuple,
        default = ()
    )
    _multiplicity = Option("multiplicity", help = "Forcibly set the molecule multiplicity. Leave blank to use the multiplicity given in the input file", default = None, type = int)
    _charge = Option("charge", help = "Forcibly set the molecule charge. Leave blank to use the charge given in the input file", default = None, type = float)
    solvent = Option(help = "Name of the solvent to use for the calculation (the model used is SCRF-PCM)", default = None, type = str)
    convert_chk = Option(help = "Whether to create an .fchk file at the end of the calculation", default = True, type = bool)
    keep_chk = Option(help = "Whether to keep the .chk file at the end of the calculation. If False, the .chk file will be automatically deleted, but not before it is converted to an .fchk file (if convert_chk is True)", default = False, type = bool)
    keep_rwf = Option(help = "Whether to keep the .rwf file at the end of the calculation. If False, the .rwf file will be automatically deleted", default = False, type = bool)
    
    optimisation = Options(
        help = "Options that control optimisations.",
        on = Option(help = "Whether to perform an optimisation. If excited states are also being calculated, then the excited state given by 'root' will be optimised", default = False, type = bool),
        options = Option(help = "Additional options to specify.", type = dict, default = {})
    )
    frequency = Options(
        help = "Options that control calculation of vibrational frequencies.",
        on = Option(help = "Whether to calculation vibrational frequencies.", default = False, type = bool),
        options = Option(help = "Additional options to specify.", type = dict, default = {})
    )
    DFT_excited_states = Options(
        help = "Options for calculation of excited states with DFT (TDA or TD-DFT)",
        multiplicity = Option(help = "Multiplicity of the excited states to calculate.", default = None, choices = (None, "Singlet", "Triplet", "50-50")),
        nstates = Option(help = "The number of excited states to calculate. If 50-50 is given as the multiplicity, this is the number of each multiplicity to calculate. If 0 is given, no excited states will be calculated", type = int, default = 0),
        root = Option(help = "The 'state of interest', the meaning of which changes depending on what type of calculation is being performed. For example, if an optimisation is being performed, this is the excited state to optimise", type = int, default = None),
        TDA = Option(help = "Whether to use the Tamm–Dancoff approximation.", type = bool, default = False),
        options = Option(help = "Additional options to specify.", type = dict, default = {})
    )
    keywords = Option(help = "Additional Gaussian route keywords. These are written to the input file with only minor modification ('keyword : value' becomes 'keyword=value'), so any option valid to Gaussian can be given here", default = {'Population': 'Regular', 'Density': 'Current'}, type = dict)            
    
    @property
    def charge(self):
        """
        The molecule/system charge that we'll actually be using in the calculation.
        
        Unlike the charge attribute, this property will translate "auto" to the actual charge to be used.
        """
        return int(self._charge if self._charge is not None else self.input_coords.implicit_charge)
    
    @property
    def multiplicity(self):
        """
        The molecule/system multiplicity that we'll actually be using in the calculation.
        
        Unlike the multiplicity attribute, this property will translate "auto" to the actual multiplicity to be used.
        """
        return int(self._multiplicity if self._multiplicity is not None else self.input_coords.implicit_multiplicity)
    
    @property
    def external_ECPs(self):
        """
        The list of basis set Configurable objects we'll be using in the calculation for effective core potentials.
        
        This property will translate the names of the basis sets, under self._extended_ECPs, to the actual objects.
        """
        known_ECPs = self.silico_options.effective_core_potentials
        try:
            return [known_ECPs.get_config(ECP) for ECP in self._external_ECPs]
        except Exception:
            raise Configurable_exception(self, "could not load external ECP")
        
    @property
    def external_basis_sets(self):
        """
        The list of basis set Configurable objects we'll be using in the calculation.
        
        This property will translate the names of the basis sets, under self._extended_basis_sets, to the actual objects.
        """
        try:
            return [self.silico_options.basis_sets.get_config(basis_set) for basis_set in self._external_basis_sets]
        except Exception:
            raise Configurable_exception(self, "could not load external basis set")
        
    @property
    def model_chemistry(self):
        """
        The 'model chemistry' to be used by the calculation, this is a string containing both the functional and basis set (separated by /).
        """    
        model = ""
        # Add the functional.
        if self.functional is not None:
            model += self.functional
        # Add a slash if we have both functional and basis set.
        if self.functional is not None and self.basis_set is not None:
            model += "/"
        # And finally the basis set.
        if self.basis_set is not None:
            model += self.basis_set
            
        return model
    
    @property
    def calculation_keywords(self):
        """
        Get a string containing the Gaussian calculation keywords and associated options for this calculation.
        """
        keyword_sections = []
        
        # Optimisations.
        if self.optimisation['on']:
            keyword_sections.append(str(Keyword("Opt", self.optimisation['options'])))
            
        # Frequencies.
        if self.frequency['on']:
            keyword_sections.append(str(Keyword("Freq", self.frequency['options'])))
            
        # Excited states.
        if self.DFT_excited_states['nstates'] != 0:
            # First, build our options.
            options = {
                'nstates': self.DFT_excited_states['nstates'],
                'multiplicity': self.DFT_excited_states['multiplicity'],
                'root': self.DFT_excited_states['root']
            }
                
            # Add our keyword, which changes base on whether we're using TDA or TD-DFT.
            keyword_sections.append(str(Keyword("TDA" if self.DFT_excited_states['TDA'] else "TD", options, self.DFT_excited_states['options'])))
            
        # Merge and return.
        return " ".join(keyword_sections)

    @property
    def route_section(self):
        """
        Get a Gaussian input file route section from this calculation target.
        """
        # Assemble our route line.
        # Add calc keywords.
        route_parts = list(self.calculation_keywords)
        
        # Model chemistry
        route_parts.append(self.model_chemistry)
        
        # Solvent.
        if self.solvent is not None:
            route_parts.append(self.keyword_to_string("SCRF", {"Solvent": self.solvent}))
        
        # Finally, add any free-form options.
        for keyword in self.keywords:
            if self.options[option] == "":
                # Blank option, just add the keyword.
                route_parts.append(option)
            
            # Skip None options.
            elif self.options[option] is not None: 
                # Option with options, add both.
                route_parts.append("{}={}".format(option, self.options[option]))
                
        # Convert to string and return.
        return " ".join(route_parts)
    
    
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Concrete_calculation._actual):
        """
        Inner class for calculations.
        """
        
        def __init__(self, *args, **kwargs):
            """
            Constructor for Gaussian calcs.
            """
            super().__init__(*args, **kwargs)
            self.chk_file_name = None
            self.com_file_name = None
            self.rwf_file_name = None
            self.com_file_body = None
        
        def prepare(self, output, input_coords):
            """
            Prepare this calculation for submission.
            
            :param output: Path to a directory to perform the calculation in.
            :param input_coords: A Silico_input object containing the coordinates on which the calculation will be performed.
            """            
            # Call parent.
            super().prepare(output, input_coords)
            
            # Decide on our file names.
            self.chk_file_name = self.safe_name(self.molecule_name + ".chk")
            self.rwf_file_name = self.safe_name(self.molecule_name + ".rwf")
            self.com_file_name = self.safe_name(self.molecule_name + ".com")
            
            # Get and load our com file template.
            self.com_file_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/gaussian/input_file.mako").render_unicode(calculation = self)
                
