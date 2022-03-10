from mako.lookup import TemplateLookup
import deepmerge

import silico
from silico.config.configurable.option import Option
from silico.submit.calculation import Concrete_calculation
from silico.config.configurable.options import Options
from silico.submit.basis import BSE_basis_set
from silico.file.input.directory import Calculation_directory_input
from silico.submit.calculation.base import AI_calculation_mixin


class Keyword():
    """
    A class that represents a Gaussian keyword.
    
    Keywords appear in the route section of a Gaussian input file and instruct Gaussian which calculation to perform. They have any of the following general syntax:
        keyword
        keyword=option
        keyword=(option=value)
        keyword=(option1, option2=value)
    """
    
    # TODO: Add better support for gaussian keyword options that don't contain a sub option (eg, TDA=(50-50))
    def __init__(self, keyword, *options):
        """
        Constructor for Gaussian Keyword objects.
        
        :param keyword: The Gaussian keyword.
        :param options: An optional number of options for the keyword. If a dict is given, name: value will be converted to keyword=(name=value) (unless value is None, in which case the option is ignored), otherwise each option will be converted to keyword=(name)
        """
        self.keyword = keyword
        self.options = {}
        
        # Deal with options, which could contain several different types.
        for option in options:
            if isinstance(option, dict):
                # Dictionary, merge.
                deepmerge.always_merger.merge(self.options, option)
                
            elif option is not None:
                # Convert to string, assume this is an option without a value.
                self.options[str(option)] = ""
        
    def __str__(self):
        """
        Stringify this Keyword.
        
        The returned string is suitable for inclusion in a Gaussian input file.
        """
        options_strings = []
        for option_name, option_value in self.options.items():
            if option_value == "":
                # 'Blank' option value, just add the option name.
                options_strings.append(option_name)
            
            # TODO: This check may be unnecessary.
            elif option_value != None:
                # Add the name and the value, unless the value is None (in which case we add neither).
                options_strings.append("{}={}".format(option_name, option_value))
        
        if len(options_strings) == 0:
            return self.keyword
        
        else:
            return "{}=({})".format(self.keyword, ", ".join(options_strings))


class Gaussian(Concrete_calculation, AI_calculation_mixin):
    """
    Calculations with Gaussian.
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
    DFT = Options(help = "Options for DFT.",
        functional = Option(help = "The DFT functional to use. Specifying an option here will enable DFT.", type = str),
        empirical_dispersion = Option(help = "Optional empirical dispersion method to use, note that not all methods are suitable with all functions.", choices = ("PFD", "GD2", "GD3", "GD3BJ", None), default = None)
    )
    post_HF_method = Option(help = "The name of a post-HF, calculation method to perform.", choices = ("MP2", "MP3", "MP4", "MP4(DQ)", "MP4(SDQ)", "MP5", "CCD", "CCSD", None), default = None)
    unrestricted = Option(help = "Whether to perform an unrestricted calculation", type = bool, default = False)
    basis_set = Options(help = "The basis set to use.",
        internal = Option(help = "The name of a basis set built in to Gaussian, see Gaussian manual for allowed values.", type = str, exclude = "exchange"),
        exchange = Option(help = "The definition of a (number of) basis sets to use from the Basis Set Exchange (BSE), in the format 'basis set name': 'applicable elements' (for example: '6-31G(d,p)': '1,3-4,B-F')", type = BSE_basis_set, dump_func = lambda option, configurable, value: value.definition if value is not None else {}, exclude = "internal", edit_vtype = "dict")
    )
#     _multiplicity = Option("multiplicity", help = "Forcibly set the molecule multiplicity. Leave blank to use the multiplicity given in the input file", default = None, type = int)
#     _charge = Option("charge", help = "Forcibly set the molecule charge. Leave blank to use the charge given in the input file", default = None, type = float)
    solvent = Option(help = "Name of the solvent to use for the calculation (the model used is SCRF-PCM)", default = None, type = str)
    convert_chk = Option(help = "Whether to create an .fchk file at the end of the calculation", default = True, type = bool)
    keep_chk = Option(help = "Whether to keep the .chk file at the end of the calculation. If False, the .chk file will be automatically deleted, but not before it is converted to an .fchk file (if convert_chk is True)", default = False, type = bool)
    keep_rwf = Option(help = "Whether to keep the .rwf file at the end of the calculation. If False, the .rwf file will be automatically deleted", default = False, type = bool)
    
    optimisation = Options(
        help = "Options that control optimisations.",
        calculate = Option(help = "Whether to perform an optimisation. If excited states are also being calculated, then the excited state given by 'root' will be optimised", default = False, type = bool),
        options = Option(help = "Additional options to specify.", type = dict, default = {})
    )
    frequency = Options(
        help = "Options that control calculation of vibrational frequencies.",
        calculate = Option(help = "Whether to calculate vibrational frequencies.", default = False, type = bool),
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
    keywords = Option(help = "Additional Gaussian route keywords. These are written to the input file with only minor modification ('keyword: option' becomes 'keyword=(option)' and 'keyword: {option: value}' becomes 'keyword=(option=value)'), so any option valid to Gaussian can be given here", default = {'Population': 'Regular', 'Density': 'Current'}, type = dict)            
    
    @property
    def method(self):
        """
        A string describing the method to use (either a DFT functional or a post-HF method).
        """
        if self.DFT['functional'] is not None:
            return self.DFT['functional']
        
        else:
            return self.post_HF_method
    
    @property
    def basis_set_name(self):
        """
        A descriptive label of the basis set being used for the calculation.
        """
        if self.basis_set['builtin'] is not None:
            return self.basis_set['builtin']
        
        elif self.basis_set['external'] is not None:
            return str(self.basis_set['external'])
        
        else:
            return None
    
    @property
    def calculation_keywords(self):
        """
        Get a string containing the Gaussian calculation keywords and associated options for this calculation.
        """
        keyword_sections = []
        
        # Optimisations.
        if self.optimisation['calculate']:
            keyword_sections.append(str(Keyword("Opt", self.optimisation['options'])))
            
        # Frequencies.
        if self.frequency['calculate']:
            keyword_sections.append(str(Keyword("Freq", self.frequency['options'])))
            
        # Excited states.
        if self.DFT_excited_states['nstates'] != 0:
            # First, build our options.
            options = {
                self.DFT_excited_states['multiplicity']: "",
                'nstates': self.DFT_excited_states['nstates'],
                'root': self.DFT_excited_states['root']
            }
                
            # Add our keyword, which changes base on whether we're using TDA or TD-DFT.
            keyword_sections.append(str(Keyword("TDA" if self.DFT_excited_states['TDA'] else "TD", options, self.DFT_excited_states['options'])))
            
        # Merge and return.
        return " ".join(keyword_sections)
    
    
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
            
        @property
        def model_chemistry(self):
            """
            The 'model chemistry' to be used by the calculation, this is a string containing both the method/functional and basis set (separated by /).
            """
            # First, determine the basis set being used.
            # This label is distinct from both the basis_set property and the basis_set.internal option.
            if self.basis_set['internal'] is not None:
                basis_set = self.basis_set['internal']
                
            elif self.basis_set['exchange'] is not None:
                # Our basis set label will be either gen or genECP (depending on whether we have any ECPs).
                basis_set = "gen" if not self.basis_set['exchange'].has_ECPs(self.input_coords.elements) else "genECP"
                
            else:
                basis_set = None
            
            model = ""
            method = self.method
            # Add the functional.
            if method is not None:
                model += "{}{}".format("u" if self.unrestricted else "", method)
            # Add a slash if we have both functional and basis set.
            if method is not None and basis_set is not None:
                model += "/"
            # And finally the basis set.
            if basis_set is not None:
                model += basis_set
                
            return model
        
        @property
        def route_section(self):
            """
            Get a Gaussian input file route section from this calculation target.
            """
            # Assemble our route line.
            # Add calc keywords.
            route_parts = [self.calculation_keywords]
            
            # Model chemistry
            route_parts.append(self.model_chemistry)
            
            # Solvent.
            if self.solvent is not None:
                #route_parts.append(self.keyword_to_string("SCRF", {"Solvent": self.solvent}))
                route_parts.append(str(Keyword("SCRF", {"Solvent": self.solvent})))
                
            # Empirical dispersion
            if self.DFT['empirical_dispersion'] is not None:
                route_parts.append(str(Keyword("EmpiricalDispersion", self.DFT['empirical_dispersion'])))
            
            # Finally, add any free-form options.
            for keyword_str in self.keywords:
                route_parts.append(str(Keyword(keyword_str, self.keywords[keyword_str])))
                    
            # Convert to string and return.
            return " ".join(route_parts)
        
        def prepare(self, output, input_coords):
            """
            Prepare this calculation for submission.
            
            :param output: Path to a directory to perform the calculation in.
            :param input_coords: A Silico_coords object containing the coordinates on which the calculation will be performed.
            """            
            if isinstance(input_coords, Calculation_directory_input):
                # Not supported ATM.
                raise NotImplementedError("Gaussian calculations cannot currently be prepared from directories")
            
            # Call parent.
            super().prepare(output, input_coords)
            
            # Decide on our file names.
            self.chk_file_name = self.safe_name(self.molecule_name + ".chk")
            self.rwf_file_name = self.safe_name(self.molecule_name + ".rwf")
            self.com_file_name = self.safe_name(self.molecule_name + ".com")
            
            # Get and load our com file template.
            self.com_file_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/gaussian/input_file.mako").render_unicode(calculation = self)
                
