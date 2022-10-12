from mako.lookup import TemplateLookup
import deepmerge

import silico
from silico.config.configurable.option import Option
from silico.submit.calculation import Concrete_calculation
from silico.config.configurable.options import Options
from silico.input.directory import Calculation_directory_input
from silico.submit.calculation.base import AI_calculation_mixin
from silico.input.silico import Silico_coords


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
            self.add_option(option)
                
    def add_option(self, option):
        """
        Add an option to this keyword.
        
        Options can either by a dict, in which case they represent an option of the form keyword=(option=value)
        Or a non dict, in which case they represent options of the form: keyword=(option)
        """
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
    # Gaussian allows the exact CPUs to be specifed in a list.
    # NOTE: This list form is generally poorly tested...
    performance = Options(
        cpu_list = Option(help = "A list of integers specifying specific CPUs to use for the calculation, starting at 0. cpu_list and num_cpu are mutually exclusive.", exclude = ("num_cpu",), default = (), type = int, list_type = tuple),
        num_cpu = Option(help = "An integer specifying the number of CPUs to use for this calculation. cou_list and num_cpu are mutually exclusive.", exclude = ("cpu_list",))
    )
    
    method = Options(
        dft = Options(
            # Add the specific grid options Gaussian supports.
            grid = Option(choices = ("Coarse", "SG1", "Fine", "Ultrafine", "Superfine", None), default = "Ultrafine"),
            # And dispersion.
            dispersion = Option(choices = ("PFD", "GD2", "GD3", "GD3BJ", None))
        )
    )
    
    scf = Options(
        method = Option(choices = (None, "DM", "SD", "SSD", "QC", "XQC", "YQC")),
        damping = Options(
            iterations = Option(help = "Allow SCF damping for the first N cycles (where N is the value of this option).", type = int, default = None)
        )
    )
    
    properties = Options(
        opt = Options(
            options = Option(help = "Additional options to specify.", type = dict, default = {})
        ),
        freq = Options(
            options = Option(help = "Additional options to specify.", type = dict, default = {})
        ),
        es = Options(
            options = Option(help = "Additional options to specify.", type = dict, default = {})
        )
    )
    
    solution = Options(
        # Add Gaussian specific solvent methods.
        model = Option(choices = ("PCM", "CPCM", "Dipole", "IPCM", "SCIPCM", "SMD"), default = "PCM"),
    )
    
    post = Options(
        convert_chk = Option(help = "Whether to create an .fchk file at the end of the calculation", default = True, type = bool),
        keep_chk = Option(help = "Whether to keep the .chk file at the end of the calculation. If False, the .chk file will be automatically deleted, but not before it is converted to an .fchk file (if convert_chk is True)", default = False, type = bool),
        keep_rwf = Option(help = "Whether to keep the .rwf file at the end of the calculation. If False, the .rwf file will be automatically deleted", default = False, type = bool)
    )
    
    keywords = Option(help = "Additional Gaussian route keywords. These are written to the input file with only minor modification ('keyword: option' becomes 'keyword=(option)' and 'keyword: {option: value}' becomes 'keyword=(option=value)'), so any option valid to Gaussian can be given here", default = {'Population': 'Regular', 'Density': 'Current'}, type = dict)            
    
    @property
    def method_keyword(self):
        """
        A string describing the method to use (either a DFT functional or a post-HF method).
        """
        if self.method['hf']['calc']:
            return "HF"
        
        elif self.method['dft']['calc']:
            return self.method['dft']['functional'].to_gaussian()
        
        elif self.method['mp']['calc']:
            return self.method['mp']['level']
        
        elif self.method['cc']['calc']:
            # TODO: Need to check what happens when EOMCCSD is given and a CC method (probably nothing good?).
            return self.method['cc']['level']
        
        else:
            # Should be unreachable.
            raise NotImplementedError("Unrecognised method?")
    
    @property
    def calculation_keywords(self):
        """
        Get a string containing the Gaussian calculation keywords and associated options for this calculation.
        """
        keywords = []
        
        # Single point.
        if self.properties['sp']['calc']:
            keywords.append(Keyword("SP"))
        
        # Optimisations.
        if self.properties['opt']['calc']:
            # Make an opt keyword.
            opt_keyword = Keyword("Opt", self.optimisation['options'])
            
            # Add max number of steps if given.
            if self.properties['opt']['iterations'] is not None:
                opt_keyword.add_option({"MaxCycles": self.properties['opt']['iterations']})
            
            keywords.append(opt_keyword)
            
        # Frequencies.
        if self.properties['freq']['calc']:
            keywords.append(Keyword("Freq", self.properties['freq']['options']))
            
        # Excited states.
        if self.properties['es']['calc']:
            # First, build our options.
            options = {
                self.properties['es']['multiplicity']: "",
                'nstates': self.properties['es']['num_states'],
                'root': self.properties['es']['state_of_interest']
            }
            
            # Decide which method we're using to calculate excited states.
            es_method = "TD" if self.properties['es']['method'] == "TD-DFT" else self.properties['es']['method']
            
            # Add our keyword, which changes base on whether we're using TDA or TD-DFT.
            keywords.append(Keyword(es_method, options, self.DFT_excited_states['options']))
            
        # Merge and return.
        return " ".join([str(keyword) for keyword in keywords])
    
    
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
                basis_set = self.basis_set['internal'].to_gaussian()
                
            elif len(self.basis_set['exchange']) > 0:
                # Our basis set label will be either gen or genECP (depending on whether we have any ECPs).
                basis_set = "gen" if not self.basis_set['exchange'].has_ECPs(self.input_coords.elements) else "genECP"
                
            else:
                basis_set = None
            
            model = ""
            method = self.method_keyword
            # Add the functional.
            if method is not None:
                model += "{}{}".format("u" if self.electron['unrestricted'] else "", method)
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
            
            # SCF options.
            scf_key = Keyword("SCF")
            # If the SCF method is non default, add that.
            if self.scf['method'] is not None:
                scf_key.add_option(self.scf['method'])
                
            # Set max iterations.
            if self.scf['iterations'] is not None:
                scf_key.add_option({"MaxCycle": self.scf['iterations']})
                
            # Set convergence.
            if self.scf['convergence'] is not None:
                scf_key.add_option({"Conver": self.scf['convergence']})
                
            # Damping options.
            if self.scf['damping']['calc']:
                scf_key.add_option("Damp")
                
                if self.scf['damping']['iterations'] is not None:
                    scf_key.add_option({"NDamp": self.scf['damping']['iterations']})
                    
            # Add SCF if there's at least one option.
            if len(scf_key.options) > 0:
                route_parts.append(str(scf_key))
            
            # Solvent.
            if self.solution['calc']:
                route_parts.append(str(Keyword("SCRF", {self.solution['model']: "" ,"Solvent": self.solution['solvent']})))
            
            # Additional DFT options.
            if self.method['dft']['calc']:
                # Empirical dispersion.
                if self.method['dft']['dispersion'] is not None:
                    route_parts.append(str(Keyword("EmpiricalDispersion", self.method['dft']['dispersion'])))
                
                # Integration grid.
                if self.method['dft']['grid'] is not None:
                    route_parts.append(str(Keyword("Integral", self.method['dft']['grid'])))
            
            # Finally, add any free-form options.
            for keyword_str in self.keywords:
                options = (self.keywords[keyword_str],) if not isinstance(self.keywords[keyword_str], list) and not isinstance(self.keywords[keyword_str], tuple) else self.keywords[keyword_str]
                route_parts.append(str(Keyword(keyword_str, *options)))
                    
            # Convert to string and return.
            return " ".join(route_parts)
        
        def prepare(self, output, input_coords, *args, **kwargs):
            """
            Prepare this calculation for submission.
            
            :param output: Path to a directory to perform the calculation in.
            :param input_coords: A Silico_coords object containing the coordinates on which the calculation will be performed.
            """            
            if isinstance(input_coords, Calculation_directory_input):
                # Not supported ATM.
                raise NotImplementedError("Gaussian calculations cannot currently be prepared from directories")
            
            # Call parent.
            super().prepare(output, input_coords, *args, **kwargs)
            
            # Decide on our file names.
            self.chk_file_name = self.safe_name(self.molecule_name + ".chk")
            self.rwf_file_name = self.safe_name(self.molecule_name + ".rwf")
            self.com_file_name = self.safe_name(self.molecule_name + ".com")
            
            # Get and load our com file template.
            self.com_file_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/gaussian/input_file.mako").render_unicode(calculation = self, write_geom = isinstance(input_coords, Silico_coords))
        
        def NTO_calc(self, transition):
            """
            Return a new calculation that can create a chk file containing NTOs for a given transition.
            """
            return make_NTO_calc(name = self.meta['name'], memory = self.performance['memory'], num_cpu = self._num_cpu, transition = transition, scratch_options = self.scratch_options)


def make_NTO_calc(*, name, memory, num_cpu, transition, scratch_options = None, scratch_path = None):
    """
    Create a Gaussian calculation object that can be use to create natural transition orbitals from an existing calculation.
    
    :param name: A name to give the calculation.
    :param memory: The amount of memory to use for the calculation (note it's not clear if this option will be respected by ricc2 or not).
    :param num_cpu: The number of CPUs to use for the calculation (note it's not clear if this option will be respected by ricc2 or not).
    :param transition: The integer number of the transition to calculate NTOs for.
    """
    if scratch_options is None and scratch_path is None:
        raise ValueError("One of either scratch_options or scratch_path must be given")
    
    # Setup scratch options if not given.
    if scratch_options is None:
        scratch_options = {
            # Gaussian is basically broken without use of a scratch directory.
            # When not given, Gaussian appears to use the current directory as the scratch,
            # but gaussian can't handle whitespace in the scratch path, which is hard to control
            # if the cwd is used (which often has whitespace in it).
            #
            # Hence we have to turn scratch on.
            # However, we can't use the default scratch location (because it might not exist).
            # Likewise, we still can't use any whitespace in the path dir, so we have to be 
            # careful about where we choose.
            #
            # The current solution is to force the user to give a location, which is fine
            # but doesn't feel very satisfactory.
            "use_scratch": True,
            "path": scratch_path
            
        }
    
    calc_t = Gaussian(
        meta = {'name': "NTOs for {}".format(name)},
        performance = {
            "memory": str(memory),
            "num_cpu": num_cpu,
        },
        keywords = {
            "Geom": "AllCheck",
            "Guess": ("Read", "Only"),
            "Density": ("Check", {"Transition": transition}),
            "Population": ("Minimal", "NTO", "SaveNTO")
        },
        # We don't need these.
        post = {
            "write_summary": False,
            "write_report": False,
            "convert_chk": False,
            "keep_chk": True,
        },
        scratch_options = scratch_options
    )
    
    # Prepare it for making new classes.
    calc_t.finalize()
    
    # Done.
    return calc_t
    
    
    