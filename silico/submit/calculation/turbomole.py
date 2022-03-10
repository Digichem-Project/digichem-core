# General imports.
from mako.lookup import TemplateLookup
from pathlib import Path

# Silico imports.
from silico.submit.calculation.base import Concrete_calculation,\
    AI_calculation_mixin
from silico.config.configurable.option import Option
from silico.submit.base import Memory
import silico
from silico.config.configurable.options import Options
from silico.file.input.directory import Calculation_directory_input

                
class Turbomole_memory(Memory):
    """
    A class for representing memory quantities in formats suitable for turbomole.
    """
    
    # When outputting units, whether to separate the number and unit with a space.
    SPACE_UNIT = False

class Turbomole(Concrete_calculation):
    """
    Abstract top-class for turbomole calcs.
    """
    
    # Identifying handle.
    CLASS_HANDLE = ()
    
    # A list of strings describing the expected input file types (file extensions) for calculations of this class. The first item of this list will be passed to obabel via the -o flag. 
    INPUT_FILE_TYPES = ["tmol"]
    
    # The format of the output file containing coordinates.
    OUTPUT_COORD_TYPE = "tmol"
    
    # Configurable options.
    memory = Option(help = "The amount of memory to use for the calculation.", required = True, type = Turbomole_memory, dump_func = lambda option, configurable, value: str(value))
    modules = Option(help = "A list of turbomole commands/programs to execute in order.", type = tuple, required = True)
    parallel_mode = Option(help = "The type of parallelization to use. SMP uses shared memory and therefore is only suitable for parallelization across a single node, while MPI uses message-passing between processes and so can be used across multiple nodes. Use 'linear' to disable parallelization.", default = "SMP", choices = ("SMP", "MPI", "linear"), type = str)
    

class Turbomole_AI(Turbomole, AI_calculation_mixin):
    """
    Ab Initio calculations with Turbomole.
    """
    
    # Identifying handle.
    CLASS_HANDLE = ("Turbomole",)
    
    # Configurable options.
    basis_set = Option(help = "The basis set to use.", required = False, type = str)
    force_unrestricted = Option(help = "Whether to force use of unrestricted HF. This option only has an effect if multiplicity is 1; as all other multiplicities will use unrestricted HF by necessity.", type = bool, default = False)
    redundant_internal_coordinates = Option(help = "Whether to use redundant internal coordinates", type = bool, default = True)
    methods = Option(help = "Method keywords and options from the define general menu, including scf, mp2, cc etc.", type = dict, default = {})
    define_timeout = Option(help = "The amount of time (s) to allow define to run for. After the given timeout, define will be forcibly stopped if it is still running, which normally occurs because something went wrong and define froze.", type = int, default = 60)
    
    @property
    def unrestricted_HF(self):
        """
        Whether this calculation will be using unrestricted HF.
        """
        return self.multiplicity != 1 or self.force_unrestricted
    
    # General turbomole options.
    scf = Options(
        help = "Options for SCF.",
        iter = Option(help = "Maximum number of SCF iterations. Set to None for Turbomole defaults.", type = int, default = None),
        scfconv = Option(help = "SCF convergence criteria will be set to 1.0x10^-n where n is the value of this option (at least between 4-12).", type = int, default = None),
        denconv = Option(help = "SCF density convergence criteria will be set to 1.0x10^-n where n is the value of this option.", type = int, default = None),
        )
    scf_damp = Options(
        help = "Options for SCF damping which can help reduce SCF oscillations, please see the turbomole manual for the option '$scfdamp'.",
        start = Option(help = "Starting weight for SCF damping procedure.", type = str, default = None),
        step = Option(help = "Amount to decrease weight by in each iteration.", type = str, default = None),
        min = Option(help = "The minimal damping weight.", type = str, default = None)
        )
    dft = Options(
        help = "Options for DFT.",
        functional = Option(help = "The DFT functional to use. If None is given, DFT will not be used.", type = str, default = None),
        grid = Option(help = "The grid size to use.", type = str, default = None),
        )
    dft_exci = Options(
        help = "Options for calculation of excited states with DFT (TDA or TD-DFT)",
        symmetry = Option(help = "Symmetry of the excited states to calculate.", type = str, default = "a"),
        multiplicity = Option(help = "Multiplicity of the excited states to calculate.", type = int, default = 1, choices = (None, 1, 3)),
        nexc = Option(help = "The number of excited states to calculate.", type = int, default = 0),
        TDA = Option(help = "Whether to use the Tamm–Dancoff approximation.", type = bool, default = False)
        )
    optimisation_state = Option(help = "The excited state to optimise. This should not exceed the number of excited states to calculate", type = int, default = None)
    dispersion = Options(
        help = "Options for dispersion correction.",
        dsp = Option(help = "Dispersion correction to use.", choices = ("GD2", "GD3", "GD3BJ", "GD4", None), default = None),
        abc = Option(help = "Whether to switch on the three-body term. This option is ignored if no dsp is to be used.", type = bool, default = False)
        )
    cc = Options(
        help = "Options for ricc2 (approximate coupled-cluster calculations).",
        cbas = Option(help = "The optional cbas auxiliary basis set to use. Specify 'auto' to use the same as the main basis set.", type = str, default = "auto"),
        cabs = Option(help = "The optional cabs complementary auxiliary basis set to use. Specify 'auto' to use the same as the main basis set.", type = str, default = None),
        jkbas = Option(help = "The optional jkbas basis set to use for Fock Matrices. Specify 'auto' to use the same as the main basis set.", type = str, default = None),
        )
    ricc2 = Options(
        help = "Options for ricc2.",
        model = Option(help = "The computational model to use (CC2, MP2 etc.). If no model is specified, no cc options will be used.", choices = ("ccs" , "cis", "mp2", "cid(d)", "adc(2)", "cc2", "ccsd", "ccsd(t)", "mp4", None), default = None),
        scs = Option(help = "Whether to use spin-component scaling", choices = ("scs", "sos", None), default = None),
        # See issue #36 for more information on this.
        intel_AVX2_fix = Option(help = "Whether to disable AVX2 CPU extensions by setting MKL_ENABLE_INSTRUCTIONS=SSE4_2 in the program environment. This option is required to fix a bug in some turbomole modules (ricc2 at least) which occurs on some intel CPUs. Setting this option when it is not required is likely to incur a performance penalty.", type = bool, default = False)
        )
    cc_geoopt = Options(
        help = "Options for geometry optimisations.",
        wavefunction = Option(help = "Wavefunction (model) to optimise (CC2, MP2 etc.)", type = str, default = None),
        state = Option(help = "Molecular state to optimise, eg. (a 1). The option may be left blank to use Turbomole defaults", type = str, default = None)
        )
    ricc2_exci = Options(
        help = "Options for calculations of excited states with ricc2.",
        symmetry = Option(help = "Symmetry of the excited states to calculate.", type = str, default = "a"),
        multiplicity = Option(help = "Multiplicity of the excited states to calculate.", type = int, default = 1),
        nexc = Option(help = "The number of excited states to calculate.", type = int, default = 0),
        oscillator_strengths = Option(help = "The operators with which to calculate oscillator strengths.", type = str, default = None),
        gradients = Option(help = "Whether to calculate excited state gradients.", type = bool, default = True)
        )
    analysis = Options(help = "Options that control the use of Turbomole in analysis mode, for example density plotting. Generally these options should not be modified for typical Turbomole calculations.",            
        plt = Options(
            help = "Options for orbital and density grid plotting.",
            calculate = Option(help = "Whether to plot densities.", type = bool, default = False),
            density = Option(help = "Whether to plot electron/spin density.", type = bool, default = False),
            orbitals = Option(help = "List of orbitals to plot for. Orbitals are identified by their 'irrep', a combination of symmetry and number.", type = tuple, default = []),
            format = Option(help = "The format to write to.", default = "cub", choices = ("cub", "plt", "map", "xyz", "plv"), type = str)
            ),
        anadens = Options(
            help = "Options for $anadens data group (calculating differential density plots etc).",
            calculate = Option(help = "Whether to calculate $anadens.", type = bool, default = False),
            first_density = Option(help = "One of the two density files (.cao) to calculate from, relative to the calculation directory."),
            second_density = Option(help = "One of the two density files (.cao) to calculate from, relative to the calculation directory."),
            operator = Option(help = "The operator to use on the two density files, eg '-' (subtraction).", default = "-"),
            output = Option(help = "Name of the file to write to, relative to the calculation directory.")
        )
    )
    
    @property
    def pretty_functional(self):
        """
        The name of the functional to use in a pretty format.
        An empty string if no functional is to be used.
        """
        if self.dft['functional'] is None:
            return ""
        else:
            return self.dft['functional'].upper()
        
    @property
    def func(self):
        """
        The name of the functional to use for this calculation in a format acceptable to turbomole.
        
        If no functional is to be used (because this is not a DFT calculation), then None is returned.
        """
        func = self.dft['functional']
        
        # If none, giveup.
        if func is None:
            return None
        
        # Convert to lowercase, but the G in Gaussian is upper case (sigh).
        return func.lower().replace("gaussian", "Gaussian")
        
    
    @property
    def unpaired_electrons(self):
        """
        The number of unpaired electrons at this given multiplicity.
        """
        return self.multiplicity -1
    
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Concrete_calculation._actual):
        """
        Inner class for calculations.
        """
        
        def __init__(self, *args, **kwargs):
            """
            Constructor for calculation objects.
            
            :param output: Path to a directory to perform the calculation in.
            :param input_coords: A string containing input coordinates in a format suitable for this calculation type.
            :param molecule_name: A name that refers to the system under study (eg, Benzene etc).
            """
            super().__init__(*args, **kwargs)
    
            self.define_input = None
            
        def prepare(self, output, input_coords, *args, **kwargs):
            """
            Prepare this calculation for submission.
            
            :param output: Path to a directory to perform the calculation in.
            :param input_coords: A string containing input coordinates in a format suitable for this calculation type.
            :param molecule_name: A name that refers to the system under study (eg, Benzene etc).
            """
            super().prepare(output, input_coords, *args, **kwargs)
            
            # Get and load our define template.
            # The type of template to use depends on whether we are restarting a previous calculation or not.
            template = "/submit/turbomole/define/restart.mako" if isinstance(input_coords, Calculation_directory_input) else "/submit/turbomole/define/driver.mako"
            self.define_input = TemplateLookup(directories = str(silico.default_template_directory())).get_template(template).render_unicode(calculation = self)
                
        def orbital_calc(self, orbitals, density):
            """
            Return a new calculation that can create orbital cube files from this calculation.
            
            :param orbitals: A list of orbital irreps to make cubes for.
            """
            return make_orbital_calc(name = self.name, memory = self.memory, num_CPUs = self._num_CPUs, orbitals = orbitals, density = density, options = self.silico_options)
        
        def anadens_calc(self, first_density, second_density, file_name, operator = "-"):
            """
            Return a new calculation that can create cube files generated with the $anadens data group.
            """
            return make_anadens_calc(name = self.name, memory = self.memory, num_CPUs = self._num_CPUs, first_density = first_density, second_density = second_density, file_name = file_name, operator = operator)
            
    
def make_orbital_calc(*, name, memory, num_CPUs, orbitals = [], density = False, modules = None, options):
    """
    Get a calculation template that can be used to create orbital objects.
    
    :param name: The name of the calculation.
    :param memory: The amount of memory to use for the calculation.
    :param num_CPUs: The number of CPUs to use for the calculation.
    :param orbitals: List of orbital indexes to create cubes for.
    :param density: Whether to create density cubes.
    :param modules: Turbomole modules to run.
    """
    if modules is None:
        modules = ["dscf"]
        
    # Append -proper flag to skip recalc.
    modules = ["{} -proper".format(module) for module in modules]
    
    # First generate our calculation.
    calc_t = Turbomole_AI(
        name = "Orbital Cubes for {}".format(name),
        #"programs": [program_name],
        memory = str(memory),
        num_CPUs = num_CPUs,
        analysis = {
            'plt': {
                "calculate": True,
                "density": density,
                "orbitals": orbitals,
                "format": "cub"
            },
        },
        modules = modules,
        # We don't need these.
        write_summary = False,
        write_report = False,
        scratch_options = {
            "use_scratch": False
        }
    )
    
    # Prepare it for making new classes.
    calc_t.finalize()
    
    # Done.
    return calc_t


def make_anadens_calc(*, name, memory, num_CPUs, first_density, second_density, file_name, operator = "-"):
    """
    Create a Turbomole calculation object that can be use to create differential density plots from an existing calculation (differential density plots are similar to NTOs plots but apply specifically to calculations performed with ricc2(?))
    
    :param name: A name to give the calculation.
    :param memory: The amount of memory to use for the calculation (note it's not clear if this option will be respected by ricc2 or not).
    :param num_CPUs: The number of CPUs to use for the calculation (note it's not clear if this option will be respected by ricc2 or not).
    :param first_density: The name of one of the two density files (.cao) to calculate the difference from.
    :param second_density: The name of the other of the two density files (.cao) to calculate the difference from.
    :param file_name: The name of the file to write to. Like the density files, this should be relative to the calc dir (and be careful of weird characters...). Because Turbomole forces this extension, the file must end in .cub
    :param operator: The operator to use on the two densities, (try) and see the Turbomole manual for allowed options. Subtraction ('-') is allowed at the very least, possibly also addition ('+'). The calculated density will be calculated according to = first_density operator second_density.
    """
    file_name = Path(file_name)
    # Panic if our file has the wrong extension.
    if file_name.suffix != ".cub":
        raise ValueError("The file extension/suffix of file_name must be .cub")
    
    # First generate our calculation.
    calc_t = Turbomole_AI(
        name = "Anadens for {}".format(name),
        memory = str(memory),
        num_CPUs = num_CPUs,
        analysis = {
            'plt': {
                "calculate": True,
                "format": "cub"
            },
            'anadens': {
                "calculate": True,
                "first_density": first_density,
                "second_density": second_density,
                "operator": operator,
                "output": file_name.stem
            }
        },
        modules = ['ricc2 -fanal'],
        # We don't need these.
        write_summary = False,
        write_report = False,
        scratch_options = {
            "use_scratch": False
        }
    )
    
    # Prepare it for making new classes.
    calc_t.finalize()
    
    # Done.
    return calc_t


class Turbomole_UFF(Turbomole):
    """
    Universal force-field calculations with Turbomole.
    """
    
    # Identifying handle.
    CLASS_HANDLE = ("Turbomole-UFF",)
    
    # A list of strings describing the expected input file types (file extensions) for calculations of this class. The first item of this list will be passed to obabel via the -o flag. 
    INPUT_FILE_TYPES = ["tmol"]
    
    # The format of the output file containing coordinates.
    OUTPUT_COORD_TYPE = "tmol"
    
    # Configurable options.
    maxcycle = Option(help = "Maximum number of UFF iterations. Set to 1 for a single point calculation.", type = int, default = 50)
    modus = Option(help = "If false, only topology will be be calculated.", type = bool, default = True)
    nqeq = Option(help = "The number of cycles to calculate partial charges.", type = int, default = 0)
    iterm = Option(help = "Switches controlling force field terms; please see the Turbomole manual for more details.", type = str, default = "111111")
    econv = Option(help = "Energy convergence criteria.", type = str, default = "0.10D-07")
    gconv = Option(help = "Gradient convergence criteria.", type = str, default = "0.10D-04")
    _qtot = Option("qtot", help = "The molecular charge. Leave blank to use the charge given in the input file.", default = None, type = float)
    dfac = Option(help = "Multiplication factor to determine bonds between atoms.", type = str, default = "1.10")
    epssteep = Option(help = "Criteria for determining whether to perform a deepest-descent-step.", type = str, default = "0.10D+03")
    epssearch = Option(help = "Criteria for performing a line-search step.", type = str, default = "0.10D-04")
    dqmax = Option(help = "Maximum displacement of a coordinate (in a.u.).", type = str, default = "0.30")
    mxls = Option(help = "Number of energy calculations for linesearch.", type = int, default = 25)
    dhls = Option(help = "Increment value for linesearch.", type = str, default = "0.10")
    ahls = Option(help = "Start value for linesearch.", type = str, default = "0.00")
    alpha = Option(help = "Alpha parameter; please see the Turbomole manual for more details.", type = str, default = "1.00")
    beta = Option(help = "Beta parameter; please see the Turbomole manual for more details.", type = str, default = "0.00")
    gamma = Option(help = "Gamma parameter; please see the Turbomole manual for more details.", type = str, default = "0.00")
    transform = Option(help = "Whether to perform the transformation in the principle axis system.", type = bool, default = False)
    lnumhess = Option(help = "Whether to calculate a numerical Hessian", type = bool, default = False)
    lmd = Option(help = "Whether to perform an MD calculation.", type = bool, default = False)
    
    
    @property
    def qtot(self):
        """
        The molecule/system charge that we'll actually be using in the calculation.
        
        Unlike the _qtot attribute, this property will translate "auto" to the actual charge to be used.
        """
        return self._qtot if self._qtot is not None else self.input_coords.charge
        
    
    # We only have one module to run.
    modules = ("uff",)
    