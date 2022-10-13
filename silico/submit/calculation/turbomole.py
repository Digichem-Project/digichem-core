# General imports.
from mako.lookup import TemplateLookup
from pathlib import Path

# Silico imports.
from silico.submit.calculation.base import Concrete_calculation,\
    AI_calculation_mixin
from silico.config.configurable.option import Option
from silico.submit.memory import Memory
import silico
from silico.config.configurable.options import Options
from silico.input.directory import Calculation_directory_input
from silico.exception.configurable import Configurable_option_exception,\
    Configurable_exception
from silico.submit.translate import Basis_set
from silico.submit.calculation.base import validate_method as parent_validate_method
from silico.exception.base import Silico_exception

                
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
    performance = Options(
        memory = Option(type = Turbomole_memory,),
        parallel_mode = Option(help = "The type of parallelization to use. SMP uses shared memory and therefore is only suitable for parallelization across a single node, while MPI uses message-passing between processes and so can be used across multiple nodes. Use 'linear' to disable parallelization.", default = "SMP", choices = ("SMP", "MPI", "linear"), type = str)
    )
    modules = Option(help = "A list of turbomole commands/programs to execute in order.", list_type = tuple, type = str, required = True)
    

def validate_bse(option, owning_obj, value):
    """Function to check that a basis set exchange basis set has not been chosen (these are not yet supported for Turbomoles)."""
    if not option.default(owning_obj) == value:
        raise Configurable_option_exception(owning_obj, option, "BSE basis sets are not currently supported for Turbomole.")
    
    return True

def validate_method(option, owning_obj, value):
    """Function to check that appropriate RI options have been chosen."""
    # Call parent func first.
    if not parent_validate_method(option, owning_obj, value):
        return False
    
    # Panic if we're not using RI and we are using ADC(2) or CC2 (where it is mandatory).
    if not value['ri']['correlated']['calc'] and (\
            (value['mp']['calc'] and value['mp']['method'] in ["ADC(2)", "MP3", "MP4"]) or \
            value['cc']['calc']
        ):
        raise Configurable_option_exception(owning_obj, option, "The resolution of the identity approximation for correlated methods (method:ri:correlated:calc) must be turned on to use ADC(2), MP3, MP4, or any coupled cluster method.")
     
    # If SCS has been chosen, make sure we can use it.
    # TODO: The excited state methods CIS(D) and CIS(Dinf) can also use SCS, so this check is disabled for now.
#     if value['spin_component_scaling']['calc'] and (\
#             (not value['mp']['calc'] or not value['mp']['method'] in ["MP2", "ADC(2)"]) and 
#             (not value['cc']['calc'] or not value['cc']['method'] in ["CC2"]) and
#             (not value['properties']['es']['calc'] or not va)
#         ):
#         raise Configurable_option_exception(owning_obj, option, "Spin-component scaling can only be turned on for MP2, ADC(2) or CC2.")
        
    return True


class Turbomole_AI(Turbomole, AI_calculation_mixin):
    """
    Ab Initio calculations with Turbomole.
    """
    
    # Identifying handle.
    CLASS_HANDLE = ("Turbomole",)
    
    # Configurable options.
    redundant_internal_coordinates = Option(help = "Whether to use redundant internal coordinates", type = bool, default = True)
    #_define_options = Option(help = "Method keywords and options from the define general menu, including scf, mp2, cc etc. This option is currently unused.", type = dict, default = {})
    
    basis_set = Options(
        exchange = Option(validate = validate_bse)
    )
    
    performance = Options(
        # See issue #36 for more information on this.
        intel_AVX2_fix = Option(help = "Whether to disable the AVX2 CPU extensions by setting MKL_ENABLE_INSTRUCTIONS=SSE4_2 in the program environment. This option is required to fix a bug in some turbomole modules (ricc2 at least) which occurs on some intel CPUs. Setting this option when it is not required is likely to incur a performance penalty.", type = bool, default = False),
        define_timeout = Option(help = "The amount of time (s) to allow define to run for. After the given timeout, define will be forcibly stopped if it is still running, which normally occurs because something went wrong and define froze.", type = int, default = 60)
    )
    
    # General turbomole options.
    method = Options(validate = validate_method,
        dft = Options(
            dispersion = Option(choices = ("GD2", "GD3", "GD3BJ", "GD4", None)),
            dispersion_abc = Option(help = "Whether to switch on the three-body term.", type = bool, default = False)
        ),
        mp = Options(
            level = Option(help = "The MP order to use.", choices = ("MP2", "CIS(D)", "ADC(2)", "MP3", "MP4"), default = "MP2")
        ),
        cc = Options(
            level = Option(help = "The CC level to use.", type = str, choices = ("CCS", "CC2", "CCSD", "CCSD(T)"), default = "CC2")
        ),
        ri = Options(help = "Options for controlling various resolution of the identity (RI) approximations.",
            coulomb = Options(help = "Options for controlling the resolution of the identity approximation for calculating the electronic Coulomb interaction (RI-J).",
                calc = Option(help = "Whether to use the RI-J approximation (for HF or DFT methods).", type = bool, default = False),
                basis_set = Option(help = "The auxiliary basis set to use for the RI-J approximation. This option corresponds to the $jbas option. If not specified, an auxiliary basis set matching the main basis set will be chosen.", type = Basis_set, default = Basis_set("auto")),
            ),
            hartree_fock = Options(help = "Options for controlling the resolution of the identity approximation for calculating the Hartree-Fock exchange (RI-HF).",
                calc = Option(help = "Whether to use the RI-HF approximation.", type = bool, default = False),
                basis_set = Option(help = "The auxiliary basis set to use for the RI-HF approximation. This option corresponds to the $jkbas option. If not specified, an auxiliary basis set matching the main basis set will be chosen.", type = Basis_set, default = Basis_set("auto")),
            ),
            correlated = Options(help = "Options for controlling the resolution of the identity approximation for methods at the correlated second-order level (RI-C).",
                calc = Option(help = "Whether to use the RI-C approximation (for MP and CC methods). For MP2, this is optional, for all other MP and CC methods it is mandatory.", type = bool, default = False),
                basis_set = Option(help = "The auxiliary basis set to use for the RI-C approximation. This option corresponds to the $cbas option. If not specified, an auxiliary basis set matching the main basis set will be chosen.", type = Basis_set, default = Basis_set("auto")),
            ),
        ),
        scs = Options(help = "Options for controlling spin-component scaling (SCS). SCS is available for MP2, CIS(D), CIS(Dinf), ADC(2) and CC2.",
            calc = Option(help = "Whether to calculate using SCS.", type = bool, default = False),
            method = Option(help = "The SCS method to use, either SCS for conventional spin-component scaling, or SOS for spin-opposite scaling.", choices = ("SCS", "SOS"), default = "SCS"),
            opposite = Option(help = "The scaling factor for the opposite-spin component (COS).", type = float, default = None),
            same = Option(help = "The scaling factor for the same-spin component (CSS), only valid for the SCS method.", type = float, default = None),
        )
    )
    
    properties = Options(
        opt = Options(
            # TODO: These options will move somewhere better.
            ricc2 = Options(help = "Options specific to calculating optimised geometries with the ricc2 module.",
                optimise_symmetry = Option(help = "The symmetry of the state to optimise. Only meaningful for optimisations.", default = "a"),
                optimise_multiplicity = Option(help = "The multiplicity of the state to optimise. Only meaningful for optimisations. If not given, the ground state multiplicity will be used.", type = int, default = None)
            )
        ),
        es = Options(
            method = Option(choices = ("TD-HF", "CIS", "CIS(D)", "CIS(Dinf)", "TD-DFT", "TDA", "MP2", "ADC(2)", "CC2", "CCS"), default = "TD-HF"),
            multiplicity = Option(choices = ("Singlet", "Triplet"), default = "Singlet"),
            symmetry = Option(help = "The symmetry of the excited states to calculate.", default = "a"),
            # TODO: These options will probably move to a more general place in the future.
            ricc2 = Options(help = "Options specific to calculating excited states with the ricc2 module.",
                spectrum_operators = Option(help = "The operators with which to calculate oscillator strengths.", type = str, default = None),
                gradients = Option(help = "Whether to calculate excited state gradients.", type = bool, default = True)
            )
        )
    )
    
    scf = Options(
        # Turbomole does not seem to offer different SCF procedures.
        method = Option(choices = (None,)),
        density_convergence = Option(help = "The SCF density convergence criteria, expressed as 10^n (where n is the value of this option).", type = int, default = None),
        damping = Options(help = "Options for damping SCF oscillations, see the Turbomole manual for the option '$scfdamp'.",
            # SCF damping cannot be easily turned off in Turbomole.
            calc = Option(choices = (True,), type = bool, default = True),
            weight = Option(help = "Starting weight for the SCF damping procedure.", type = float, default = None),
            step = Option(help = "Amount to decrease weight by in each iteration.", type = float, default = None),
            min = Option(help = "The minimal damping weight.", type = float, default = None) 
        )
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
            help = "Options for $anadens data group (calculating difference density plots etc).",
            calculate = Option(help = "Whether to calculate $anadens.", type = bool, default = False),
            first_density = Option(help = "One of the two density files (.cao) to calculate from, relative to the calculation directory."),
            second_density = Option(help = "One of the two density files (.cao) to calculate from, relative to the calculation directory."),
            operator = Option(help = "The operator to use on the two density files, eg '-' (subtraction).", default = "-"),
            output = Option(help = "Name of the file to write to, relative to the calculation directory.")
        )
    )
    
    _modules = Option("modules", help = "Manually set the Turbomole modules to call. If not specified, these will be determined automatically.", type = str, list_type = list, default = ())
    
    @property
    def modules(self):
        """
        Get a list of modules (Turbomole sub-programs) to call.
        """
        # If given explicitly, use that.
        if len(self._modules) != 0:
            return self._modules
            
        # Otherwise, build them up.
        modules = []
        # If this is an opt, we use jobex to call our modules.
        if self.properties['opt']['calc']:
            # Build options for jobex.
            jobex_options = ["jobex -keep"]
            
            # Add max cycles.
            if self.properties['opt']['iterations'] is not None:
                jobex_options.append("-c " + str(self.properties['opt']['iterations']))
                
            # Add the method.
            if self.method['hf']['calc'] or self.method['dft']['calc']:
                jobex_options.append("-level scf")
                
            elif self.method['mp']['calc'] and self.method['mp']['level'] == "MP2":
                jobex_options.append("-level mp2")
                
            else:
                # ricc2 opt.
                jobex_options.append("-level cc2")
                
            # Add appropriate RI-J and RI-HF switches.
            # This option might not be necessary/do anything for ricc2.
            if self.method['ri']['coulomb']['calc'] or self.method['ri']['correlated']['calc']:
                jobex_options.append("-ri")
            
            # This option might only be meaningful for ricc2.
            if self.method['ri']['hartree_fock']['calc']:
                jobex_options.append("-rijk")
                
            # If we're doing opt and es, add that option.
            # This probably only works for TD-HF and TD-DFT
            if self.properties['es']['calc']:
                jobex_options.append("-ex")
                
            # Combine and set.
            modules.append(" ".join(jobex_options))
        
        else:
            # TOOD: What about gradients with grad or rdgrad?
            # All methods start with a HF or DFT run, using dscf or ridft.
            # dscf is used for 'conventional' calcs, those not using any RI-J or RI-HF approximations.
            if not self.method['ri']['coulomb']['calc'] and not self.method['ri']['hartree_fock']['calc']:
                modules.append("dscf")
            
            else:
                modules.append("ridft")
                
            # If we're calculating excited states with HF or DFT, call escf.
            if self.properties['es']['calc'] and (self.method['hf']['calc'] or self.method['dft']['calc']):
                modules.append("escf")
                
            # If we're using non-RI MP2, add mpgrad.
            if self.method['mp']['calc'] and self.method['mp']['level'] == "MP2" and not self.method['ri']['correlated']['calc']:
                modules.append("mpgrad")
                
            # If we're using RI-MP2, CSS, CIS(D), ADC(2) or CC2, add ricc2.
            # TODO: Need to check CCS is run with ricc2...
            elif (self.method['mp']['calc'] and self.method['mp']['level'] in ("MP2", "CIS(D)", "ADC(2)")) or (self.method['cc']['calc'] and self.method['cc']['level'] in ("CCS", "CC2")):
                modules.append("ricc2")
            
            # Finally, if we're using a higher order method, use ccsdf12.
            if (self.method['mp']['calc'] and self.method['mp']['level'] in ("MP3", "MP4")) or (self.method['cc']['calc'] and self.method['cc']['level'] in ("CCSD", "CCSD(T)")):
                modules.append("ccsdf12")
                
        # Finally, if we've been asked for frequencies, add aoforce or NumForce as appropriate.
        if self.properties['freq']['calc']:
            # Work out if we're doing numerical or analytical.
            if self.properties['freq']['numerical'] is not None:
                numerical = self.properties['freq']['numerical']
            
            # Not given explicitly, use analytical if we can, otherwise numerical.
            elif self.method['hf']['calc'] or self.method['dft']['calc']:
                numerical = False
            
            else:
                numerical = True
                
            # aoforce requires less setup than NumForce.
            if not numerical:
                modules.append("aoforce")
                
            else:
                # These options look similar to those for jobex, but are subtly different.
                numforce_options = ["NumForce"]
                
                # Add -ri flag if appropriate.
                if self.methods['ri']['coulomb']['calc'] or self.methods['ri']['correlated']['calc']:
                    numforce_options.append("-ri")
                
                # Choose level.
                if self.method['hf']['calc'] or self.method['dft']['calc']:
                    # No option.
                    pass
                
                elif self.method['mp']['calc'] and self.method['mp']['level'] == "mp2":
                    numforce_options.append("-level mp2")
                
                else:
                    # Assume ricc2 level.
                    numforce_options.append("-level cc2")
                    
                    # Add -rijk flag if appropriate.
                    if self.methods['ri']['hartree_fock']['calc']:
                        numforce_options.append("-rijk")
                        
                modules.append(" ".join(numforce_options))
        
        if len(modules) == 0:
            raise Configurable_exception(self, "The modules to run could not be automatically determined.")
        
        return modules
    
    @modules.setter
    def modules(self, value):
        self._modules = value
        
    @property
    def ricc2_method(self):
        """
        Return the ricc2 method to be used in this calculation, or None if ricc2 is not being used.
        """
        if self.method['mp']['calc'] and (self.method['mp']['level'] == "MP2" and self.method['ri']['correlated']['calc']):
            return "MP2"
        
        elif self.method['cc']['calc'] and self.method['cc']['level'] == "CC2":
            return "CC2"
        
        else:
            return None
        
    @property
    def post_HF_method(self):
        """
        Return the post HF method to be used in this calculation, or None if one is not being used.
        
        The method returned here is used to set the wavefunction in the ricc2 menu.
        """
        if self.method['mp']['calc']:
            return self.method['mp']['level']
        
        elif self.method['cc']['calc']:
            return self.method['cc']['level']
        
        else:
            return None
        
    @property
    def scs_line(self):
        """
        Get the SCS (spin-component scaling) input line for define.
        """
        line_parts = []
        # First, pick our method.
        line_parts.append(self.method['scs']['method'])
        
        # If we've got a specific COS factor, use that.
        if self.method['scs']['oppsoite'] is not None:
            line_parts.append("cos=" + self.method['scs']['opposite'])
            
        # And CSS factor.
        if self.method['scs']['method'] == "scs" and self.method['scs']['same'] is not None:
            line_parts.append("css=" + self.method['scs']['same'])
            
        return " ".join(line_parts)
    
    @property
    def optimise_multiplicity(self):
        """
        The multiplicity of the state to optimise.
        """
        if self.properties['opt']['ricc2']['optimise_multiplicity'] is not None:
            return self.properties['opt']['ricc2']['optimise_multiplicity']
        
        else:
            return self.multiplicity
    
    @property
    def exci_irrep_line(self):
        """
        Get the ricc2 excited states irrep line for define.
        """
        line_parts = ["irrep=" + self.properties['es']['symmetry']]
        
        # Add multiplicity if we're not unrestricted.
        if not self.is_unrestricted:
            line_parts.append("multiplicity=" + "1" if self.properties['es']['multiplicity'] == "Singlet" else "3")
           
        # Number of states.
        line_parts.append("nexc=" + str(self.properties['es']['num_states']))
        
        return " ".join(line_parts)
    
    @property
    def exci_gradient_line(self):
        """
        Get the ricc2 excited states gradients line for deifne.
        """
        irrep = self.properties['es']['symmetry']
        
        # Add mult if we are not unrestricted.
        if not self.is_unrestricted:
            irrep += "\{{}\}".format(self.properties['es']['multiplicity'])
            
        return "xgrad states=({} 1-{})".format(irrep, self.properties['es']['num_states'])
    
    
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
            return make_orbital_calc(name = self.meta['name'], memory = self.performance['memory'], num_cpu = self.performance['memory'], orbitals = orbitals, density = density, options = self.silico_options)
        
        def anadens_calc(self, first_density, second_density, file_name, operator = "-"):
            """
            Return a new calculation that can create cube files generated with the $anadens data group.
            """
            return make_anadens_calc(name = self.meta['name'], memory = self.performance['memory'], num_cpu = self.performance['memory'], first_density = first_density, second_density = second_density, file_name = file_name, operator = operator)


# TODO: This class inherits far too much from its parents, inheritance needs tweaking here.
class Turbomole_analysis(Turbomole_AI):
    """
    Calculations for performing analysis on an existing calculation.
    """
    
    # Identifying handle.
    CLASS_HANDLE = ("Turbomole-Analysis",)
    
    # No need to check for methods or properties.
    method = Options(validate = None)
    properties = Options(validate = None)
    
    analysis = Options(help = "Options that control the use of Turbomole in analysis mode, for example density plotting. Generally these options should not be modified for typical Turbomole calculations.",            
        plt = Options(
            help = "Options for orbital and density grid plotting.",
            calculate = Option(help = "Whether to plot densities.", type = bool, default = False),
            density = Option(help = "Whether to plot electron/spin density.", type = bool, default = False),
            orbitals = Option(help = "List of orbitals to plot for. Orbitals are identified by their 'irrep', a combination of symmetry and number.", type = tuple, default = []),
            format = Option(help = "The format to write to.", default = "cub", choices = ("cub", "plt", "map", "xyz", "plv"), type = str)
            ),
        anadens = Options(
            help = "Options for $anadens data group (calculating difference density plots etc).",
            calculate = Option(help = "Whether to calculate $anadens.", type = bool, default = False),
            first_density = Option(help = "One of the two density files (.cao) to calculate from, relative to the calculation directory."),
            second_density = Option(help = "One of the two density files (.cao) to calculate from, relative to the calculation directory."),
            operator = Option(help = "The operator to use on the two density files, eg '-' (subtraction).", default = "-"),
            output = Option(help = "Name of the file to write to, relative to the calculation directory.")
        )
    )


def make_orbital_calc(*, name, memory, num_cpu, orbitals = [], density = False, modules = None, options):
    """
    Get a calculation template that can be used to create orbital objects.
    
    :param name: The name of the calculation.
    :param memory: The amount of memory to use for the calculation.
    :param num_cpu: The number of CPUs to use for the calculation.
    :param orbitals: List of orbital indexes to create cubes for.
    :param density: Whether to create density cubes.
    :param modules: Turbomole modules to run.
    """
    if modules is None:
        modules = ["dscf"]
        
    # Append -proper flag to skip recalc.
    modules = ["{} -proper".format(module) for module in modules]
    
    # First generate our calculation.
    calc_t = Turbomole_analysis(
        meta = {"name": "Orbital Cubes for {}".format(name)},
        performance = {
            "memory": str(memory),
            "num_cpu": num_cpu,
        },
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
        post = {
            "write_summary": False,
            "write_report": False,
        },
        scratch_options = {
            "use_scratch": False
        }
    )
    
    # Prepare it for making new classes.
    calc_t.finalize()
    
    # Done.
    return calc_t


def make_anadens_calc(*, name, memory, num_cpu, first_density, second_density, file_name, operator = "-"):
    """
    Create a Turbomole calculation object that can be use to create difference density plots from an existing calculation (difference density plots are similar to NTOs plots but apply specifically to calculations performed with ricc2(?))
    
    :param name: A name to give the calculation.
    :param memory: The amount of memory to use for the calculation (note it's not clear if this option will be respected by ricc2 or not).
    :param num_cpu: The number of CPUs to use for the calculation (note it's not clear if this option will be respected by ricc2 or not).
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
    calc_t = Turbomole_analysis(
        meta = {'name': "Anadens for {}".format(name)},
        performance = {
            "memory": str(memory),
            "num_cpu": num_cpu,
        },
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
        post = {
            "write_summary": False,
            "write_report": False,
        },
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
    
    method = Options(
        mm = Options(help = "Options for molecular mechanics calculations (MM).",
            calc = Option(help = "Whether to use the MM method.", type = bool, default = True),
            force_field = Option(help = "The force field to use. Currently, only UFF is available.", choices = ("UFF", ), default = "UFF"),
            #maxcycle = Option(help = "Maximum number of UFF iterations. Set to 1 for a single point calculation.", type = int, default = 50),
            modus = Option(help = "If false, only topology will be be calculated.", type = bool, default = True),
            nqeq = Option(help = "The number of cycles to calculate partial charges.", type = int, default = 0),
            iterm = Option(help = "Switches controlling force field terms; please see the Turbomole manual for more details.", type = str, default = "111111"),
            econv = Option(help = "Energy convergence criteria.", type = str, default = "0.10D-07"),
            gconv = Option(help = "Gradient convergence criteria.", type = str, default = "0.10D-04"),
            #_qtot = Option("qtot", help = "The molecular charge. Leave blank to use the charge given in the input file.", default = None, type = float),
            dfac = Option(help = "Multiplication factor to determine bonds between atoms.", type = str, default = "1.10"),
            epssteep = Option(help = "Criteria for determining whether to perform a deepest-descent-step.", type = str, default = "0.10D+03"),
            epssearch = Option(help = "Criteria for performing a line-search step.", type = str, default = "0.10D-04"),
            dqmax = Option(help = "Maximum displacement of a coordinate (in a.u.).", type = str, default = "0.30"),
            mxls = Option(help = "Number of energy calculations for linesearch.", type = int, default = 25),
            dhls = Option(help = "Increment value for linesearch.", type = str, default = "0.10"),
            ahls = Option(help = "Start value for linesearch.", type = str, default = "0.00"),
            alpha = Option(help = "Alpha parameter; please see the Turbomole manual for more details.", type = str, default = "1.00"),
            beta = Option(help = "Beta parameter; please see the Turbomole manual for more details.", type = str, default = "0.00"),
            gamma = Option(help = "Gamma parameter; please see the Turbomole manual for more details.", type = str, default = "0.00"),
            transform = Option(help = "Whether to perform the transformation in the principle axis system.", type = bool, default = False),
            lnumhess = Option(help = "Whether to calculate a numerical Hessian", type = bool, default = False),
            lmd = Option(help = "Whether to perform an MD calculation.", type = bool, default = False),
        ),
        # TODO: Disable other options.
    )
    
    electron = Options(help = "Electron occupancy options.",
            charge = Option(help = "The molecular charge. Leave blank to use the charge given in the input file.", default = None, type = int)
        )
    
    @property
    def max_cyle(self):
        """The maximum number of UFF cycles."""
        if self.properties['sp']['calc']:
            return 1
        
        else:
            return self.properties['opt']['iterations']
    
    @property
    def qtot(self):
        """
        The molecule/system charge that we'll actually be using in the calculation.
        
        Unlike the _qtot attribute, this property will translate "auto" to the actual charge to be used.
        """
        # TODO: Should be implicit_charge?
        return self.electron['charge'] if self.electron['charge'] is not None else self.input_coords.charge
        
    
    # We only have one module to run.
    modules = ("uff",)
    