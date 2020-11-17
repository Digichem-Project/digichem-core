from silico.submit.calculation.base import Concrete_calculation
from silico.config.configurable.option import Option
from silico.submit.base import Memory
from mako.lookup import TemplateLookup
import silico
from silico.config.configurable.options import Options

				
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
	memory = Option(help = "The amount of memory to use for the calculation.", required = True, type = Turbomole_memory, rawtype = str)
	modules = Option(help = "A list of turbomole commands/programs to execute in order.", type = tuple, required = True)
	parallel_mode = Option(help = "The type of parallelization to use. SMP uses shared memory and therefore is only suitable for parallelization across a single node, while MPI uses message-passing between processes and so can be used across multiple nodes. Use 'linear' to disable parallelization.", default = "SMP", choices = ("SMP", "MPI", "linear"), type = str)
	
class Turbomole_AI(Turbomole):
	"""
	Ab Initio calculations with Turbomole.
	"""
	
	# Identifying handle.
	CLASS_HANDLE = ("Turbomole",)
	
	# Configurable options.
	basis_set = Option(help = "The basis set to use.", required = True, type = str)
	_charge = Option("charge", help = "Forcibly set the molecule charge. Leave blank to use the charge given in the input file.", default = None, type = float)
	_multiplicity = Option("multiplicity", help = "Forcibly set the molecule multiplicity. Leave blank to use the multiplicity given in the input file. A value of 1 will request turbomole defaults, which will be RHF singlet in most cases. Any other multiplicity will request UHF.", default = None, type = int)
	force_unrestricted = Option(help = "Whether to force use of unrestricted HF. This option only has an effect if multiplicity is 1; as all other multiplicities will use unrestricted HF by necessity.", type = bool, default = False)
	redundant_internal_coordinates = Option(help = "Whether to use redundant internal coordinates", type = bool, default = True)
	methods = Option(help = "Method keywords and options from the define general menu, including scf, mp2, cc etc.", type = dict, default = {})
	define_timeout = Option(help = "The amount of time (s) to allow define to run for. After the given timeout, define will be forcibly stopped if it is still running, which normally occurs because something went wrong and define froze.", type = int, default = 15)	
	
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
		model = Option(help = "The computational model to use (CC2, MP2 etc.). If no model is specified, no cc options will be used.", type = str, default = None),
		scs = Option(help = "Whether to use spin-component scaling", choices = ("scs", "sos", None), default = None)
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
			
			# If we've been asked to write everything to scratch, print a warning.
			#if self.scratch_options.all_output:
			#	logging.getLogger(silico.logger_name).warning("'{}';'scratch_options: all_output: True' is not supported for Turbomole; this setting will be ignored".format(self.NAME))
			
		def prepare(self, *args, **kwargs):
			"""
			Prepare this calculation for submission.
			
			:param output: Path to a directory to perform the calculation in.
			:param input_coords: A string containing input coordinates in a format suitable for this calculation type.
			:param molecule_name: A name that refers to the system under study (eg, Benzene etc).
			"""
			super().prepare(*args, **kwargs)
			
			# Get and load our define template.
			self.define_input = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/define/driver.mako").render_unicode(calculation = self)
				
	
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
	
# 	############################
# 	# Class creation mechanism #
# 	############################
# 	
# 	class _actual(Turbomole_UFF._actual):
# 		"""
# 		Inner class for calculations.
# 		"""
		