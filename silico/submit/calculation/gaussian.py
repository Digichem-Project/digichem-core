from mako.lookup import TemplateLookup

import silico
from silico.exception import Configurable_exception
from silico.exception.base import Submission_error, Silico_exception
from silico.config.configurable.option import Option
from silico.misc.base import is_int
from silico.submit.calculation import Concrete_calculation
from logging import getLogger

class Gaussian_DFT(Concrete_calculation):
	"""
	DFT (density functional theory) calculations with Gaussian.
	"""
	# Identifying handle.
	CLASS_HANDLE = ("Gaussian-DFT",)
	
	# A list of strings describing the expected input file types (file extensions) for calculation's of this class. The first item of this list will be passed to obabel via the -o flag. 
	INPUT_FILE_TYPES = ["gau", "com", "gjf", "gjc"]

	# Configurable options.
	CPU_list = Option(help = "A list of integers specifying specific CPUs to use for the calculation, starting at 0. CPU_list and num_CPUs are mutually exclusive", exclude = ("num_CPUs",), default = (), type = tuple)
	num_CPUs = Option(help = "An integer specifying the number of CPUs to use for this calculation. CPU_list and num_CPUs are mutually exclusive", exclude = ("CPU_list",), default = 1, type = int)
	calculation_keywords = Option(help = "A list of Gaussian keywords specifying the calculation(s) to perform. Options, where appropriate, can also be included (eg, TDA=(NStates=10) )", default = (), type = tuple)
	functional = Option(help = "The DFT functional to use", type = str)
	basis_set = Option(help = "The basis set to use. 'Gen' or 'GenECP' should be given here if an external basis set is to be used", type = str)
	_external_basis_sets = Option(
		"external_basis_sets",
		help = "A list of external basis sets to use. The order given here is the order the basis sets will be appended to the input file",
		choices = lambda option, configurable: [name for basis_set in configurable.available_basis_sets for name in basis_set.NAMES],
		type = tuple,
		default = ()
	)
	_external_ECPs = Option(
		"external_ECPs",
		help = "A list of external ECPs (effective core potentials) to use",
		choices = lambda option, configurable: [name for basis_set in configurable.external_ECPs for name in basis_set.NAMES],
		type = tuple,
		default = ()
	)
	_multiplicity = Option("multiplicity", help = "Forcibly set the system multiplicity. Use 'auto' to use the multiplicity given in the input file", default = 'auto', validate = lambda option, configurable, value: value == "auto" or is_int(value))
	_charge = Option("charge", help = "Forcibly set the system charge. Use 'auto' to use the charge given in the input file", default = 'auto', validate = lambda option, configurable, value: value == "auto" or is_int(value))
	solvent = Option(help = "Name of the solvent to use for the calculation (the model used is SCRF-PCM)", default = None, type = str)
	options = Option(help = "Additional Gaussian route options. These are written to the input file with only minor modification ('name : value' becomes 'name=value'), so any option valid to Gaussian can be given here", default = {'Population': 'Regular', 'Density': 'Current'}, type = dict)
	convert_chk = Option(help = "Whether to create an .fchk file at the end of the calculation", default = True, type = bool)
	keep_chk = Option(help = "Whether to keep the .chk file at the end of the calculation. If False, the .chk file will be automatically deleted, but not before it is converted to an .fchk file (if convert_chk is True)", default = False, type = bool)
	
	@property
	def charge(self):
		"""
		The molecule/system charge that we'll actually be using in the calculation.
		
		Unlike the charge attribute, this property will translate "auto" to the actual charge to be used.
		"""
		return self._charge if self._charge != "auto" else self.input_file.charge
	
	@property
	def multiplicity(self):
		"""
		The molecule/system multiplicity that we'll actually be using in the calculation.
		
		Unlike the multiplicity attribute, this property will translate "auto" to the actual multiplicity to be used.
		"""
		return self._multiplicity if self._multiplicity != "auto" else self.input_file.multiplicity
	
	@property
	def external_ECPs(self):
		"""
		The list of basis set Configurable objects we'll be using in the calculation for effective core potentials.
		
		This property will translate the names of the basis sets, under self._extended_ECPs, to the actual objects.
		"""
		try:
			return [self.available_basis_sets.get_config(extended_ECP) for extended_ECP in self._external_ECPs]
		except Exception:
			raise Configurable_exception(self, "could not load external ECP")
		
	@property
	def external_basis_sets(self):
		"""
		The list of basis set Configurable objects we'll be using in the calculation.
		
		This property will translate the names of the basis sets, under self._extended_basis_sets, to the actual objects.
		"""
		try:
			return [self.available_basis_sets.get_config(extended_basis_set) for extended_basis_set in self._external_basis_sets]
		except Exception:
			raise Configurable_exception(self, "could not load external basis set")
		
	@property
	def model_chemistry(self):
		"""
		The 'model chemistry' to be used by the calculation, this is a string containing both the functional and bases set (separated by /).
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
			route_parts.append("SCRF=(Solvent={})".format(self.solvent))
		
		# Finally, add any free-form options.
		for option in self.options:
			if self.options[option] == "":
				# Blank option, just add the keyword.
				route_parts.append(option)
			
			# Skip None options.
			elif self.options[option] is not None: 
				# Option with options, add both.
				route_parts.append("{}={}".format(option, self.options[option]))
				
		# Convert to string and return.
		return " ".join(route_parts)
		
	@classmethod
	def safe_name(self, file_name):
		"""
		Get a filename safe for Gaussian.
		
		What constitutes a safe name from Gaussian's point of view is not entirely clear, to play it safe we'll only allow alpha-numeric characters, dots and underscores.
		
		:param file_name: The file name to make safe, note that a path (containing path separators) will not maintain its structure after a call to safe_name().
		:return: The safe path name.
		"""
		# Adapted from https://stackoverflow.com/questions/7406102/create-sane-safe-filename-from-any-unsafe-string
		safe_chars = "._"
		return "".join([char if char.isalnum() or char in safe_chars else "_" for char in file_name])
	
	
	def _submit_init(self, output, input_str, name):
		"""
		Step 1/4 of the submission process, this method is called to set-up submission.
		
		Performs init for Gaussian calculations.
		
		:param output: Path to perform the calculation in.
		:param input_str: String containing a calculation file to submit.
		:param name: Name of the submitted file (should include extension if any). This is used eg as the base name for files created during the calculation.
		"""
		super()._submit_init(output, input_str, name)
		# Load our input file.
		try:
			self.input_file = Gaussian_input_parser(input_str)
		except Exception:
			raise Submission_error(self, "Failed to read input file")
		
	def _submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		"""
		super()._submit_pre()
		
		# Decide on our file names.
		self.chk_file_name = self.safe_name(self.name + ".chk")
		self.com_file_name = self.safe_name(self.name + ".com")
		
		# Get and load our com file template.
		self.com_file_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/gaussian_input.mako").render_unicode(calculation = self)
		
		
class Gaussian_input_parser():
	"""
	Class for parsing Gaussian input files.
	"""
	
	# Each 'section' in Gaussian is separated by a double newline; this seems to be the only consistent delimiter.
	SECTION_SEPARATOR = "\n\n"
	
	def __init__(self, file_str = None):
		"""
		Constructor for Gaussian input files.
		
		:param file_str: String containing a Gaussian input file.
		"""
		# Link 0 commands that appear at the top of the input file, we currently do no parsing on these (because we ignore them anyway)
		# Link 0 commands are stored as they appear (as a string, possibly containing newlines).
		self.link_0 = None
		
		# The route section, describes the calculation to be performed.
		## Most of the route is ignored as we set it ourselves as part of the submission process, but some options (geom=connectivity, genECP etc) control the format of the rest of the file, so we do parse these.
		## This is a dictionary of options; where 'key: value' translates to 'key=value' in the gaussian file. If value is None, then the translation is to 'key'.
		# The route section is currently not parsed (it is just a string).		
		self.route = None
		
		# The title of the input file.
		self.title = None
		
		# The multiplicity and charge of the molecule, appears at the top of the geometry section as 'charge, mult' or 'charge mult'.
		self.multiplicity = None
		self.charge = None
		
		# The geometry (atoms and charge etc) section as a string (almost certainly containing newlines).
		self.geometry = None
		
		# These are additional sections that can optionally appear in some input files (connectivity, basis set etc).
		self.additional_sections = []
		# If we've been given a file, load it now.
		if file_str is not None:
			self.load(file_str)
		
	@property
	def title(self):
		"""
		Get the title section of this input file.
		
		None and empty string values are translated to a single whitespace character (because otherwise they will be interpreted wrong).
		"""
		return self._title if self._title is not None and self._title is not "" else " "
	
	@title.setter
	def title(self, value):
		"""
		Set the title section of this input file.
		"""
		self._title = value
	
	def load(self, file_str):
		"""
		Load a Gaussian input file.
		
		:param file_str: String containing a Gaussian input file.
		
		"""
		# First, split on our delimeter.
		sections = file_str.split("\n\n")
		
		link_0_lines = []
		route_lines = []
		# Split the first section into link 0 and route (link 0 starts with %, the first line not to start without % is route).
		if len(sections) > 0:
			in_route = False
			for line in sections[0].split("\n"):
				# Check first char.
				if not in_route and line[:1] != "%":
					in_route = True
					
				# Add to one of our two lists.
				if not in_route:
					link_0_lines.append(line)
				else:
					route_lines.append(line)
		
		# Now set.
		self.link_0 = "\n".join(link_0_lines) if len(link_0_lines) > 0 else None
		self.route = "\n".join(route_lines) if len(route_lines) > 0 else None
		
		# Set title.
		self.title = sections[1] if len(sections) > 1 else None
		
		# Geometry (the first line contains charge and mult).
		geometry_section = (sections[2] if len(sections) > 2 else "").split("\n", 1)
		
		# We'll first try to split on comma (,) for charge,mult.
		charge_mult = geometry_section[0].split(",", 1)
		# If we didn't get what we want, try on whitespace.
		if len(charge_mult) != 2:
			charge_mult = geometry_section[0].split(" ", 1)
		
		# Now try and set.
		try:
			self.charge = int(charge_mult[0])
			self.multiplicity = int(charge_mult[1])
		except Exception:
			raise Silico_exception("Unable to determine charge and multiplicity from '{}'".format(charge_mult))
		
		# The rest of the geometry section contains atoms.
		if len(geometry_section) != 2 or geometry_section[1] == "":
			raise Silico_exception("Gaussian input file does not appear to contain any atoms")
		
		self.geometry = geometry_section[1]
		
		# And anything else.
		self.additional_sections = sections[3:]
		