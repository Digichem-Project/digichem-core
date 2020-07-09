from mako.lookup import TemplateLookup

import silico
from silico.exception import Configurable_exception
from silico.exception.base import Submission_error, Silico_exception
from silico.submit.calculation import Calculation_target

class Gaussian_DFT(Calculation_target):
	"""
	DFT (density functional theory) calculations with Gaussian.
	"""
	# Identifying handle.
	CLASS_HANDLE = ("Gaussian-DFT",)
	
	# A list of strings describing the expected input file types (file extensions) for calculation's of this class. The first item of this list will be passed to obabel via the -o flag. 
	INPUT_FILE_TYPES = ["gau", "com", "gjf", "gjc"]
		
	def _post_init(self,
		*args,
		use_chk = None,
		CPU_list = None,
		num_CPUs = None,
		calculation_keywords = None,
		functional = None,
		basis_set = None,
		extended_basis_sets = None,
		extended_ECPs = None,
		convert_chk = None,
		keep_chk = None,
		multiplicity = None,
		charge = None,
		options = None,
		**kwargs
	):
		"""
		Constructor for Gaussian DFT calculations.	
		"""
		super()._post_init(*args, **kwargs)
		
		self.use_chk = use_chk if use_chk is not None else True
		
		self.CPU_list = CPU_list if CPU_list is not None else []
		# If both CPU_list and num_CPUs are empty, set to auto.
		self.num_CPUs = "auto" if num_CPUs is None and len(self.CPU_list) == 0 else num_CPUs
				
		self.calculation_keywords = calculation_keywords if calculation_keywords is not None else []
		self.functional = functional
		self.basis_set = basis_set
		self.options = options if options is not None else {}
		self.keep_chk = keep_chk if keep_chk is not None else True
		self.convert_chk = convert_chk if convert_chk is not None else True
		
		# Charge and multiplicity, auto means to take the value from the input file.
		self.multiplicity = multiplicity if multiplicity is not None else "auto"
		self.charge = charge if charge is not None else "auto"
		
		# Save our basis set.
		self.extended_basis_sets = [self.available_basis_sets.get_config(extended_basis_set) for extended_basis_set in extended_basis_sets] if extended_basis_sets is not None else []
		self.extended_ECPs = [self.available_basis_sets.get_config(extended_ECP) for extended_ECP in extended_ECPs] if extended_ECPs is not None else []
	
	@property
	def real_charge(self):
		"""
		The molecule/system charge that we'll actually be using in the calculation.
		
		Unlike the charge attribute, this property will translate "auto" to the actual charge to be used.
		"""
		return self.charge if self.charge != "auto" else self.input_file.charge
	
	@property
	def real_multiplicity(self):
		"""
		The molecule/system multiplicity that we'll actually be using in the calculation.
		
		Unlike the multiplicity attribute, this property will translate "auto" to the actual multiplicity to be used.
		"""
		return self.multiplicity if self.multiplicity != "auto" else self.input_file.multiplicity
		
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
		route_parts = self.calculation_keywords
		
		# Model chemistry
		route_parts.append(self.model_chemistry)
		
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
		
	@property
	def CPU_list(self):
		"""
		A list of integer indices identifying specific CPUs to use for the calculation.
		"""
		return self._CPU_list
	
	@CPU_list.setter
	def CPU_list(self, value):
		"""
		Set the list of integers indices identifying specific CPUs to use for the calculation.
		"""
		if len(value) != 0 and self._num_CPUs is not None:
			raise Configurable_exception(self, "'CPU_list' and 'num_CPUs' cannot be specified simultaneously")

		# Set.
		self._CPU_list = value
	
	@property
	def num_CPUs(self):
		"""
		The number of CPUs to use for the calculation.
		"""
		if len(self._CPU_list) == 0:
			return super().num_CPUs
		else:
			return len(self._CPU_list)
	
	@num_CPUs.setter
	def num_CPUs(self, value):
		"""
		Set the number of CPUs to use for the calculation. In addition to an exact integer amount, the string "auto" can also be supplied, in which case all available CPUs will be used.
		"""
		if value is not None and len(self.CPU_list) != 0:
			raise Configurable_exception(self, "'CPU_list' and 'num_CPUs' cannot be specified simultaneously")

		# Set.
		super(Gaussian_DFT, self.__class__).num_CPUs.fset(self, value)
	
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
		self.chk_file_name = self.safe_name(self.name + ".chk") if self.use_chk else None
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
		