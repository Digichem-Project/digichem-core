from silico.submit.calculation.base import Concrete_calculation
from silico.config.configurable.option import Option
from silico.submit.base import Memory
from mako.lookup import TemplateLookup
import silico


				
class Turbomole_memory(Memory):
	"""
	A class for representing memory quantities in formats suitable for turbomole.
	"""
	
	# The units that we know about.
	UNITS = {
		' TB': 1000000000000,
		' GB': 1000000000,
		' MB': 1000000,
		' KB': 1000,
		' B': 1
		}

class Turbomole_CC(Concrete_calculation):
	"""
	Coupled Cluster calculations with Turbomole.
	"""
	
	# Identifying handle.
	CLASS_HANDLE = ("turbomole",)
	
	# A list of strings describing the expected input file types (file extensions) for calculation's of this class. The first item of this list will be passed to obabel via the -o flag. 
	INPUT_FILE_TYPES = ["tmol"]
	
	# Configurable options.
	memory = Option(help = "The amount of memory to use for the calculation", required = True, type = Turbomole_memory, rawtype = str)
	basis_set = Option(help = "The basis set to use", type = str)
	charge = Option(help = "Set the molecule charge", default = 0, type = int)
	redundant_internal_coordinates = Option(help = "Whether to use redundant internal coordinates", type = bool, default = True)
	methods = Option(help = "Method keywords and options from the define general menu, including scf, mp2, cc etc.", type = dict, default = {})
	
	def _submit_pre(self):
		"""
		Step 2/4 of the submission process, this method is called before submission begins.
		"""
		super()._submit_pre()		
		# Get and load our com file template.
		self.com_file_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/turbomole/define_driver.mako").render_unicode(calculation = self)