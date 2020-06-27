from silico.submit.calculation.base import Calculation_target

class Calculation_series(Calculation_target):
	"""
	A pseudo calculation that automatically runs a number of calculations in series.
	"""
	
	CLASS_HANDLE = ["Series"]
	
	def __init__(self, 
		*args,
		calculations,
		configs,
		silico_options,
		available_basis_sets,
		**kwargs
	):
		"""
		Constructor for Calculation_series objects.
		
		These Configurable_targets represent a number of calculations to run in series.
		
		:param calculations: A list of config IDs which are the calculations this series represents.
		"""
		super().__init__(*args, configs = configs, silico_options = silico_options, available_basis_sets = available_basis_sets, **kwargs)
		self.calculations = [Calculation_target.from_name_in_configs(calculation_name, configs, silico_options = silico_options, available_basis_sets = available_basis_sets) for calculation_name in calculations]
		
		
	def prepare(self):
		"""
		Prepare this calculation target for submission.
		
		For Calculation_series, we return the calcs we represent.
		
		:return: A list of ready-to-go calculation targets.
		"""
		return self.calculations