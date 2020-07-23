from silico.submit import Configurable_target
from silico.config.configurable.option import Option


class Extended_basis_set(Configurable_target):
	"""
	Top-level class for basis set targets.
	"""
	
	CLASS_HANDLE = ("basis_set",)
	
	basis_set = Option(help = "Basis set data", required = True, type = lambda data: str(data).strip())
	