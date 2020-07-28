from silico.submit.basis import Extended_basis_set
from silico.config.configurable.option import Option


class Gaussian_basis_set(Extended_basis_set):
	"""
	External basis sets (including possibly ECP) for Gaussian.
	"""
	
	CLASS_HANDLE = ("gaussian_basis",)
	
	_ECP = Option("ECP", help = "ECP (effective core potential) data", default = None, type = str)
	
	@property
	def ECP(self):
		return self._ECP.strip() if self._ECP is not None else None