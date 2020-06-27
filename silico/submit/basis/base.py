from silico.submit import Configurable_target


class Extended_basis_set(Configurable_target):
	"""
	Top-level class for basis set targets.
	"""
	
	def __init__(self, *, basis_set, **kwargs):
		super().__init__(**kwargs)
		self.basis_set = basis_set.strip()