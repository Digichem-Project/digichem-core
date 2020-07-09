from silico.submit import Configurable_target


class Extended_basis_set(Configurable_target):
	"""
	Top-level class for basis set targets.
	"""
	
	def _post_init(self, *, basis_set, **kwargs):
		super()._post_init(**kwargs)
		self.basis_set = basis_set.strip()