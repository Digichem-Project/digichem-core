from silico.submit.basis import Extended_basis_set


class Gaussian_basis_set(Extended_basis_set):
	"""
	External basis sets (including possibly ECP) for Gaussian.
	"""
	
	CLASS_HANDLE = ("gaussian_basis",)
	
	def _post_init(self, *, ECP = None, **kwargs):
		super()._post_init(**kwargs)
		self.ECP = ECP.strip() if ECP is not None else None