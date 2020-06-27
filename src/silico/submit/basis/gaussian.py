from silico.submit.basis import Extended_basis_set


class Gaussian_basis_set(Extended_basis_set):
	"""
	External basis sets (including possibly ECP) for Gaussian.
	"""
	
	CLASS_HANDLE = "Gaussian"
	
	def __init__(self, *, ECP = None, **kwargs):
		super().__init__(**kwargs)
		self.ECP = ECP.strip() if ECP is not None else None