import math

class Spectroscopy_graph():
	"""
	Top level class for graphing spectroscopy results (energy vs intensity).
	
	For generating pictures of these graphs, see silico.image.spectroscopy
	"""
	
	def __init__(self, coordinates):
		"""
		Constructor for Spectroscopy_graph objects
		
		:param coordinates: A list of (energy, intensity) tuples to plot. The units of energy and intensity are irrelevant here (but should be consistent.
		"""
		self.coordinates = coordinates
		
	@classmethod
	def gaussian(self, a, b, c, x):
		"""
		An implementation of the gaussian function.
		
		:param a: The maximum height of the peak.
		:param b: The x position of the center of the peak.
		:param c: The width of the peak.
		:param x: The x-value to plot for.
		:return: The corresponding y value.
		"""
		return a * math.exp(-( (x - b)**2 / (2 * c**2) ))
	
	@classmethod
	def gaussian_x(self, a, b, c, y):
		"""
		An implementation of the gaussian function rearranged to find x
		
		:param a: The maximum height of the peak.
		:param b: The x position of the center of the peak.
		:param c: The width of the peak.
		:param y: The y-value to plot for.
		:return: A tuple of the two corresponding x values.
		"""
		# We have two solutions of the form x = ± d + b
		# Calculate d first.
		d = math.sqrt( -math.log( y/a ) * 2 * c **2  )
		
		# Now return our two solutions.
		return (-d + b, d + b)
		
		
	@classmethod
	def fwhm_to_c(self, fwhm):
		"""
		Convert a full width at half maximum to the corresponding c-value in the Gaussian function.
		
		:param fwhm: The desired full width at half maximum (in nm).
		:return: The c-value.
		"""
		return fwhm / (2 * math.sqrt(2 * math.log(2)))