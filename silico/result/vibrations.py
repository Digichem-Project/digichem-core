from silico.result import Result_container
from silico.result.base import Result_object
from silico.image.spectroscopy import Frequency_graph_maker
from pathlib import Path
from itertools import zip_longest

class Vibration_list(Result_container):
	"""
	Class for representing a group of molecular vibrations.
	"""
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
	@property
	def negative_frequencies(self):
		"""
		Get a Vibration_list object of the vibrations in this list that have negative frequencies (they are imaginary).
		"""
		return type(self)([vibration for vibration in self if vibration.frequency < 0])
	
	def set_file_options(self, output_dir, output_name, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# First our states diagram.
		self._files['simulated_IR_graph'] = Frequency_graph_maker.from_image_options(
			Path(output_dir, output_name + ".simulated_frequencies.png"),
			vibrations = self,
			**kwargs
		)
		
	@property
	def simulated_IR_graph(self):
		return self.get_file('simulated_IR_graph')

	
	@classmethod
	def from_parser(self, parser):
		"""
		Get an Vibration_list object from an output file parser.
		
		:param parser: An output file parser.
		:return: An Vibration_list object. The list will be empty if no vibration frequency data is available.
		"""
		return self(Vibration.list_from_parser(parser))
	
class Vibration(Result_object):
	"""
	Class for representing vibrational frequencies.
	"""
	
	def __init__(self, level, frequency, intensity, symmetry = None):
		"""
		Vibration object constructor.
		
		:param level: The level of the vibration.
		:param frequency: Frequency of the vibration (cm-1).
		:param intensity: Intensity of the vibration in IR (km mol-1)
		:param symmetry: Symmetry term of the vibration (string).
		"""
		self.level = level
		self.symmetry = symmetry
		self.frequency = frequency
		self.intensity = intensity
		
		
	@classmethod
	def list_from_parser(self, parser):
		"""
		Create a list of Vibration objects from an output file parser.
		
		:param parser: An output file parser.
		:return: A list of Vibration objects. The list will be empty if no frequency data is available.
		"""
		try:
			return [self(index+1, symmetry, frequency, intensity) for index, (symmetry, frequency, intensity) in enumerate(zip_longest(getattr(parser.data, 'vibsyms', []), parser.data.vibfreqs, parser.data.vibirs, fillvalue = None))]
		except AttributeError:
			return []