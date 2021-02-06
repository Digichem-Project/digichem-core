from silico.file import File_maker
from PIL import Image

class Image_maker(File_maker):
	"""
	Top level class for image maker objects.
	"""
	
	# Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
	output_file_type = "image"
	
	def __init__(self, *args, **kwargs):
		"""
		General constructor for Image_maker objects.
		
		These object are used to represent and make images. Note that a single Image_maker object can represent several image_files
		
		:param output: A path to an output file to write to. How exactly this operates depends on the inheriting class, it is often only used as a base file name. See the class you are using.
		:param dont_modify: Flag that modifies how image creation works. If True, no new images will be written to file.
		:param use_existing: Flag that modifies how image creation works. If True, existing files will be preferentially used if available (set to False to force overwriting existing files).
		"""
		super().__init__(*args, existing_file = None, **kwargs)
		
	def get_image(self, name = 'file'):
		"""
		Get the path to one of the images that this class represents, rendering the image to file first if necessary.
		
		The functioning of this method is controlled by the dont_modify & use_existing flags.
		
		You can also use the normal python attribute mechanism (either through getattr() or dot notation) to get these paths.
		
		:raises KeyError: If name is not the name of one of the images this class represents.
		:param name: The name of an image to get. Depends on the images in self.file_path.
		:return: A pathlib Path object pointing to the image represented by name, or None if no image could be created.
		"""
		return self.safe_get_file(name)
	
	@property
	def creation_message(self):
		"""
		A short message that may (depending on log-level) be printed to the user before make_files() is called.
		"""
		return "Rendering {} to file(s)".format(self.output)
	
	
		
	def get_constrained_size(self, max_width, max_height):
		"""
		Get the maximum possible dimensions of this image that retain the same aspect ratio that don't exceed a set of dimensions.
		
		:param max_width: The maximum width to resize to; set to math.inf for no max.
		:param max_heigh: The maximum height to resize to; set to math.inf for no max.
		"""
		# Weasyprint hack for broken max-width.
				
		# Now we open our diagram with pillow.
		im = Image.open(self.get_image())
		
		width, height = im.size
		
		scale_factor = min(max_width/width, max_height/height)
		
		return (width *scale_factor, height * scale_factor)
	
		