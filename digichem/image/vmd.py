
# General imorts.
from pathlib import Path
import subprocess
import math
from logging import getLogger
from PIL import Image
import numpy as np
import logging
import pkg_resources
from uuid import uuid4
import os
from math import fabs

# Silico imports.
from silico.exception.base import File_maker_exception
from silico.file import File_converter
import silico

class VMD_image_maker(File_converter):
	"""
	Class for generating image files from Gaussian outputs using VMD.
	
	Nearly all of the work here is done by other programs, primarily VMD (https://www.ks.uiuc.edu/Research/vmd/). These classes are mostly just wrappers.
	"""
	
	# The name of the vmd/tcl script (relative to the vmd script folder) used to render images. Inheriting classes should set this to an appropriate script file.
	vmd_script = ""
	# The filename extension to use for molecule scene files produced by vmd.
	scene_file_extension = ".scene"
	# 'Path' to the vmd executable.
	vmd_execuable = "vmd"
	# 'Path' to the tachyon executable.
	tachyon_executable = "tachyon"
	# The initial resolution at which an image is rendered. This test image will then be discarded once relevant info has been extracted.
	test_resolution = 300
	
	# Text description of our input file type, used for error messages etc. This can be changed by inheriting classes.
	input_file_type = "cube"
	# Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
	output_file_type = "image"
	
	# Name of the section where we get some specific configs.
	options_name = "orbital"
	
	def __init__(self, *args, cube_file = None, rotations = None, auto_crop = True, rendering_style = "pastel", resolution = 1024, also_make_png = True, isovalue = 0.2, **kwargs):
		"""
		Constructor for Image_maker objects.
		
		:param output: The path to write image files to. As more than one image can be created by this class, this name is used as the base name and is modified for each image. The suffix will be honoured (rendered imges will use it as a hint for their output format.
		:param cube_file: The path to a cube_file to use to render images.
		:param rotations: A list of tuples of rotations, where the first index in the tuple specifies the axis to rotate about and the second is the angle to rotate (in radians).
		:param auto_crop: If False, images will not have excess white space cropped.
		:param rendering_style: A string describing the rendering style to use, either 'silico' or 'gaussian'.
		:param resolution: The max width or height of the rendered images in pixels.
		:param also_make_png: If True, additional images will be rendered in PNG format. This option is useful to generate higher quality images alongside more portable formats. If 'output' is a .png file, then it is wise to set this option to False (otherwise two png files will be rendered, which is a waste).
		:param isovalue: The isovalue to use for rendering isosurfaces. Has no effect when rendering only atoms.
		"""
		super().__init__(*args, input_file = cube_file, **kwargs)
		# Save our translations list.
		#self.translations = translations if translations is not None else (0,0,0)
		self.translations = (0,0,0)
		# Save our rotations list.
		# We save the negative of the rotations because VMD rotates the opposite way to us.
		self.rotations = [(axis, math.degrees(-theta)) for axis, theta in rotations] if rotations is not None else []
		
		# Some options that control how we function.
		self.auto_crop = auto_crop
		self.rendering_style = rendering_style
		self.target_resolution = resolution
		self.also_make_png = also_make_png
		self.isovalue = isovalue		
				
		# The colours we will be using for rendering (depends on rendering_stype and is only for reference; the actual definition of each style is in VMD).
		if self.rendering_style == "gaussian":
			self.primary_colour = "green"
			self.secondary_colour = "red"
		elif self.rendering_style == "vesta":
			self.primary_colour = "blue"
			self.secondary_colour = "yellow"
		else:
			self.primary_colour = "blue"
			self.secondary_colour = "red"
				
		# These 4 attributes are file paths to the four images we create.
		# We'll keep the same file extension as the was given to us.
		self.file_path = {
			'x0y0z0': self.output.with_suffix(".x0y0z0" + self.output.suffix),
			'x90y0z0': self.output.with_suffix(".x90y0z0" + self.output.suffix),
			'x0y90z0': self.output.with_suffix(".x0y90z0" + self.output.suffix),
			'x45y45z45': self.output.with_suffix(".x45y45z45"  + self.output.suffix),
			# There are higher quality PNG versions.
 			'x0y0z0_big': self.output.with_suffix(".x0y0z0.png"),
 			'x90y0z0_big': self.output.with_suffix(".x90y0z0.png"),
 			'x0y90z0_big': self.output.with_suffix(".x0y90z0.png"),
 			'x45y45z45_big': self.output.with_suffix(".x45y45z45.png")
			}
	
	@classmethod
	def from_options(self, output, *, cube_file = None, rotations = None, options, **kwargs):
		"""
		Constructor that takes a dictionary of config like options.
		"""		
		return self(
			output,
			cube_file = cube_file,
			rotations = rotations,
			auto_crop = options['molecule_image']['auto_crop'],
			rendering_style = options['molecule_image']['rendering_style'],
			resolution = options['molecule_image']['resolution'],
			isovalue = options['molecule_image'][self.options_name]['isovalue'],
			use_existing = options['molecule_image']['use_existing'],
			dont_modify = options['image']['dont_modify'],
			**kwargs
		)
		
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
	def vmd_script_path(self):
		"""
		Get the file path to the VMD script to be used by this class to render images (as a pathlib Path).
		"""
		return Path(pkg_resources.resource_filename('silico', 'data/vmd/{}'.format(self.vmd_script)))
	
	@property
	def tcl_common_path(self):
		"""
		Get the file path to the Tcl script that is common to all of our other Tcl scripts (as a pathlib Path).
		"""
		return Path(pkg_resources.resource_filename('silico', 'data/vmd/common.tcl'))
	
	@property
	def creation_message(self):
		"""
		A short message that may (depending on log-level) be printed to the user before make_files() is called.
		"""
		return "Rendering {} to file(s)".format(self.output)
		
	def make_files(self):
		"""
		Make the image files referenced by this object.
		
		The new image will be written to file.
		"""
		# Run VMD, which writes a plain text description of our scene to file.
		self.run_VMD_script()
		
		# Next run tachyon for each of our 4 scene files, which creates tga files.
		for image_name in ['x0y0z0', 'x90y0z0', 'x0y90z0', 'x45y45z45']:
			image_path = self.file_path[image_name]
			try:
				# First we'll render a test image at a lower resolution. We'll then crop it, and use the decrease in final resolution to know how much bigger we need to render in our final image to hit our target resolution.
				# Unless of course auto_crop is False, in which case we use our target resolution immediately.
				resolution = self.test_resolution if self.auto_crop else self.target_resolution
				self.run_tachyon_renderer(image_path.with_suffix(self.scene_file_extension), image_path.with_suffix(".tga"), resolution)
				
				if self.auto_crop:
					# Load the test image and autocrop it.
					with Image.open(image_path.with_suffix(".tga"), "r") as test_im:
						small_test_im = self.auto_crop_image(test_im)
						
					# Get the cropped size. We're interested in the largest dimension, as this is what we'll output as.
					cropped_resolution = max(small_test_im.size)
					
					# From this we can work out the ratio between our true resolution and the resolution we've been asked for.
					resolution_ratio = cropped_resolution / self.test_resolution
					
					# Now we can re-render asking for a larger resolution so that when we crop again, we should get our target resolution.
					# Note that resolution ratio is a fraction, so self.target_resolution / resolution_ratio results in an increase.
					self.run_tachyon_renderer(image_path.with_suffix(self.scene_file_extension), image_path.with_suffix(".tga"), self.target_resolution / resolution_ratio)
				
				# Delete the now unneeded scene file.
				os.remove(image_path.with_suffix(self.scene_file_extension))
			except Exception:
				raise File_maker_exception(self, "Error in tachyon rendering")
			
			# Convert to a better set of formats.
			# Open the file we just rendered.
			with Image.open(image_path.with_suffix(".tga"), "r") as im:
				
				# If we've been asked to autocrop, do so.
				if self.auto_crop:
					try:
						cropped_image = self.auto_crop_image(im)
					except Exception:
						raise File_maker_exception(self, "Error in post-rendering auto-crop")
				else:
					cropped_image = im
				
				# Now save in out main output format.
				cropped_image.save(image_path)
				
				# And also as a higher quality png if we've been asked to.
				if self.also_make_png:
					cropped_image.save(self.file_path[image_name + "_big"])
				
				# And delete the old .tga (unless we've just rendered to .tga!).
				if image_path.suffix != ".tga":
					os.remove(image_path.with_suffix(".tga"))
		
	@property
	def prepared_translations(self):
		"""
		Our list of translations in a from ready for VMD/Tcl.
		"""
		return "{},{},{}".format(self.translations[0], self.translations[1], self.translations[2])
	
	@property
	def prepared_rotations(self):
		"""
		Our list of rotations in a form ready for VMD/Tcl.
		"""
		# This is a stringified form of our list of tuples that we'll pass to VMD (or rather the TCL interpreter).
		rot_string = " ".join(["{},{}".format(axis, angle) for axis, angle in self.rotations])
		return rot_string
		
	
	@property
	def VMD_signature(self):
		"""
		The signature pass to subprocess.run used to call VMD. Inheriting classes should write their own implementation.
		"""
		return ""
	
	
	def run_VMD_script(self):
		"""
		Called as part of the make() method, inheriting classes can implement this method if they have a different VMD call signature.
		
		Also see VMD_signature, which may do what you need.
		
		:return: The CompletedProcess object. 
		"""
		# Setting these variables greatly improves VMD startup time.
		os.environ["VMDNOOPTIX"] = "1"
		os.environ["VMDNOOSPRAY"] = "1"
		
		# Run VMD, which renders our image for us.
		try:
			return subprocess.run(
				self.VMD_signature,
				# We normally capture and discard stdout (because VMD is VERY verbose), but if we're at a suitable log level, we'll print it.
				# Nothing useful appears to be printed to stderr, so we'll treat it the same as stdout.
				stdout = subprocess.DEVNULL if getLogger(silico.logger_name).level > logging.DEBUG else None,
				stderr = subprocess.STDOUT,
				universal_newlines = True,
				# VMD has a tendency to sigsegv when closing with VMDNOOPTIX set to on (even tho everything is fine) so we can't check retval sadly.
				#check = True
			)
		except FileNotFoundError:
			raise File_maker_exception(self, "Could not locate vmd executable '{}'".format(self.vmd_execuable))
		
	def run_tachyon_renderer(self, scene_file, tga_file, resolution):
		"""
		Called as part of the make() method.
		
		Take a molecule 'scene' file produced by VMD and render to an image with tachyon.
		
		:param scene_file: A VMD scene file to read in as input.
		:param tga_file: The filename to write to.
		:param resolution: The resolution to render with.
		"""
		# Sadly, tachyon is old and decrepit and doesn't parse the output filename correctly (but weirdly the input file is fine).
		# To get around this, we'll use a semi-random output name that we know is safe, and change-directory to the output directory immediately before calling tachyon.
		# This way we can ensure there will be no nasty file-names for tachyon to get stuck on."
		working_directory = tga_file.parent
		tmpfile_name = ".tachyon_output_" + uuid4().hex + ".tga"
		tmpfile_full_path = Path(tga_file.parent, tmpfile_name)
		try:
			try:
				# Now we can run tachyon.
				subprocess.run(
					[
						"{}".format(self.tachyon_executable),
						scene_file.relative_to(working_directory),
						"-aasamples", "12",
						"-res", "{}".format(resolution), "{}".format(resolution),
						"-o", tmpfile_name
					],
					stdout = subprocess.DEVNULL if getLogger(silico.logger_name).level > logging.DEBUG else None,
					stderr = subprocess.STDOUT,
					universal_newlines = True,
					check = True,
					cwd = working_directory
				)
			except FileNotFoundError:
				# THINK this can only occur if the command itself cannot be found, but I might be wrong.
				raise File_maker_exception(self, "Could not locate tachyon executable '{}'".format(self.tachyon_executable))
			
			# Render complete (sadly we don't appear to have any ability to check for errors).
			# Now rename our tmpfile.
			os.rename(tmpfile_full_path, tga_file)
		except Exception:
			# Something went wrong, remove our tmpfile.
			try:
				os.remove(tmpfile_full_path)
			except Exception:
				pass
			raise
				
	
	def auto_crop_image(self, im):
		"""
		Automatically remove excess white pixels around the outside of our image.
		
		:param im: A PIL Image object.
		:returns: A new PIL Image object with white-space removed.
		"""
		# (try to) load our image from file.
		#im = Image.open(file_path, "r")
		im.load()
				
		# This solution was inspired by https://stackoverflow.com/questions/14211340/automatically-cropping-an-image-with-python-pil
		image_data = np.asarray(im)
		# Axis 2 is the tuple/array of RGB values for each pixel. We want to know which ones are fully white (255,255,255), so we'll take the min of the 3 values and see if this is equal to 255. If it is, then all the pixels are 255		
		image_data_bw = image_data.min(axis = 2)
		
		# Now we just check which rows and columns are not pure white.
		non_empty_columns = np.where(image_data_bw.min(axis = 0) != 255)[0]
		non_empty_rows = np.where(image_data_bw.min(axis = 1) != 255)[0]
		
		# Get the bounding box of non-white stuff.
		cropBox = (min(non_empty_rows), max(non_empty_rows), min(non_empty_columns), max(non_empty_columns))
		
		# Copy the image data we want.
		image_data_new = image_data[cropBox[0]:cropBox[1]+1, cropBox[2]:cropBox[3]+1 , :]
		
		# And save over our old file.
		new_image = Image.fromarray(image_data_new)
		return new_image
		#new_image.save(file_path)
		

class Structure_image_maker(VMD_image_maker):
	"""
	Class for creating structure images.
	"""
	
	vmd_script ="generate_structure_images.tcl"
	
	@property
	def VMD_signature(self):
		"""
		The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
		"""
		return [
				"{}".format(self.vmd_execuable),
				"-dispdev", "none",
				"-e", "{}".format(self.vmd_script_path),
				"-args",
				"{}".format(self.input_file),
				"{}".format(self.tcl_common_path),
				"{}".format(self.rendering_style),
				"{}".format(self.prepared_translations),
				"{}".format(self.prepared_rotations),
				"{}".format(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension))
			]
		
		
class Orbital_image_maker(Structure_image_maker):
	"""
	Class for creating orbital images.
	"""
	
	vmd_script ="generate_orbital_images.tcl"
	
	@property
	def VMD_signature(self):
		"""
		The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
		"""
		return [
				"{}".format(self.vmd_execuable),
				"-dispdev", "none",
				"-e", "{}".format(self.vmd_script_path),
				"-args",
				"{}".format(self.input_file),
				"{}".format(self.tcl_common_path),
				"{}".format(self.rendering_style),
				"{}".format(self.isovalue),
				"{}".format(self.prepared_translations),
				"{}".format(self.prepared_rotations),
				"{}".format(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension))
			]
	
	
class Spin_density_image_maker(Orbital_image_maker):
	"""
	Class for creating spin density images.
	"""
		
	vmd_script ="generate_spin_images.tcl"
	
	# Name of the section where we get some specific configs.
	options_name = "spin"
	
	def __init__(self, *args, spin = "both", **kwargs):
		"""
		Constructor for Spin_density_image_maker objects.
		
		See Orbital_image_maker for a complete constructor.
		:param spin: A string indicating which net-spins to render, either 'positive', 'negative' or 'both'.
		"""
		super().__init__(*args, **kwargs)
		self.spin = spin
				
	@property
	def VMD_signature(self):
		"""
		The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
		"""
		return [
				"{}".format(self.vmd_execuable),
				"-dispdev", "none",
				"-e", "{}".format(self.vmd_script_path),
				"-args",
				"{}".format(self.input_file),
				"{}".format(self.tcl_common_path),
				"{}".format(self.rendering_style),
				"{}".format(fabs(self.isovalue)),
				"{}".format(self.spin),
				"{}".format(self.prepared_translations),
				"{}".format(self.prepared_rotations),
				"{}".format(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension)),
				"{}".format(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension))
			]
	
class Alpha_orbital_image_maker(Orbital_image_maker):
	pass
	#cubegen_type = "AMO"
	
	
class Beta_orbital_image_maker(Orbital_image_maker):
	pass
	#cubegen_type = "BMO"
	
	
class Combined_orbital_image_maker(VMD_image_maker):
	"""
	Class for creating images with both the HOMO and LUMO shown together.
	"""
	
	vmd_script ="generate_combined_orbital_images.tcl"
	
	def __init__(self, *args, HOMO_cube_file = None, LUMO_cube_file = None, **kwargs):
		"""
		Constructor for combined orbital image maker objects.
		
		:param output: Path to write to. See the constructor for VMD_image_maker for how this works.
		:param HOMO_cube_file: Path to the HOMO cube file.
		:param LUMO_cube_file: Path to the LUMO cube file.
		:param *args: See the constructor for VMD_image_maker for further options.
		:param **kwargs: See the constructor for VMD_image_maker for further options.
		"""
		super().__init__(*args, cube_file = None, **kwargs)
		self.HOMO_cube_file = HOMO_cube_file
		self.LUMO_cube_file = LUMO_cube_file
		
	@classmethod
	def from_options(self, output, *, HOMO_cube_file = None, LUMO_cube_file = None, rotations = None, options, **kwargs):
		"""
		Constructor that takes a dictionary of config like options.
		"""		
		return self(
			output,
			HOMO_cube_file = HOMO_cube_file,
			LUMO_cube_file = LUMO_cube_file,
			rotations = rotations,
			auto_crop = options['molecule_image']['auto_crop'],
			rendering_style = options['molecule_image']['rendering_style'],
			resolution = options['molecule_image']['resolution'],
			isovalue = options['molecule_image'][self.options_name]['isovalue'],
			use_existing = options['molecule_image']['use_existing'],
			dont_modify = options['image']['dont_modify'],
			**kwargs
		)
	
	def check_can_make(self):
		"""
		Check whether it is feasible to try and render the image(s) that we represent.
		
		Reasons for rendering not being possible are varied and are up to the inheriting class, but include eg, a required input (cube, fchk) file not being given.
		
		This method returns nothing, but will raise a File_maker_exception exception if the rendering is not possible.
		"""
		if self.HOMO_cube_file is None or self.HOMO_cube_file.safe_get_file() is None:
			raise File_maker_exception(self, "No HOMO cube file is available")
		if self.LUMO_cube_file is None  or self.LUMO_cube_file.safe_get_file() is None:
			raise File_maker_exception(self, "No LUMO cube file is available")
	
	@property
	def VMD_signature(self):
		"""
		The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
		"""
		return [
			"{}".format(self.vmd_execuable),
			"-dispdev", "none",
			"-e", "{}".format(self.vmd_script_path),
			"-args",
			"{}".format(self.HOMO_cube_file),
			"{}".format(self.LUMO_cube_file),
			"{}".format(self.tcl_common_path),
			"{}".format(self.rendering_style),
			"{}".format(fabs(self.isovalue)),
			"{}".format(self.prepared_translations),
			"{}".format(self.prepared_rotations),
			"{}".format(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension)),
			"{}".format(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension)),
			"{}".format(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension)),
			"{}".format(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension))
		]

	
class Dipole_image_maker(Structure_image_maker):
	"""
	Class for creating dipole images.
	"""
	
	vmd_script ="generate_dipole_images.tcl"
	
	def __init__(self, *args, dipole_moment = None, **kwargs):
		"""
		Constructor for Dipole_image_maker objects.
		
		:param output: Path to write to. See the constructor for VMD_image_maker for how this works.
		:param cube_file: A Gaussian cube file to use to render the new images.
		:param dipole_moment: A Dipole_moment object that will be rendered as a red arrow in the scene.
		:param *args: See the constructor for VMD_image_maker for further options.
		:param **kwargs: See the constructor for VMD_image_maker for further options.
		"""
		super().__init__(*args, **kwargs)
		self. dipole_moment = dipole_moment
		
	@classmethod
	def from_image_options(self, output, *, cube_file, dipole_moment = None, existing_file = None, rotations = None, options = None, **kwargs):
		"""
		An alternative constructor that discards any additional keyword arguments.
		"""
		return self(
			output, 
			cube_file = cube_file,
			dipole_moment = dipole_moment,
			existing_file = existing_file,
			rotations = rotations,
			**self._get_config(options['molecule_image']),
			**options['image']
		)
		
	def check_can_make(self):
		"""
		Check whether it is feasible to try and render the image(s) that we represent.
		
		Reasons for rendering not being possible are varied and are up to the inheriting class, but include eg, a required input (cube, fchk) file not being given.
		
		This method returns nothing, but will raise a File_maker_exception exception if the rendering is not possible.
		""" 
		super().check_can_make()
		
		# Also make sure we have a dipole.
		if self.dipole_moment is None:
			raise File_maker_exception(self, "No dipole moment is available.")
		
		
	@property
	def VMD_signature(self):
		"""
		The arguments which we'll pass to VMD, inheriting classes can implement this method if they have a different VMD call signature.
		"""
		return [
			"{}".format(self.vmd_execuable),
			"-dispdev", "none",
			"-e", "{}".format(self.vmd_script_path),
			"-args",
			"{}".format(self.input_file),
			"{}".format(self.tcl_common_path),
			"{}".format(self.rendering_style),
			"{}".format(self.prepared_translations),
			"{}".format(self.prepared_rotations),
			# We don't use the normal origin_coords/vector_coords because these are already rotated, while we want/need to do this rotation with the camera in VMD.
			# Hence use _origin_coords/_vector_coods, which aren't rotated.
			"{} {} {}".format(self.dipole_moment._origin_coords[0], self.dipole_moment._origin_coords[1], self.dipole_moment._origin_coords[2]),
			"{} {} {}".format(self.dipole_moment._vector_coords[0], self.dipole_moment._vector_coords[1], self.dipole_moment._vector_coords[2]),
			"{}".format(self.file_path['x0y0z0'].with_suffix(self.scene_file_extension)),
			"{}".format(self.file_path['x90y0z0'].with_suffix(self.scene_file_extension)),
			"{}".format(self.file_path['x0y90z0'].with_suffix(self.scene_file_extension)),
			"{}".format(self.file_path['x45y45z45'].with_suffix(self.scene_file_extension))
		]