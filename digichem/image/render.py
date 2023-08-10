"""
Classes for rendering 3D images of molecules and densities, primarily with blender and beautiful atoms.

Also see vmd.py for an older render engine.
"""
from pathlib import Path

from silico.file.base import File_converter
from silico.image.base import Cropable_mixin
import pkg_resources

class Render_maker(File_converter, Cropable_mixin):
    """
    ABC for classes that make 3D renders from cube files.
    """
        
    # Text description of our input file type, used for error messages etc. This can be changed by inheriting classes.
    input_file_type = "cube"
    # Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
    output_file_type = "render"
    
    def __init__(
            self,
            *args,
            cube_file = None,
            rotations = None,
            auto_crop = True,
            resolution = 1024,
            also_make_png = True,
            isovalue = 0.2,
            **kwargs):
        """
        Constructor for Image_maker objects.
        
        :param output: The path to write image files to. As more than one image can be created by this class, this name is used as the base name and is modified for each image. The suffix will be honoured (rendered imgaes will use it as a hint for their output format).
        :param cube_file: The path to a cube_file to use to render images.
        :param rotations: A list of tuples of rotations, where the first index in the tuple specifies the axis to rotate about and the second is the angle to rotate (in radians).
        :param auto_crop: If False, images will not have excess white space cropped.
        :param resolution: The max width or height of the rendered images in pixels.
        :param also_make_png: If True, additional images will be rendered in PNG format. This option is useful to generate higher quality images alongside more portable formats.
        :param isovalue: The isovalue to use for rendering isosurfaces. Has no effect when rendering only atoms.
        :param blender_executable:
        """
        super().__init__(*args, input_file = cube_file, **kwargs)
        # Save our translations list.
        #self.translations = translations if translations is not None else (0,0,0)
        # This is not used, why not?
        self.translations = (0,0,0)
        # Save our rotations list.
        self._rotations = rotations if rotations is not None else []
        
        # Some options that control how we function.
        self.auto_crop = auto_crop
        self.target_resolution = resolution
        self.also_make_png = also_make_png
        self.isovalue = isovalue
                
        # These 4 attributes are file paths to the four images we create.
        # We'll keep the same file extension as the was given to us.
        self.file_path = {
            'x0y0z0': self.output.with_suffix(".x0y0z0" + self.output.suffix),
            'x90y0z0': self.output.with_suffix(".x90y0z0" + self.output.suffix),
            'x0y90z0': self.output.with_suffix(".x0y90z0" + self.output.suffix),
            'x45y45z45': self.output.with_suffix(".x45y45z45"  + self.output.suffix),
            # TODO: This should respect also_make_png, but works fine for now (these are actually unused).
            # There are higher quality PNG versions.
             'x0y0z0_big': self.output.with_suffix(".x0y0z0.png"),
             'x90y0z0_big': self.output.with_suffix(".x90y0z0.png"),
             'x0y90z0_big': self.output.with_suffix(".x0y90z0.png"),
             'x45y45z45_big': self.output.with_suffix(".x45y45z45.png")
            }
        
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
    def rotations(self):
        # Silico rotates the wrong way round for some reason, reverse for our rendering engine.
        return [(axis, -theta) for axis, theta in self._rotations]


class Batoms_renderer(Render_maker):
    """
    Class for rendering images with Beautiful Atoms and Blender.
    """
        
    # Text description of our input file type, used for error messages etc. This can be changed by inheriting classes.
    input_file_type = "cube"
    # Text description of our output file type, used for error messages etc. This can be changed by inheriting classes.
    output_file_type = "render"
    
    def __init__(
            self,
            *args,
            cube_file = None,
            rotations = None,
            auto_crop = True,
            resolution = 1024,
            render_samples = 256,
            also_make_png = True,
            isovalue = 0.02,
            blender_executable = None,
            cpus = 1,
            perspective = "perspective",
            **kwargs):
        """
        Constructor for Image_maker objects.
        
        :param output: The path to write image files to. As more than one image can be created by this class, this name is used as the base name and is modified for each image. The suffix will be honoured (rendered imgaes will use it as a hint for their output format).
        :param cube_file: The path to a cube_file to use to render images.
        :param rotations: A list of tuples of rotations, where the first index in the tuple specifies the axis to rotate about and the second is the angle to rotate (in radians).
        :param auto_crop: If False, images will not have excess white space cropped.
        :param resolution: The max width or height of the rendered images in pixels.
        :param also_make_png: If True, additional images will be rendered in PNG format. This option is useful to generate higher quality images alongside more portable formats.
        :param isovalue: The isovalue to use for rendering isosurfaces. Has no effect when rendering only atoms.
        :param blender_executable:
        """
        super().__init__(
            *args,
            cube_file = cube_file,
            rotations = rotations,
            auto_crop = auto_crop,
            resolution = resolution,
            also_make_png = also_make_png,
            isovalue = isovalue,
            **kwargs
        )
        
        # Blender specific options.
        self.render_samples = render_samples
        self.cpus = cpus
        self.perspective = perspective
        
        # Use explicit blender location if given.
        if blender_executable is not None:
            self.blender_executable = blender_executable
        
        else:
            # Otherwise, use a default location.
            self.blender_executable = Path(pkg_resources.resource_filename('silico', 'data/blender/blender'))
            
    @classmethod
    def from_options(self, output, *, cube_file = None, rotations = None, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            cube_file = cube_file,
            rotations = rotations,
            auto_crop = options['rendered_image']['auto_crop'],
            rendering_style = options['rendered_image']['rendering_style'],
            resolution = options['rendered_image']['resolution'],
            isovalue = options['rendered_image'][self.options_name]['isovalue'],
            use_existing = options['rendered_image']['use_existing'],
            dont_modify = not options['rendered_image']['enable_rendering'],
            vmd_executable = options['external']['vmd'],
            tachyon_executable = options['external']['tachyon'],
            vmd_logging = options['logging']['vmd_logging'],
            **kwargs
        )
