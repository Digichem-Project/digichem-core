# Silico logging.
import silico.logging

# General imports.
try:
    import rdkit.Chem.Draw
    import rdkit.Chem.AllChem
    HAVE_RDKIT = True
    
except Exception:
    # Missing rdkit module.
    silico.logging.get_logger().debug("Failed to import module 'rdkit', falling back to openbabel for 2D structure images")
    HAVE_RDKIT = False

# Silico imports.
from silico.image.base import Image_maker


class Skeletal_image_maker(Image_maker):
    """
    A class for rendering skeletal-style molecule structure images.
    """
    
    def __init__(self, output, coords, *args, resolution = 1024, **kwargs):
        """
        Constructor for Structure_image_maker objects.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg: png, jpeg).
        :param coords: The coordinates of the molecule/system to render as a Silico_coords object.
        :param resolution: The width and height of the rendered image.
        """
        self.coords = coords
        self.resolution = resolution
        super().__init__(output, *args, **kwargs)
        
    @classmethod
    def from_options(self, output, *, coords, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """    
        return self(
            output,
            coords = coords,
            **kwargs
        )
        
    def make_files(self):
        """
        Make the image described by this object.
        """
        # There are two libraries we can use for rendering structure images.
        # The first is rdkit, which is generally superior but is only used in this one part of silico,
        # hence it would be unreasonable to demand rdkit as a hard dependency.
        # As an alternative, pybel (from openbabel) can also be used, wich has the advantage that
        # openbabel is used in various parts of silico for file conversion and so is already a hard dependency.
        # However, the images rendered by pybel are typically inferior (especially where sterics are tight),
        # so rdkit will be preferred if it is available.
        if HAVE_RDKIT:
            # We have rdkit.
            # First, convert our geometry to a format that can be used by rdkit.
            molecule = rdkit.Chem.MolFromMolBlock(self.coords.to_format("mol"))
            rdkit.Chem.AllChem.Compute2DCoords(molecule)
            
            # Then write the file.
            rdkit.Chem.Draw.MolToFile(molecule, self.output, (self.resolution, self.resolution))
            
        else:
            # No rdkit, convert directly.
            self.coords.to_format("png", self.output)
            
            
            
