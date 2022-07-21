# Silico logging.
import silico.logging

# General imports.
# rdkit is an optional, hidden import

# Silico imports.
from silico.image.base import Image_maker


class Skeletal_image_maker(Image_maker):
    """
    A class for rendering skeletal-style molecule structure images.
    """
    
    def __init__(self, output, coords, *args, resolution = 1024, render_backend = "rdkit", **kwargs):
        """
        Constructor for Structure_image_maker objects.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg: png, jpeg).
        :param coords: The coordinates of the molecule/system to render as a Silico_coords object.
        :param resolution: The width and height of the rendered image.
        :param render_backend: The library to prefer to use for rendering, possible options are 'rdkit' or 'obabel'. If 'rdkit' is chosen but is not available, obabel will be used as a fallback. 
        """
        self.coords = coords
        self.resolution = resolution
        self.render_backend = render_backend
        assert self.render_backend == "rdkit" or self.render_backend == "obabel"
        super().__init__(output, *args, **kwargs)
        
    @classmethod
    def from_options(self, output, *, coords, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """    
        return self(
            output,
            coords = coords,
            resolution = options['skeletal_image']['resolution'],
            render_backend = options['skeletal_image']['render_backend'],
            **kwargs
        )
        
    def make_files(self):
        """
        Make the image described by this object.
        """
        render_backend = self.render_backend
        if render_backend == "rdkit":
            # Try and import rdkit.
            try:
                import rdkit.Chem.Draw
                import rdkit.Chem.AllChem
                #import rdkit.Chem.rdchem
                import rdkit.RDLogger
                
            except ModuleNotFoundError:
                # Missing rdkit module.
                silico.logging.get_logger().warning("Failed to import module 'rdkit', falling back to openbabel for 2D structure images", exc_info = True)
                render_backend = 'obabel'
            
            except Exception:
                # Something went wrong.
                silico.logging.get_logger().error("An error occurred trying to import module 'rdkit', falling back to openbabel for 2D structure images", exc_info = True)
                render_backend = 'obabel'


        if render_backend == "rdkit":
            # We've been asked to use rdkit.
            # First, convert our geometry to a format that can be used by rdkit.
            # WrapLogs() outputs rdkit logging to python's stderr (which might be redirected to an urwid widget).
            # If/when rdkit is further intergrated into silico, this call will likely be moved elsewhere. 
            #rdkit.Chem.rdchem.WrapLogs()
            # Sadly the behaviour of WrapLogs() is a bit bizzare, although we do get redirection to our custom widgets etc,
            # logs are also still dumped to screen...
            # for now, disable logging...
            rdkit.RDLogger.DisableLog('rdApp.*')
            
            molecule = rdkit.Chem.MolFromMolBlock(self.coords.to_format("mol"))
            # rdkit will silently return None if parsing fails, best to a check.
            if molecule is None:
                raise Exception("Failed to parse coordinates with rdkit")
            rdkit.Chem.AllChem.Compute2DCoords(molecule)
            
            # Then write the file.
            rdkit.Chem.Draw.MolToFile(molecule, str(self.output), (self.resolution, self.resolution))
            
        else:
            # We've been asked to use obabel
            self.coords.to_format("png", self.output)



