from silico.report.image.gaussian import Gaussian_setup
from silico.report.image.turbomole import Turbomole_setup
from silico.report.image.main import Cube_setup

def class_from_metadata(metadata):
    """
    Get an appropriate image setup class based on a metadata object.
    """
    if metadata.package == "Turbomole":
        return Turbomole_setup
    elif metadata.package == "Gaussian":
        return Gaussian_setup
    else:
        return Cube_setup
    
def from_metadata(report, *, metadata, do_orbitals, do_spin, options, calculation = None):
    """
    Constructor for Image_setup objects.
    
    :param report: The report object we will setup images for.
    :param metadata: The metadata corresponding to the (sub) calculation we will make images from.
    :param do_orbitals: Whether to generate orbitals.
    :param do_spin: Whether to generate spin density plots.
    :param options: Dictionary of config options.
    :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
    """
    return class_from_metadata(metadata)(report, metadata = metadata, options = options, do_orbitals = do_orbitals, do_spin = do_spin, calculation = calculation)