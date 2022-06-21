from silico.report.image.gaussian import Gaussian_setup
from silico.report.image.turbomole import Turbomole_setup
from silico.report.image.main import Cube_setup

def class_from_result(result):
    """
    Get an appropriate image setup class based on a metadata object.
    """
    if result.metadata.package == "Turbomole":
        return Turbomole_setup
    elif result.metadata.package == "Gaussian":
        return Gaussian_setup
    else:
        return Cube_setup
    
def from_result(report, *, result, do_orbitals, do_spin, options, calculation = None):
    """
    Constructor for Image_setup objects.
    
    :param report: The report object we will setup images for.
    :param result: The result corresponding to the (sub) calculation we will make images from.
    :param do_orbitals: Whether to generate orbitals.
    :param do_spin: Whether to generate spin density plots.
    :param options: Dictionary of config options.
    :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
    """
    return class_from_result(result)(report, result = result, options = options, do_orbitals = do_orbitals, do_spin = do_spin, calculation = calculation)