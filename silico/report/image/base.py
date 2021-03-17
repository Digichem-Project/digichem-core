from silico.report.image.gaussian import Gaussian_setup
from silico.report.image.turbomole import Turbomole_setup
from silico.report.image.main import Image_setup

def class_from_metadata(metadata):
    """
    Get an appropriate image setup class based on a metadata object.
    """
    if metadata.package == "Turbomole":
        return Turbomole_setup
    elif metadata.package == "Gaussian":
        return Gaussian_setup
    else:
        return Image_setup
    
def from_metadata(report, metadata, options, calculation = None):
    """
    Constructor for Image_setup objects.
    
    :param report: The Report object we'll create images for.
    :param metadata: The metadata object to use to setup images.
    :param options: Options to use to create images.
    :param calculation: Optional calculation which will be used as a template for new calculations to generate new images.
    """
    return class_from_metadata(metadata)(report, metadata, options, calculation = calculation)