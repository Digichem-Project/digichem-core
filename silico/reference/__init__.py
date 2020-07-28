from .base import RSC_format
import pkg_resources

silico_references = {}

# Load our references from file.
with open(pkg_resources.resource_filename('silico', 'data/references.bib'), "rt") as reference_file:
	silico_references = RSC_format.dict_from_file(reference_file)