# General imports.

# Silico imports.
from silico.parser.main import Parser


class Turbomole_parser(Parser):
    """
    Top level class for parsing output from Turbomole files.
    """
    
        
    def parse_metadata(self):
        """
        Parse additional calculation metadata.
        """
        # Add name.
        if self.name is not None:
            self.data.metadata['name'] = self.name
                

    