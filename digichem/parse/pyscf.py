from cclib.bridge.cclib2pyscf import _makecclib

from digichem.parse.base import Parser_abc
import digichem.log
import digichem.file.types as file_types

class Pyscf_parser(Parser_abc):
    """
    Top level class for parsing output from pyscf data.
    """
    
    def __init__(self, methods, **kwargs):
        self.methods = methods
        super().__init__(**kwargs)

    def _parse(self):
        self.data = _makecclib(**self.methods)
        self.data._id = "pass"
            
            