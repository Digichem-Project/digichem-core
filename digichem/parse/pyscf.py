from cclib.bridge.cclib2pyscf import _makecclib

from digichem.parse.base import Parser_abc
import digichem.log
import digichem.file.types as file_types

class Pyscf_parser(Parser_abc):
    """
    Top level class for parsing output from pyscf data.
    """
    
    def __init__(self, mol_name, methods, **kwargs):
        self.methods = methods
        self.mol_name = mol_name
        super().__init__(**kwargs)

    def _parse(self):
        self.data = _makecclib(**self.methods)
        self.data._id = "pass"
        self.data.metadata['name'] = self.mol_name
        self.data._aux = {'methods': self.methods}
            
            