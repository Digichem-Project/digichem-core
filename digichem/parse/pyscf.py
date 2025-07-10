import hashlib
import json
from uuid import uuid4

from cclib.bridge.cclib2pyscf import cclibfrommethods

from digichem.parse.base import Parser_abc
import digichem.log

class Pyscf_parser(Parser_abc):
    """
    Top level class for parsing output from pyscf data.
    """
    
    def __init__(self, mol_name, methods, **kwargs):
        self.methods = methods
        self.mol_name = mol_name
        super().__init__(**kwargs)

    def _parse(self):
        self.data = cclibfrommethods(**self.methods)
        
        try:
            # Try to generate a checksum from metadata.
            self.data._id = hashlib.sha1(json.dumps(self.data.metadata, sort_keys = True).encode('utf-8')).hexdigest()

        except Exception:
            # No luck, something in metadata must be unhashable.
            digichem.log.get_logger().error("Unable to generate hash ID from calculation metadata, using random ID instead", exc_info = True)
            # TODO: Think of a better way to do this.
            self.data._id = hashlib.sha1(uuid4().hex.encode('utf-8')).hexdigest()
        
        self.data.metadata['name'] = self.mol_name
        self.data._aux = {'methods': self.methods}
            