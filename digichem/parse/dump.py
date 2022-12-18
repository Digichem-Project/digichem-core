"""Parser(s) for reading from dumped silico results sets"""
import yaml
from silico.parse.base import Parser_abc
from silico.result.tdm import Transition_dipole_moment

class Dump_parser(Parser_abc):
    """
    Parser for reading from dumped result sets in yaml format.
    """
    
    def __init__(self, yaml_file):
        """
        Top level constructor for calculation parsers.
        
        :param yaml_file: The path to a yaml file to read.
        """
        super().__init__(yaml_file)
    
    def parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        with open(self.log_file_paths[0], "rt") as yaml_file:
            self.data = yaml.safe_load(yaml_file)
            
    def load_result_part(self, result_cls, *, data = None, **kwargs):
        """
        Get part of a result file.
        
        For most parsers, this will simply call from_parser() of the given class, but some parsers do something more interesting.
        Any arguments other than cls will be parsed to the underlying function.
        """
        data = self.data if data is None else data
                
        if result_cls != Transition_dipole_moment:
            return result_cls.from_dump(data, **kwargs)
        
        else:
            # This is a bit hacky...
            return result_cls.list_from_dump(data, **kwargs)