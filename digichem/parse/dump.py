"""Parser(s) for reading from dumped silico results sets"""
import yaml
from silico.parse.base import Parser_abc

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
            
    def load_result_part(self, result_cls, *args, **kwargs):
        """
        Get part of a result file.
        
        For most parsers, this will simply call from_parser() of the given class, but some parsers do something more interesting.
        Any arguments other than cls will be parsed to the underlying function.
        """
        return result_cls.from_dump(self.data, *args, **kwargs)