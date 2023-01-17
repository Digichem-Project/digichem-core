"""Parser(s) for reading from dumped silico results sets"""
import yaml
import json

from silico.parse.base import Parser_abc
from silico.result.tdm import Transition_dipole_moment
from silico.result.result import Result_set
from silico.result.metadata import Metadata
from silico.result.orbital import Molecular_orbital_list
from silico.result.atom import Atom_list
from silico.result.excited_state import Excited_state_list
from silico.result.energy import Energies
from silico.result.ground_state import Ground_state
from silico.result.soc import SOC_list
from silico.result.dipole_moment import Dipole_moment
from silico.result.vibration import Vibrations_list
from silico.result.emission import Relaxed_excited_state
import silico.log


class Dump_multi_parser_abc(Parser_abc):
    """
    ABC for classes that can read multiple result sets from dumped data.
    """
    
    def __init__(self, input_file, *other_log_files, raw_data = None, **aux_files):
        """
        Top level constructor for calculation parsers.
        
        :param input_file: The path to read from.
        """
        # Neither other_log_files nor aux_files are currently used...
        super().__init__(input_file, raw_data = raw_data)
        self.all_results = []
        
    @classmethod
    def from_data(self, input_file, data):
        return self(input_file, raw_data = data)
            
    @property
    def results(self):
        return self.all_results[0]
    
    @results.setter
    def results(self, value):
        pass
    
    def get_sub_parser(self):
        raise NotImplementedError("Implement in subclass")
            
    def process_all(self, alignment_class):
        """
        Get all the Result set objects produced by this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms. If not specified the Minimal alignment method will be used by default.
        :return: A list of the populated result sets.
        """
        # Unlike most other parsers, our data can actually contain lots of results.
        self.all_results = [self.get_sub_parser()(self.log_file_path, raw_data = data).process(alignment_class) for data in self.data]
        
        return self.all_results
    
    def process(self, alignment_class):
        """
        Get a Result set object from this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms. If not specified the Minimal alignment method will be used by default.
        :return: The populated result set.
        """
        self.process_all(alignment_class)
        return self.results[0]
    
    
class Yaml_multi_parser(Dump_multi_parser_abc):
    """
    Parser for reading from dumped result sets in yaml format.
    """
    
    def get_sub_parser(self):
        return Yaml_parser
    
    def parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        silico.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        with open(self.log_file_path, "rt") as yaml_file:
            self.data = list(yaml.safe_load_all(yaml_file))
    

class Json_multi_parser(Dump_multi_parser_abc):
    """
    Parser for reading from dumped result sets in JSON format.
    """
    
    def get_sub_parser(self):
        return Json_parser
    
    def parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        silico.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        with open(self.log_file_path, "rt") as json_file:
            self.data = json.load(json_file)


class Dump_parser_abc(Parser_abc):
    """
    ABC for parsers that read dumped data.
    """
    
    def __init__(self, input_file, *other_log_files, raw_data = None, **aux_files):
        """
        Top level constructor for calculation parsers.
        
        :param input_file: The path to the input file to read.
        """
        # Neither other_log_files nor aux_files are currently used...
        super().__init__(input_file, raw_data = raw_data)
        
    @classmethod
    def from_data(self, input_file, data):
        return self(input_file, raw_data = data)
    
    def process_all(self, alignment_class):
        """
        Get all the Result set objects produced by this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms. If not specified the Minimal alignment method will be used by default.
        :return: A list of the populated result sets.
        """
        self.process(alignment_class)
        return [self.results]
            
    def process(self, alignment_class):
        """
        Get a Result set object from this parser.
        
        :param: alignment_class: An alignment class object to use to reorientate atoms. If not specified the Minimal alignment method will be used by default.
        :return: The populated result set.
        """
        
        # Get our result set.
        self.results = Result_set(Metadata.from_dump(self.data['metadata'], self.results))
        
        # First get our list of MOs (because we need them for excited states too.)
        self.results.orbitals = Molecular_orbital_list.from_dump(self.data['orbitals'], self.results)
        self.results.beta_orbitals = Molecular_orbital_list.from_dump(self.data['beta_orbitals'], self.results)
        
        # Assign total levels.
        #self.results.orbitals.assign_total_level(self.results.beta_orbitals)
        
        # Our alignment orientation data.
        # The constructor for each alignment class automatically performs realignment.
        self.results.atoms = alignment_class.from_dump(self.data['atoms'], self.results)
        self.results.raw_atoms = Atom_list.from_dump(self.data['raw_atoms'], self.results)
        
        # TEDM and TMDM.
        self.results.transition_dipole_moments = Transition_dipole_moment.list_from_dump(self.data['excited_states']['values'], self.results)
        
        # Excited states.
        self.results.excited_states = Excited_state_list.from_dump(self.data['excited_states'], self.results)
        
        # Energies.
        self.results.energies = Energies.from_dump(self.data['energies'], self.results)
        
        # Our ground state.
        self.results.ground_state = Ground_state.from_dump(self.data['ground_state'], self.results)
        
        # And a similar list but also including the ground.
        self.results.energy_states = Excited_state_list()
        self.results.energy_states.append(self.results.ground_state)
        self.results.energy_states.extend(self.results.excited_states)
        
        # SOC.
        self.results.soc = SOC_list.from_dump(self.data['soc'], self.results)
        
        # PDM
        self.results.pdm = Dipole_moment.from_dump(self.data['pdm'], self.results)
        
        # Frequencies.
        self.results.vibrations = Vibrations_list.from_dump(self.data['vibrations'], self.results)
        
        # Finally, try and set emission.
        self.results.emission.vertical, self.results.emission.adiabatic = Relaxed_excited_state.guess_from_results(self.results)
        
        # Return the populated result set for convenience.
        return self.results


class Yaml_parser(Dump_parser_abc):
    """
    Parser for reading from dumped result sets in yaml format.
    """
    
    def parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        silico.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        with open(self.log_file_paths[0], "rt") as yaml_file:
            self.data = yaml.safe_load(yaml_file)
            

class Json_parser(Dump_parser_abc):
    """
    Parser for reading from TinyDB json files.
    """
    
    def parse(self):
        """
        Extract results from our output files.
        """
        raise NotImplementedError("Databases contain many results so it does not make sense to return only one. Try Tinydb_multi_parser instead.")
                
        
    
