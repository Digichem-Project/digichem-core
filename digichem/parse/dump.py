"""Parser(s) for reading from dumped digichem results sets"""
import yaml
import json

from digichem.parse.base import Parser_abc
from digichem.result.tdm import Transition_dipole_moment
from digichem.result.result import Result_set
from digichem.result.metadata import Metadata
from digichem.result.orbital import Molecular_orbital_list
from digichem.result.atom import Atom_list
from digichem.result.excited_state import Excited_state_list
from digichem.result.energy import Energies
from digichem.result.ground_state import Ground_state
from digichem.result.soc import SOC_list
from digichem.result.dipole_moment import Dipole_moment
from digichem.result.vibration import Vibrations_list
from digichem.result.emission import Relaxed_excited_state
import digichem.log
from digichem.result.alignment.base import Alignment, Minimal
from digichem.result.nmr import NMR_list


class Dump_multi_parser_abc(Parser_abc):
    """
    ABC for classes that can read multiple result sets from dumped data.
    """
    
    def __init__(self, input_file, *other_log_files, raw_data = None, **auxiliary_files):
        """
        Top level constructor for calculation parsers.
        
        :param input_file: The path to read from.
        """
        # Neither other_log_files nor auxiliary_files are currently used...
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
    
    def post_parse(self):
        """
        Perform any required operations after line-by-line parsing.
        """
        for result in self.data:
            # Add current username.
            # TODO: It would probably be better if we used the name of the user who owns the output file, rather than the current user...
            result['metadata']['user'] = self.get_current_username()
            
            # Set our file paths.
            result['metadata']['log_files'] = self.log_file_paths
            result['metadata']['auxiliary_files'] = {}
    
    def get_sub_parser(self):
        raise NotImplementedError("Implement in subclass")
            
    def process_all(self, options):
        """
        Get all the Result set objects produced by this parser.
        
        :param options: A Digichem options nested dictionary containing options to control parsing.
        :return: A list of the populated result sets.
        """
        # Unlike most other parsers, our data can actually contain lots of results.
        self.all_results = [self.get_sub_parser()(self.log_file_path, raw_data = data).process(options) for data in self.data]
        
        return self.all_results
    
    def process(self, options):
        """
        Get a Result set object from this parser.
        
        :param options: A Digichem options nested dictionary containing options to control parsing.
        :return: The populated result set.
        """
        self.process_all(options)
        return self.results
    
    
class Yaml_multi_parser(Dump_multi_parser_abc):
    """
    Parser for reading from dumped result sets in yaml format.
    """
    
    def get_sub_parser(self):
        return Yaml_parser
    
    def _parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        digichem.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        with open(self.log_file_path, "rt") as yaml_file:
            self.data = list(yaml.safe_load_all(yaml_file))
    

class Json_multi_parser(Dump_multi_parser_abc):
    """
    Parser for reading from dumped result sets in JSON format.
    """
    
    def get_sub_parser(self):
        return Json_parser
    
    def _parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        digichem.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        with open(self.log_file_path, "rt") as json_file:
            self.data = json.load(json_file)


class Dump_parser_abc(Parser_abc):
    """
    ABC for parsers that read dumped data.
    """
    
    def __init__(self, input_file, *other_log_files, raw_data = None, **auxiliary_files):
        """
        Top level constructor for calculation parsers.
        
        :param input_file: The path to the input file to read.
        """
        # Neither other_log_files nor auxiliary_files are currently used...
        super().__init__(input_file, raw_data = raw_data)
        
    @classmethod
    def from_data(self, input_file, data):
        return self(input_file, raw_data = data)
    
    def process_all(self, options):
        """
        Get all the Result set objects produced by this parser.
        
        :param options: A Digichem options nested dictionary containing options to control parsing.
        :return: A list of the populated result sets.
        """
        self.process(options)
        return [self.results]
            
    def process(self, options):
        """
        Get a Result set object from this parser.
        
        :param options: A Digichem options nested dictionary containing options to control parsing.
        :return: The populated result set.
        """
        
        # Get our result set.
        self.results = Result_set(
            _id = self.data.get("_id"),
            metadata = Metadata.from_dump(self.data['metadata'], self.results, options)
        )
        
        # First get our list of MOs (because we need them for excited states too.)
        self.results.orbitals = Molecular_orbital_list.from_dump(self.data['orbitals'], self.results, options)
        self.results.beta_orbitals = Molecular_orbital_list.from_dump(self.data['beta_orbitals'], self.results, options)
        
        alignment_class = Alignment.from_class_handle(options['alignment']) if options['alignment'] is not None else Minimal
        
        # Our alignment orientation data.
        # The constructor for each alignment class automatically performs realignment.
        self.results.atoms = alignment_class.from_dump(self.data['atoms'], self.results, options)
        self.results.raw_atoms = Atom_list.from_dump(self.data['raw_atoms'], self.results, options)
        
        # TEDM and TMDM.
        self.results.transition_dipole_moments = Transition_dipole_moment.list_from_dump(self.data['excited_states']['values'], self.results, options)
        
        # Excited states.
        self.results.excited_states = Excited_state_list.from_dump(self.data['excited_states'], self.results, options)
        
        # Energies.
        self.results.energies = Energies.from_dump(self.data['energies'], self.results, options)
        
        # Our ground state.
        self.results.ground_state = Ground_state.from_dump(self.data['ground_state'], self.results, options)
        
        # And a similar list but also including the ground.
        self.results.energy_states = Excited_state_list()
        self.results.energy_states.append(self.results.ground_state)
        self.results.energy_states.extend(self.results.excited_states)
        
        # SOC.
        self.results.soc = SOC_list.from_dump(self.data['soc'], self.results, options)
        
        # PDM
        self.results.pdm = Dipole_moment.from_dump(self.data['pdm'], self.results, options)
        
        # Frequencies.
        self.results.vibrations = Vibrations_list.from_dump(self.data['vibrations'], self.results, options)
        
        # NMR.
        if "nmr" in self.data:
            self.results.nmr = NMR_list.from_dump(self.data['nmr'], self.results, options)
        
        else:
            self.results.nmr = NMR_list(atoms = self.results.atoms, options = options)
            
        # Finally, try and set emission.
        try:
            self.results.emission = self.results.emission.from_dump(self.data['emission'], self.results, options)
        
        except Exception:
            digichem.log.get_logger().warning("Failed to parse emission data", exc_info = True)
        
        # If we don't have explicit emission set, try and guess.
        #if len(self.results.emission.vertical) == 0 and len(self.results.emission.adiabatic) == 0:
        #    self.results.emission.vertical, self.results.emission.adiabatic = Relaxed_excited_state.guess_from_results(self.results)
        
        # Return the populated result set for convenience.
        return self.results


class Yaml_parser(Dump_parser_abc):
    """
    Parser for reading from dumped result sets in yaml format.
    """
    
    def _parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        digichem.log.get_logger().info("Parsing calculation result '{}'".format(self.description))
        with open(self.log_file_paths[0], "rt") as yaml_file:
            self.data = yaml.safe_load(yaml_file)
            

class Json_parser(Dump_parser_abc):
    """
    Parser for reading from TinyDB json files.
    """
    
    def _parse(self):
        """
        Extract results from our output files.
        """
        raise NotImplementedError("Databases contain many results so it does not make sense to return only one. Try Tinydb_multi_parser instead.")
                
        
    
