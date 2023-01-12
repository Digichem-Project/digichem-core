"""Parser(s) for reading from dumped silico results sets"""
import yaml

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

class Dump_parser(Parser_abc):
    """
    Parser for reading from dumped result sets in yaml format.
    """
    
    def __init__(self, yaml_file, *other_log_files, raw_data = None, **aux_files):
        """
        Top level constructor for calculation parsers.
        
        :param yaml_file: The path to a yaml file to read.
        """
        # Neither other_log_files nor aux_files are currently used...
        super().__init__(yaml_file, raw_data = raw_data)
        
    @classmethod
    def from_data(self, yaml_file, data):
        return self(yaml_file, raw_data = data)
    
    def parse(self):
        """
        Extract results from our output files.
        """
        # Read that file.
        with open(self.log_file_paths[0], "rt") as yaml_file:
            self.data = yaml.safe_load(yaml_file)
            
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