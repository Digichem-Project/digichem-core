"""Parser(s) for reading from dumped silico results sets"""
import yaml

from silico.parse.base import Parser_abc
from silico.result.tdm import Transition_dipole_moment
from silico.result.result import Result_set
from silico.result.metadata import Metadata
from silico.result.orbital import Molecular_orbital_list

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
        # TODO: URGENT: Should we re-align the alignment class?
        self.results.alignment = alignment_class.from_dump(self.data['atoms'])
        self.results.atoms = self.load_result_part(Atom_list)
        
        # TEDM and TMDM.
        self.results.transition_dipole_moments = self.load_result_part(Transition_dipole_moment)
        
        # Excited states.
        self.results.excited_states = self.load_result_part(Excited_state_list)
        
        # Energies.
        self.results.energies = self.load_result_part(Energies)
        
        # Our ground state.
        self.results.ground_state = self.load_result_part(Ground_state)
        
        # And a similar list but also including the ground.
        self.results.energy_states = Excited_state_list()
        self.results.energy_states.append(self.results.ground_state)
        self.results.energy_states.extend(self.results.excited_states)
        
        # SOC.
        self.results.soc = self.load_result_part(SOC_list)
        
        # PDM
        self.results.pdm = self.load_result_part(Dipole_moment)
        
        # Frequencies.
        self.results.vibrations = self.load_result_part(Vibrations_list)
        
        # Finally, try and set emission.
        self.results.emission.vertical, self.results.emission.adiabatic = Relaxed_excited_state.guess_from_results(self.results)
        
        # Return the populated result set for convenience.
        return self.results
            
    def load_result_part(self, result_cls, *, data = None, **kwargs):
        """
        Get part of a result file.
        
        For most parsers, this will simply call from_parser() of the given class, but some parsers do something more interesting.
        Any arguments other than cls will be parsed to the underlying function.
        """
        data = self.data if data is None else data
                
        if result_cls != Transition_dipole_moment:
            return result_cls.from_dump(data, self.results, **kwargs)
        
        else:
            # This is a bit hacky...
            return result_cls.list_from_dump(data, **kwargs)