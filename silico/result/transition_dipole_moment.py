# General imports.

# Silico imports.
from silico.result.dipole_moment import Electric_dipole_moment_mixin, Dipole_moment_ABC, Magnetic_dipole_moment_mixin

class Transition_dipole_moment(Dipole_moment_ABC):
    """
    ABC for electric and magnetic transition dipole moments.
    """
    
    def __init__(self, state_level, *args, **kwargs):
        """
        Electric_transition_dipole_moment constructor.
        
        :param state: The excited state level this dipole transitions to.
        """
        super().__init__(*args, **kwargs)
        
        # An int describing the excited state we belong to. This is always available.
        self.state_level = state_level
        
        # The excited state object that we belong to. This may be None. We can't set this here because then both the Excited_state and Electric_transition_dipole_moment constructors would depend on each other.
        self.excited_state = None
        
        # Save a name describing which dipole we are (permanent vs transition etc).
        self.dipole_type = "transition"
    
    @property
    def name(self):
        """
        Name that describes this dipole moment.
        """
        return "{} TDM".format(self.excited_state.state_symbol)
    
    def set_excited_state(self, excited_state):
        """
        Set the excited state object that we belong to.
        """
        # Not sure what we should do if the given excited_state has a different level to what we expect.
        if excited_state.level != self.state_level:
            pass
        self.excited_state = excited_state
        
        
class Electric_transition_dipole_moment(Transition_dipole_moment, Electric_dipole_moment_mixin):
    """
    Class that represents an electric transition dipole moment.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.electromagnetic_type = "electric"
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of TDMs from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of TDM objects. An empty list will be returned if no TDM data is available.
        """
        try:
            return [
                self(state_level = index +1,
                origin_coords = (0,0,0),
                vector_coords = (parser.au_to_debye(tdm[0]), parser.au_to_debye(tdm[1]), parser.au_to_debye(tdm[2])),
                atoms = parser.results.alignment) for index, tdm in enumerate(parser.data.etdips)]
        except AttributeError:
            return []


class Magnetic_transition_dipole_moment(Transition_dipole_moment, Magnetic_dipole_moment_mixin):
    """
    Class that represents a magnetic transition dipole moment.
    
    Unlike electric dipole moments, the units here are in a.u. (0.5 bohr magnetons).
    """
        
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.electromagnetic_type = "magnetic"
    
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of TDMs from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of TDM objects. An empty list will be returned if no TDM data is available.
        """
        try:
            return [
                self(state_level = index +1,
                origin_coords = (0,0,0),
                vector_coords = tdm,
                atoms = parser.results.alignment) for index, tdm in enumerate(parser.data.etmagdips)]
        except AttributeError:
            return []

