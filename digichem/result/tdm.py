# Silico imports.
from silico.result.dipole_moment import Electric_dipole_moment_mixin, Dipole_moment_ABC, Magnetic_dipole_moment_mixin
import itertools


class Transition_dipole_moment_ABC(Dipole_moment_ABC):
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
        
        # The excited state object that we belong to. This may be None. We can't set this here because then both the Excited_state and Transition_dipole_moment constructors would depend on each other.
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
        
        
class Electric_transition_dipole_moment(Transition_dipole_moment_ABC, Electric_dipole_moment_mixin):
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


class Magnetic_transition_dipole_moment(Transition_dipole_moment_ABC, Magnetic_dipole_moment_mixin):
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


class Transition_dipole_moment():
    """
    A compound class that represents both the electric and magnetic components of a transition dipole moment.
    
    This class can also be used in most places where an electric TDM is expected, as it will pass references to TEDM attributes to the actual TEDM class.
    """
    
    def __init__(self, electric = None, magnetic = None):
        """
        Constructor for TDM.
        
        Either or both of the substituent TDM can be None.
        """
        self.electric = electric
        self.magnetic = magnetic
        self.electromagnetic_type = "electromagnetic"
    
    def set_excited_state(self, excited_state):
        """
        Set the excited state object that we belong to.
        """
        if self.electric is not None:
            self.electric.set_excited_state(excited_state)
        
        if self.magnetic is not None:
            self.magnetic.set_excited_state(excited_state)
    
    @classmethod
    def list_from_parser(self, parser):
        """
        Get a list of Transition_dipole_moment from an output file parser.
        
        :param parser: An output file parser.
        :return: A list of Transition_dipole_moment objects. An empty list will be returned if no TDM data is available.
        """
        electric = Electric_transition_dipole_moment.list_from_parser(parser)
        magnetic = Magnetic_transition_dipole_moment.list_from_parser(parser)
        
        return [self(electric_part, magnetic_part) for electric_part, magnetic_part in itertools.zip_longest(electric, magnetic)]
        
    def __getattr__(self, name):
        # This check prevents infinite recursion when pickling, and is probably a sensible check anyway...
        if name != "electric":
            return getattr(self.electric, name)
        else:
            raise AttributeError(name)
            
    
    def angle(self, cgs = True):
        """
        Return the angle between the electric and magnetic components of this transition dipole moment.
        """
        return self.electric.angle(self.magnetic, cgs)
    
    def cos_angle(self, cgs = True):
        """
        Return the cos of the angle between the electric and magnetic components of this transition dipole moment.
        """
        return self.electric.cos_angle(self.magnetic, cgs)
    
    @property
    def g_value(self):
        """
        Return the 'g value'; the dissymmerty factor of this transition dipole moment.
        """
        try:
            return (4 * self.electric.gaussian_cgs * self.magnetic.gaussian_cgs * self.cos_angle(True)) / (self.electric.gaussian_cgs **2 + self.magnetic.gaussian_cgs **2)
        except (FloatingPointError, ZeroDivisionError):
            return 0