from pysoc.io.soc_td import SOC


class Spin_orbit_coupling(SOC):
    """
    Class that represents spin-orbit coupling between two states.
    """
    
    def __init__(self, singlet_state, triplet_state, positive_one, zero, negative_one):
        """
        Constructor for SOC objects.
        
        :param singlet_state: The singlet excited state this coupling is between.
        :param triplet_state: The triplet excited state this coupling is between.
        :param positive_one: SOC with quantum number +1.
        :param zero: SOC with quantum number 0.
        :param negative_one: SOC with quantum number -1.
        """
        self.singlet_state = singlet_state
        self.triplet_state = triplet_state
        self.positive_one = positive_one
        self.zero = zero
        self.negative_one = negative_one
        
    @classmethod
    def list_from_parser(self, parser):
        """
        Create a SOC object from an output file parser.
        """
        # Go through our data.
        try:
            return [self(
                singlet_state = parser.results.energy_states.get_state(state_symbol = singlet_symbol),
                triplet_state = parser.results.energy_states.get_state(state_symbol = triplet_symbol),
                positive_one = positive_soc,
                zero = zero_soc,
                negative_one = negative_soc
            ) for (singlet_symbol, triplet_symbol), (positive_soc, zero_soc, negative_soc) in zip(parser.data.socstates, parser.data.soc)]
        except AttributeError:
            # No data.
            return [] 
        