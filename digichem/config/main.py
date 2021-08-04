# General imports.
from copy import deepcopy

# Silico imports.
from silico.config.base import Config
from silico.config.configurable.loader import Configurable_list


class Silico_options(Config):
    """
    Class for holding main Silico options from various sources.
    """
    
    def __init__(self, options, *, methods = None, programs = None, calculations = None, basis_sets = None):
        """
        Constructor for Silico_options objects.
        """
        super().__init__(options)
            
        # Set Configurable lists.
        self.methods = Configurable_list([], "method") if methods is None else methods
        self.programs = Configurable_list([], "program") if programs is None else programs
        self.calculations = Configurable_list([], "calculation") if calculations is None else calculations
        self.basis_sets = Configurable_list([], "basis_set") if basis_sets is None else basis_sets
        
    @property
    def effective_core_potentials(self):
        """
        A list of known ECPs.
        """
        return [basis_set for basis_set in self.basis_sets if basis_set.ECP is not None]
            
    def validate(self):
        """
        Validate the configurable list that we contain.
        """
        self.methods.validate()
        self.programs.validate()
        self.calculations.validate()
        self.basis_sets.validate()
    
    def set_log_level(self, logger):
        """
        Set the logging level of a logger based on the config options in the object.
        
        :param logger: The logger to set (from logging.getLogger()).
        """        
        # Set from log_level first.
        if self['logging']['log_level'] == "OFF":
            logger.setLevel(60)
        else:
            logger.setLevel(self['logging']['log_level'])
        
        # Now adjust with verbosity.
        verbose = self['logging']['verbose']
        if verbose is not None:
            # Set from verbosity.
            new_level = logger.level - verbose * 10
            
            # Don't allow us to reach 0 (because this is actually 'UNSET').
            if new_level <= 0:
                new_level = 10
            
            # And set.
            logger.setLevel(new_level)
            
    def __deepcopy__(self, memo):
        """
        Overriding deep copy.
        
        Because Silico_options can contain large lists of Configurables (complex objects), real deepcopy can be very slow.
        We overcome this by excluding these lists.
        """
        # TODO: Might be a better way to do this...
        methods = self.methods
        self.methods = []
        programs = self.programs
        self.programs = []
        calculations = self.calculations
        self.calculations = []
        basis_sets = self.basis_sets
        self.basis_sets = []
        # TMP remove __deepcopy__ so we don't recurse.
        copyfunc = self.__deepcopy__
        self.__deepcopy__ = None
        
        new = deepcopy(self, memo)

        
        # Restore.
        self.methods = methods
        new.methods = methods
        self.programs = programs
        new.programs = programs
        self.calculations = calculations
        new.calculations = calculations
        self.basis_sets = basis_sets
        new.basis_sets = basis_sets
        self.__deepcopy__ = copyfunc
        
        # And return.
        return new
            