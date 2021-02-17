
from silico.config.base import Config
from silico.config.configurable import Configurable_list
from silico.exception.configurable import Configurable_exception
from copy import deepcopy

class Silico_options(Config):
    """
    Class for holding main Silico options from various sources.
    """
    
    def __init__(self, configs = []):
        """
        Constructor for Silico_options objects.
        """
        super().__init__()
            
        # Set Configurable lists.
        self.methods = Configurable_list()
        self.programs = Configurable_list()
        self.calculations = Configurable_list()
        self.basis_sets = Configurable_list()
            
        # Add our list.
        for config in configs:
            self.add_config(config)
            
    def resolve(self):
        """
        Resolve all Configurables which we contain.
        
        NOTE: This method can be (and is) quite expensive.
        """
        # First, resolve each list.
        self.methods = self.methods.resolve()
        self.programs = self.programs.resolve()
        self.calculations = self.calculations.resolve()
        self.basis_sets = self.basis_sets.resolve()
        
        # Now init.
        self.methods.configure()
        self.programs.configure()
        self.basis_sets.configure()
        self.calculations.configure(silico_options = self, available_basis_sets = self.basis_sets)

    def add_config(self, config):
        """
        Add a set of options to this config object.
        
        The functioning of this method is controlled by the 'TYPE' key in config. If this key is None (or is missing entirely), then config is taken to contain normal config options which are added to this object so that later options will override earlier ones.
        If 'TYPE' is not 'None' or missing, then 'TYPE' specifies the name of a key in this Config object to which the given config is appended as a whole object.
        
        :param config: The new config to add; a dict of options. For convenience, a list of dicts can also be given, and they will be added in order. Further, higher-order lists (lists of lists of dicts etc) can also be given and will be fully traversed and added in order. 
        :return: self, for convenience.
        """
        # If config is a list, we'll recursively call ourself.
        if isinstance(config, list):
            for sub_config in config:
                self.add_config(sub_config)
            # Done.
            return
        
        # Either merge or append, depending on what sort of config we have.
        if not hasattr(config, 'TYPE'):
            # Normal config, merge with ourself.
            #self.merge_configs(config, self)
            self.merge(config, none_to_old = True)
        else:
            # Configurable, add to the correct list.
            try:
                getattr(self, config.TYPE + "s").append(config)
            except AttributeError:
                # Possibly because TYPE is invalid; see if we have an attribute called TYPE.
                if hasattr(self, config.TYPE + "s"):
                    # Something else went wrong, don't handle it here.
                    raise
                else:
                    # Bad TYPE.
                    raise Configurable_exception(config, "configurable TYPE '{}' is not recognised".format(config.TYPE)) from None
        
        # Give our self back.
        return self
    
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
            