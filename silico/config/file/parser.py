#
# Code for loading/reading/parsing config files and configurable files.
#

import pkg_resources
from pathlib import Path
import yaml
import glob
import itertools

from silico.config.base import Config
from silico.config.main import Silico_options
from silico.config.configurable.loader import Partial_loader, Configurable_list,\
    Update_loader
from silico.config.configurable.loader import Single_loader
from silico.exception.configurable import Configurable_loader_exception


class Config_parser():
    """
    Class for loading standard silico options.
    """
    
    @classmethod
    def DEFAULT_CONFIG_PATH(self):
        """
        Location of the config containing default config options.
        """
        return Path(pkg_resources.resource_filename('silico', "data/config/.default.yaml"))

    @classmethod
    def MASTER_CONFIG_PATH(self):
        """
        Location of the master config path.
        """
        return Path(pkg_resources.resource_filename('silico', "data/config/silico.yaml"))
        
    @classmethod
    def SYSTEM_CONFIG_LOCATION(self):
        """
        The location of the system-wide config file (which has precedence over the master).
        """
        return Path("/etc/silico/silico.yaml")
    
    @classmethod
    def USER_CONFIG_LOCATION(self):
        """
        The location of the user specific config file (which has precedence over the system).
        """
        return Path(Path.home(), ".config/silico/silico.yaml") 
    
    
    def __init__(self, config_string):
        """
        Constructor for Options loader objects.
        
        :param config_string: A string to parse options from.
        """
        self.config_string = config_string
        self.config_path = None
    
    def parse(self, config_stream):
        """
        """
        # Load with PyYaml.
        config = yaml.safe_load(config_stream)
        
        # Check the config isn't None (which it will be if the file is empty.)
        try:
            return self.pre_process(config)
        
        except Exception as e:
            if config is None:
                # Don't panic, the file was just empty.
                return self.pre_process({})
            else:
                # Something else went wrong; panic.
                raise e
    
    def load(self):
        """
        Load the config file.
        """
        return self.parse(self.config_string)
    
    def pre_process(self, config):
        """
        Peform some pre-processing on a just-loaded config.
        
        :param config: Dictionary of config options.
        """
        return Config(config, FILE_NAME = self.config_path)
    
        
class Config_file_parser(Config_parser):
    """
    Class for loading config files.
    """
    
    def __init__(self, path):
        """
        Constructor for Options loader objects.
        
        :param path: Path to the config file to load.
        """
        self.config_path = path
        
    def load(self, not_exists_ok = False):
        """
        Load the config file.
        
        :param not_exists_ok: If False and our config file does not exist, raise an exception.
        """
        try:
            with open(self.config_path, "rt") as config_stream:
                return self.parse(config_stream)
                    
        except FileNotFoundError:
            if not_exists_ok:
                # No need to panic.
                return self.pre_process({})
            
            else:
                # Panic.
                raise
    
    @classmethod
    def silico_options(self):
        """
        Load all Silico config files.
        
        Config files are searched for in the following order (increasing precedence) 1) the silico install directory, 2) /etc/silico, 3) the user's ~/.config/silico directory
        
        Note that configurables will be added to the returned Silico_options, but they will not be resolved.
        
        :return: A Silico_options object (a fancy dict).
        """
        # First, we always load the default config object.
        config = self(self.DEFAULT_CONFIG_PATH()).load(True)
        
        # Next the 'master' (located in the silico install directory).
        config.merge(self(self.MASTER_CONFIG_PATH()).load(True))
        # Next the system (located at /etc/silico)
        config.merge(self(self.SYSTEM_CONFIG_LOCATION()).load(True))
        # Finally the user specific (located at ~/.config/silico)
        config.merge(self(self.USER_CONFIG_LOCATION()).load(True))
        
        # Get a merged version.
        options = Silico_options(config)
        
        # Add configurables.
        Configurables_parser.add_default_configurables(options)
        
        # And return.
        return options

    
class Configurables_parser():
    """
    Reads and parsers all configurable files from a location.
    """
    
    def __init__(self, *paths, TYPE):
        """
        Constructor for Configurables_loader object.
        
        :param paths: Paths to  directories to load .yaml files from. All *.yaml files under each directory will be loaded and processed.
        :param TYPE: The TYPE of the configurables we are loading; this is a string which identifies the type of the configurables (eg, Method, Calculation etc).
        """
        self.root_directories = paths
        # A type to set for all configurables we load.
        self.TYPE = TYPE
        # A list of the configurable loaders we have parsed.
        self.loaders = []
        # Configs that update other configs.
        self.updates = []
        
    def load(self):
        """
        Read, parse and process all configs at our location.
        
        :returns: A list of configurable loaders
        """
        self.parse()
        
        return self.process()
    
    
    def parse(self):
        """
        Read and parse all .yaml files within our directory.
        """
        # Get our file list; all files within our top directory which end in .yaml.
        file_list = itertools.chain( *(sorted(glob.glob(glob.escape(root_directory) + "/**/*.yaml", recursive = True)) for root_directory in self.root_directories) )
        
        # Now parse them.
        for file_name in file_list:
            try:
                with open(file_name, "rt") as file:
                    # Parse each.
                    for config in yaml.safe_load_all(file):
                        # If the config has its SUB_TYPE set to update, file it away separately.
                        if config.get('SUB_TYPE') == "update":
                            # An update, no need to pre process.
                            self.updates.append(Update_loader(file_name, self.TYPE, config))
                        
                        else:
                            # Normal config, pre process and add.
                            self.loaders.append(self.pre_process(config, file_name))
                    
            except FileNotFoundError:
                # This should be ok to ignore, just means a file was moved inbetween us looking how many files there are and reading them.
                # Possibly no need to worry about this?
                pass
            
    def process(self):
        """
        Process the config dicts that we have parsed.
        """
        # First, apply any updates.
        for update in self.updates:
            update.update(self.loaders)
        
        # We need to link any partial loaders together.
        for loader in self.loaders:
            loader.link(self.loaders)
            
        # Next we need to purge partial configurables from the top level list.
        #conf_list = Configurable_list([loader for loader in self.loaders if loader.TOP], self.TYPE)
        conf_list = [loader for loader in self.loaders if loader.TOP]
        return conf_list       
    
    def pre_process(self, config, config_path):
        """
        Convert loaded config dicts to appropriate objects
        """
        config['TYPE'] = self.TYPE
        
        # First, panic if no TAG is set.
        if "TAG" not in config:
            raise Configurable_loader_exception(config, self.TYPE, config_path, "missing required option 'TAG'")
        
        # If we have a sub type set, use that to get the name of the class.
        if 'SUB_TYPE' in config:
            if config['SUB_TYPE'] == "pseudo" or config['SUB_TYPE'] == "partial":
                return Partial_loader(config_path, self.TYPE, config, pseudo = config['SUB_TYPE'] == "pseudo")
            
            else:
                # Panic, we don't recognise this sub type.
                raise Configurable_loader_exception(config, self.TYPE, config_path, "SUB_TYPE '{}' is not recognised".format(config['SUB_TYPE']))
            
        else:
            # Use a single loader.
            return Single_loader(config_path, self.TYPE, config)        
        
    
    @classmethod
    def add_default_configurables(self, options):
        """
        Load configurables from all silico locations and add to a silico options object.
        
        :param options: A Silico_options object to populate with configurables.
        """            
        # Iterate through each configurable type.
        for configurable_name, configurable_type, TYPE in (
            ("Basis Sets", "basis_sets", "basis_set"),
            ("Calculations", "calculations", "calculation"),
            ("Programs", "programs", "program"),
            ("Methods", "methods", "method")
        ):
            root_directories = [Path(location, configurable_name) for location in (Config_parser.MASTER_CONFIG_PATH().parent, Config_parser.SYSTEM_CONFIG_LOCATION().parent, Config_parser.USER_CONFIG_LOCATION().parent)]
            
            # Load from the location and add to our silico_options object.
            # The name of the attribute we add to is the same as the location we read from, but in lower case...
            try:
                getattr(options, configurable_type).NEXT.extend(self(*root_directories, TYPE = TYPE).load())
                
            except FileNotFoundError:
                # No need to panic.
                pass
