#
# Code for loading/reading/parsing config files and configurable files.
#

import pkg_resources
from pathlib import Path
import yaml
import os
import himl
import glob

from silico.config.base import Config
from silico.config.main import Silico_options
from silico.config.configurable import Configurable
import silico.misc.directory
from silico.exception.base import Silico_exception
from silico.config.configurable.loader import Partial_loader, Configurable_list
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
        
    def load(self, not_exists_ok = True):
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
        config = self(self.DEFAULT_CONFIG_PATH()).load()
        
        # Next the 'master' (located in the silico install directory).
        config.merge(self(self.MASTER_CONFIG_PATH()).load())
        # Next the system (located at /etc/silico)
        config.merge(self(self.SYSTEM_CONFIG_LOCATION()).load())
        # Finally the user specific (located at ~/.config/silico)
        config.merge(self(self.USER_CONFIG_LOCATION()).load())
        
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
    
    def __init__(self, path, TYPE):
        """
        Constructor for Configurables_loader object.
        
        :param path: Path to a directory to load .yaml files from. All *.yaml files under this directory will be loaded and processed.
        :param TYPE: The TYPE of the configurables we are loading; this is a string which identifies the type of the configurables (eg, Method, Calculation etc).
        """
        self.root_directory = path
        # A type to set for all configurables we load (if not None). Individual configurables can also set this themselves.
        self.TYPE = TYPE
        # A list of the configurable loaders we have parsed.
        self.loaders = []
        
    def load(self):
        """
        Read, parse and process all configs at our location.
        
        :returns: A Configurable_list of loaded configs.
        """
        self.parse()
        
        return self.process()
    
    
    def parse(self):
        """
        Read and parse all .yaml files within our directory.
        """
        # Get our file list; all files within our top directory which end in .yaml.
        file_list = sorted(glob.glob(glob.escape(self.root_directory) + "/**/*.yaml", recursive = True))
        
        # Now parse them.
        for file_name in file_list:
            try:
                with open(file_name, "rt") as file:
                    self.loaders.extend(self.pre_process(config, file_name) for config in yaml.safe_load_all(file))
                    
            except FileNotFoundError:
                # This should be ok to ignore, just means a file was moved inbetween us looking how many files there are and reading them.
                # Possibly no need to worry about this?
                pass
            
    def process(self):
        """
        Process the config dicts that we have parsed.
        """
        # We need to link any partial loaders together.
        for loader in self.loaders:
            loader.link(self.loaders)
            
        # Next we need to purge partial configurables from the top level list.
        conf_list = Configurable_list([loader for loader in self.loaders if loader.TOP], self.TYPE)
        return conf_list
            
        
    def pre_process(self, config, config_path):
        """
        Convert loaded config dicts to appropriate objects
        """
        config['TYPE'] = self.TYPE
        
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
        # Load each of our locations in turn.
        for location in (Config_parser.MASTER_CONFIG_PATH().parent, Config_parser.SYSTEM_CONFIG_LOCATION().parent, Config_parser.USER_CONFIG_LOCATION().parent):
            # And each configurable type:
            for configurable_name, configurable_type, TYPE in (
                ("Basis Sets", "basis_sets", "basis_set"),
                ("Calculations", "calculations", "calculation"),
                ("Programs", "programs", "program"),
                ("Methods", "methods", "method")
            ):
                root_directory = Path(location, configurable_name)
                
                # Load from the location and add to our silico_options object.
                # The name of the attribute we add to is the same as the location we read from, but in lower case...
                try:
                    getattr(options, configurable_type).NEXT.append(self(root_directory, TYPE = TYPE).load())
                    
                except FileNotFoundError:
                    # No need to panic.
                    pass
        
        
        
    
    
    

class Hierarchy_loader_old():
    """
    Class for loading hierarchical yaml config files.
    """
    
    def __init__(self, path):
        """
        Constructor for Hierarchy_loader objects.
        :param: Path to a directory to load yaml files from. All *.yaml files, but not *.inc.yaml files, that reside under this directory will be loaded.
        """
        self.root_directory = path
        # This is one above the root directory, as we'll include the root_directory in the path we load from.
        self.top_directory = Path(self.root_directory, "..").resolve()
        
    def load(self):
        """
        Load the config files.
        """        
        config_paths = self.leaf_nodes(self.root_directory)
        
        processor = himl.ConfigProcessor()
        
        # Load each yaml file in turn.
        with silico.misc.directory.cd(self.top_directory):
            # Our loader.
            configs = [self.pre_process(processor.process(path = str(Path(config_path).relative_to(self.top_directory))), config_path, Path(config_path).relative_to(self.root_directory)) for config_path in config_paths]
        
        return configs
    
    def pre_process(self, config, config_path, relative_path):
        """
        Peform some pre-processing on a just-loaded config.
        
        :param config: Dictionary of config options.
        :param config_path: Path to the file from which this config was loaded.
        """
        # Do nothing.
        return config

    @classmethod
    def leaf_nodes(self, root_directory):
        """
        Return a list of bottom-most directories containined within a directoty.
        
        Bottom-most directories are those that don't contain any other directories.
        """
        nodes = [root for root, directories, files in os.walk(root_directory) if len(directories) == 0 and not root_directory.samefile(root)]
        nodes.sort()
        return nodes


class Configurable_loader_old(Hierarchy_loader_old):
    """
    Class for loading hierarchical configurables.
    """
    
    def __init__(self, *args, TYPE, **kwargs):
        """
        Constructor for Configurable_loader objects.
        
        :param TYPE: The TYPE of the configurables we are loading.
        """
        super().__init__(*args, **kwargs)
        self.TYPE = TYPE
        self.type_class = Configurable.from_class_handle(self.TYPE)
    
    def pre_process(self, config, config_path, relative_path):
        """
        Peform some pre-processing on a just-loaded config.
        
        :param config: Dictionary of config options.
        :param config_path: Path to the file from which this config was loaded.
        :param relative_path: Path to the file relative to the base directory.
        :param TYPE: The TYPE of the configurable.
        """
        config['TYPE'] = self.TYPE
        
        # Try and get the configurable class.
        try:
            cls = self.type_class.from_class_handle(config['CLASS'])
        except ValueError:
            raise Silico_exception("Error loading configurable of type '{}' from file '{}'; CLASS '{}' is not recognised".format(self.TYPE, config_path, config['CLASS'])) from None
        except KeyError:
            #raise Silico_exception("Error loading configurable of type '{}' from file '{}'; no CLASS set".format(self.TYPE, config_path)) from None
            # If no class set, use the top level class.
            cls = self.type_class
        
        configurable = cls(config_path, relative_path, False, **config) 
        configurable.configure_auto_name()
        return configurable
    
    def load(self, configurable_list = None):
        """
        Load the config files.
        
        :param configurable_list: Optional Configurable_list object to which the new configs will be added. If not given a new empty list will be used.
        """
        configurable_list = configurable_list if configurable_list is not None else Configurable_list()
        return configurable_list.update_from_list(super().load())
    
    @classmethod
    def silico_options(self, options):
        """
        Load configurables from all silico locations and add to a silico options object.
        
        :param options: A Silico_options object to populate with configurables.
        """
        # Load each of our locations in turn.
        for location in (Config_loader.MASTER_CONFIG_PATH().parent, Config_loader.SYSTEM_CONFIG_LOCATION().parent, Config_loader.USER_CONFIG_LOCATION().parent):
            # And each configurable type:
            for configurable_name, configurable_type, TYPE in (
                ("Basis Sets", "basis_sets", "basis_set"),
                ("Calculations", "calculations", "calculation"),
                ("Programs", "programs", "program"),
                ("Methods", "methods", "method")
            ):
                root_directory = Path(location, configurable_name)
                
                # Load from the location and add to our silico_options object.
                # The name of the attribute we add to is the same as the location we read from, but in lower case...
                try:
                    self(root_directory, TYPE = TYPE).load(getattr(options, configurable_type))
                except FileNotFoundError:
                    # No need to panic.
                    pass
            
        
    