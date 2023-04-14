"""Code for loading/reading/parsing config files and configurable files."""

import yaml

from silico.config import Config
from silico.config.silico import Silico_options
from silico.config.locations import master_config_path, system_config_location, user_config_location


class Config_parser():
    """
    Class for loading standard silico options (from a string).
    """
    
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
        try:
            config = yaml.safe_load(config_stream)
            
        except Exception as e:
            # We could ignore errors here but it might not be obvious to the user why their settings have been ignored
            raise Exception("Error parsing settings") from e
            
        # Check the config isn't None (which it will be if the file is empty.)
        try:
            if not isinstance(config, dict):
                raise Exception("Config option '{}' is formatted incorrectly".format(config_stream))
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
    def silico_options(self, extra_config_files = None, extra_config_strings = None):
        """
        Load all Silico config files.
        
        Config files are searched for in the following order (increasing precedence) 1) the silico install directory, 2) /etc/silico, 3) the user's ~/.config/silico directory
        
        Note that configurables will be added to the returned Silico_options, but they will not be resolved.
        
        :return: A Silico_options object (a fancy dict).
        """
        if extra_config_files is None:
            extra_config_files = []
            
        if extra_config_strings is None:
            extra_config_strings = []
        
        # First, we always load the 'master' (located in the silico install directory).
        config = self(master_config_path).load(True)
        
        # Next the system (located at /etc/silico)
        config.merge(self(system_config_location).load(True))
        # Finally the user specific (located at ~/.config/silico)
        config.merge(self(user_config_location).load(True))
        
        # Now load any additional config files.
        for extra_config_file in extra_config_files:
            config.merge(self(extra_config_file).load())
            
        # Then load any additional config strings.
        for extra_config_string in extra_config_strings:
            config.merge(Config_parser(extra_config_string).load())
        
        # Get a merged version.
        options = Silico_options(validate_now = True, **config)
        
        # And return.
        return options
