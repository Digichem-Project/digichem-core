from .base import Config, Auto_type

from silico.config.silico import Silico_options
from silico.config.locations import master_config_path, system_config_location, user_config_location
from silico.config.parse import Config_file_parser, Config_parser

from silico.log import get_logger
#import silico.log

# The main silico options object.
# When running as a program, this will be merged with run-time options.
options = None

def get_config(extra_config_files = None, extra_config_strings = None):
        """
        Get a Silico options object.
        
        The returned object will take options from three sources:
        1) Any config files found in the default locations (see silico.config.locations).
        2) Any additional config files specified as an argument.
        3) Any additional config strings specified as an argument.
        
        IMPORTANT: The returned options object will be merged with any additional options without prior copying.
        Hence the object returned by this function will be the same as the options attribute of this module.
        
        :return: A Silico_options object (a fancy dict).
        """
        log_level = get_logger().level
        get_logger().setLevel("DEBUG")
        
        try:
            # First, load options if not already done so.
            if options is None:
                # Load config options from file.
                # These objects are simple dicts.
                config = Config_file_parser(master_config_path).load(True)
                config.merge(Config_file_parser(system_config_location).load(True))
                config.merge(Config_file_parser(user_config_location).load(True))
                
                # No need to validate here, we're going to do it later anyway.
                globals()['options'] = Silico_options(validate_now = False, **config)
            
            if extra_config_files is None:
                extra_config_files = []
                
            if extra_config_strings is None:
                extra_config_strings = []
                
    #         if len(extra_config_files) == 0 and len(extra_config_strings) == 0:
    #             # Do nothing.
    #             return options
            
            # Load any additional config files.
            for extra_config_file in extra_config_files:
                options.deep_merge(Config_file_parser(extra_config_file).load())
                
            # Then load any additional config strings.
            for extra_config_string in extra_config_strings:
                options.deep_merge(Config_parser(extra_config_string).load())
            
            # Check everything is valid.
            options.validate()
            
            # And return.
            return options
        
        except Exception:
            raise
        
        else:
            # Only restore log level if there was no problem.
            get_logger().setLevel(log_level)
