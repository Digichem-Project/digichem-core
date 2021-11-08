# General imports.
import argparse
import logging
import numpy
import signal
import itertools

# Silico imports.
from silico.exception.uncatchable import Signal_caught
from silico.config.file.parser import Config_file_parser, Config_parser
import silico.result.angle
from silico.exception.base import Silico_exception


class Program():
    """
    A class for representing one of the (sub) programs of silico.
    
    Each program can be run either directly from the command line, or interactively in an urwid widget.
    """
    
    name = ""
    command = "implement in subclass"
    description = ""
    epilog = "{} V{}. Written by {}. Last updated {}.".format(name, silico.version, silico.author, silico.last_updated.strftime("%d/%m/%Y"))
    usage = None
    help = None
    
    def __init__(self, args, config, logger):
        """
        Constructor for Program objects.
        
        :param args: The command-line arguments the program was started with.
        :param config: A loaded Silico config object.
        :param logger: The logger to use for output.
        """
        self.args = args
        self.config = config
        self.logger = logger
        
    @classmethod
    def init_from_argparse(self, args):
        """
        Perform program setup from the data provided by argparse.
        
        :param args: Argparser object.
        :return: A tuple of (args, config, logger).
        """
        # First, sort out our logger.
        logger = logging.getLogger(silico.logger_name)
        
        # Next, load all our config files.
        config = self.get_silico_config(args)
        
        # Add any settings set on the command line to this merged config 'file'
        self.process_standard_arguments(args, config)
        
        # Call our arg_to_config function now so any other command line arguments can be added to the config object.
        self.arg_to_config(args, config)
        
        # Setup program wide stuff.
        self.program_init(args, config, logger)
        
        return (args, config, logger)
        
    @classmethod
    def load_from_argparse(self, args, config, logger):
        """
        Create a program instance.
        
        :param args: The command-line arguments the program was started with.
        :param config: A loaded Silico config object.
        :param logger: The logger to use for output.
        """
        return self(args = args, config = config, logger = logger)
    
    @classmethod
    def from_argparse(self, args):
        """
        Create a Program object from the data provided by argparse.
        """
        args, config, logger = self.init_from_argparse(args)
        return self.load_from_argparse(args, config, logger)
    
    @classmethod
    def get_silico_config(self, args):
        """
        Get the main silico config file to use.
        
        For most programs, we parse this from config files, but some programs do other more weird stuff (see resume).
        """
        return Config_file_parser.silico_options()
    
    @classmethod
    def subprocess_init(self, *args, **kwargs):
        """
        Init function for subprocess workers.
        This function disables a signal handler that silico normally monitors as this signal appears to be used by multiprocessing.pool for communication.
        """
        # We unset the silico signal handler for SIGTERM because multiprocessing.pool appears to use this for communication.
        # Perhaps we should unset all silico signal handlers in children?
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
    
    @classmethod
    def program_init(self, args, config, logger):
        """
        Perform program wide setup. This method should only be called once per invocation of the entire silico program.
        """
        # Set our log level.
        config.set_log_level(logger)
        
        # Set external programs.
        # TODO: Find a better way to do this.
        silico.image.vmd.VMD_image_maker.vmd_execuable = config['external']['vmd']
        silico.image.vmd.VMD_image_maker.tachyon_executable = config['external']['tachyon']
        silico.file.cube.Fchk_to_cube.cubegen_executable = config['external']['cubegen']
        silico.file.fchk.Chk_to_fchk.cubegen_executable = config['external']['formchk']
        
        # Set angles.
        silico.result.angle.Angle.set_default_angle_units(config['angle_units'])
            
        # Set numpy errors (not sure why this isn't the default...)
        numpy.seterr(invalid = 'raise', divide = 'raise')
        
        # Add our signal exception handler.
        self.init_signals(args, config, logger)
    
    @classmethod
    def init_signals(self, args, config, logger):
        """
        Setup signal handling.
        """
        # We override default signal handling by instead raising a type of exception.
        for signal_desc in ['SIGABRT', 'SIGBREAK', 'SIGBUS', 'SIGFPE', 'SIGHUP', 'SIGILL', 'SIGSEGV', 'SIGTERM', 'SIGUSR1', 'SIGUSR2']:
            signalnum = getattr(signal, signal_desc, None)
            if signalnum is not None:
                try:
                    signal.signal(signalnum, Signal_caught.raise_from_signal)
                except Exception:
                    # This signal cannot be caught, pass.
                    logger.debug("Cannot handle signal {} ({}); skipping".format(signalnum, Signal_caught.signal_to_name(signalnum)))
            else:
                logger.debug("Signal '{}' is not supported on this system; skipping".format(signal_desc))
        
    def main_wrapper(self):
        """
        Wrap the body of this program with exception catchers.
        """    
        try:
            retval = self.main()
            
            # Success.
            return retval if retval is not None else 0
        
        except KeyboardInterrupt:
            self.logger.info("interrupted by user (ctrl-c)", exc_info = self.logger.level == logging.DEBUG)
            return -1
        
        except Signal_caught:
            self.logger.error("stopped by signal", exc_info = True)
            return -2
        
        except Exception:
            # Something went wrong.
            self.logger.error("stopped with error", exc_info = True)
            
            return -2
    
#     def run(self):
#         """
#         Run this program.
#         
#         :returns: The return value of this program.
#         """
#         self.program_init()
#         
#         self.logger.debug("Startup completed")
#         
#         return self.main_wrapper()
    
    @classmethod
    def arg_to_config(self, args, config):
        """
        A class method that will be called to add command line arguments from 'args' to the configuration object 'config'.
        
        :param args: An argparser object.
        :param config: A Silico config object to add to.
        """
        # This default implementation does nothing.

    @classmethod
    def standard_arguments(self):
        """
        Get standard command line arguments which are common to all Silico sub programs.
        
        :returns: An ArgumentParser object.
        """
        standard_args = argparse.ArgumentParser(add_help=False)
        group = standard_args.add_argument_group("general options", "general options that control various aspects of silico")
        group.add_argument("-I", "--interactive", help = "Run this command interactively", action = 'store_true')
        group.add_argument("-V", "--verbose", help = "increase verbosity, stack with itself to further increase verbosity (this option overrides log_level)",  action='count', default = None)
        group.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'OFF'], help = "the level of messages to print", default = None)
        group.add_argument("-S", "--setting", help = "set a config option to a value. Options of this type are parsed as if they were a config file (in yaml format) and are then used to set corresponding options, eg -S \"absorption_graph: {fwhm: 100}\"", nargs = "*", default = [], action = "append")
        group.add_argument("--config_files", help = "an additional config file to read from. See the master config file for possible config options. Note that the master config file (at silico/data/config/silico.yaml) and user config file (at ~/.config/silico/silico.yaml) are always read automatically and do not need to be specified here. Multiple files may be given and will be processed in the order specified (the last having highest priority)", nargs = "*", default = [])
        return standard_args
    
    @classmethod
    def process_standard_arguments(self, args, config):
        """
        Add the standard command line arguments to a config object.
        
        :param args: An argparse Namespace object.
        :param config: The Silico_options object to add to.
        """    
        # Set log level and/or verbosity.
        if args.log_level is not None:
            config['logging']['log_level'] = args.log_level
            
        if args.verbose is not None:
            config['logging']['verbose'] = args.verbose
        
        # We will also process any additional config files given on the command line, and add them in order.
        for config_file_name in args.config_files:
            # Load the file and merge.
            config.merge(Config_file_parser(config_file_name).load())
            
        # Finally, we set any config options.
        # We do this last so that they'll have highest priority.
        # We read each individually because otherwise options that have been set twice won't merge properly (the yaml parser will overwrite the older).
        for config_file in itertools.chain(*args.setting):
            # Now we'll parse it as YAML.
            try:
                config.merge(Config_parser(config_file).load())
                
            except Exception:
                raise Silico_exception("failed to parse command-line config options")
            
        # All done.
    

    def main(self):
        """
        The actual logic of this program; this method should be implemented in a real subclass.
        """
        raise NotImplementedError("Implement in subclass to perform functionality")
    
    def get_interface(self, window, reload = False):
        """
        Return an urwid widget suitable for controlling this program.
        """
        if getattr(self, "_interface", None) is None or reload:
            self._interface = self.load_interface(window)
            
        return self._interface
    
    def load_interface(self, window):
        """
        Load the urwid widget suitable for controlling this program.
        """
        raise NotImplementedError("This program does not support an interface")
    
    @classmethod
    def arguments(self, sub_parsers_object):
        """
        Add this program's arguments to an argparser subparsers object.
        
        :param sub_parsers_object: An argparser subparsers object to add to (as returned by add_subparsers().
        :returns: The argparser object (which supports add_argument()).
        """
        sub_parser = sub_parsers_object.add_parser(
            self.command,
            description = self.description,
            parents = [self.standard_arguments()],
            usage = self.usage,
            epilog = self.epilog,
            help = self.help)
         
        # Set main function.
        sub_parser.set_defaults(prog_cls = self)
        
        return sub_parser
        
    def process_arguments(self, args):
        """
        Method called to alter the values set for command line arguments.
        
        This method is useful for eg setting complex default values.
        """
        # This default implementation does nothing.
    