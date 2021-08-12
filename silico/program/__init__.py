import logging
import numpy
import sys
import signal
import itertools
import argparse
import deepmerge

# Hidden import.
import silico.reference
from silico.exception.uncatchable import Signal_caught
from silico.exception import Silico_exception
import silico.image.vmd
import silico.file.cube
import silico.file.fchk
import silico.result.angle
import silico.config
import silico.logging
from silico.config.file.parser import Config_file_parser, Config_parser

def main_wrapper(inner_func, *, args, arg_to_config = None, logger_name = None, **kwargs):
    """
    Wrapper function for program's main function.
    
    This function calls inner_func (passing args, config and logger), watching for errors and logging appropriately to logger.
    
    :returns: 0 on success, -1 on catching KeyboardInterrupt, -2 on other error.
    """
    logger_name = logger_name if logger_name is not None else silico.logger_name
    # First, sort out our logger.
    logger = logging.getLogger(silico.logger_name)
    
    # Use our generic init function.
    try:
        args, config, logger = init_program(
            args = args,
            arg_to_config = arg_to_config,
            logger = logger)
    except Exception:
        # Couldn't start.
        logger.error("failed to start program", exc_info = True)
        return -2
    
    return run(inner_func, args = args, config = config, logger = logger, **kwargs)
    
def run(inner_func, **kwargs):
    """
    Wrapper function for program's main function.
    
    This function calls and returns the return value of inner_func, **kwargs is passed directly to inner_func.
    """
    logger = logging.getLogger(silico.logger_name)
    try:
        logger.debug("Startup completed")
        retval = inner_func(**kwargs)
        
        # Success.
        return retval if retval is not None else 0
    except KeyboardInterrupt:
        logger.info("interrupted by user (ctrl-c)", exc_info = logger.level == logging.DEBUG)
        return -1
    except Signal_caught:
        logger.error("stopped by signal", exc_info = True)
        return -2
    except Exception:
        # Something went wrong.
        logger.error("stopped with error", exc_info = True)
        return -2

def init_program(*, args, arg_to_config = None, logger):
    """
    Common initialisation routines for all programs in the silico package.
    
    :param args: An argparse namespace object.
    :param logger: The logger to innit. This logger will print to STDERR and will automatically have its log-level set based on the command-line arguments and config file.
    :param arg_to_config: An optional function (or something like a function) that will be called: arg_to_config(args, config) where 'args' will be the argparse namespace object and 'config' will be the silico config object. The intention is this function should set relevant values in the config object from the argparse arguments, if so desired.
    :return: A tuple of the argparse namespace object, silico config object and logging logger object as (args, config, logger). The logger object can of course also be obtained via logging.getLogger(logger_name).
    """
    # Next, load all our config files.
    config = Config_file_parser.silico_options()
    
    # Add any settings set on the command line to this merged config 'file'
    process_standard_args(config, args)
    
    # If we were given an arg_to_config function call it now so any other command line arguments can be added to the config object.
    if arg_to_config is not None:
        arg_to_config(args, config)
    
    # Set what settings we can.
    init_from_config(logger, config)
    
    # Set numpy errors (not sure why this isn't the default...)
    numpy.seterr(invalid = 'raise', divide = 'raise')
    
    # Add our signal exception handler.
    init_signals(logger)
        
    # Give back our command line arguments and config options.
    return (args, config, logger)

def init_signals(logger):
    """
    Setup signal handling.
    """
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


def init_from_config(logger, config):
    """
    Init program wide stuff from a loaded config object.
    
    Silico still has a few program-wide global variables (get off my back; it's in beta) which we set here.
    These globals are:
        - logging object
        - paths to external programs (VMD, cubegen etc).
        - angle units 
    
    Once these globals are refactored away, this function will disappear.
    """
    # Set what settings we can.
    # Set our log level.
    config.set_log_level(logger)
    
    # Set external programs.
    silico.image.vmd.VMD_image_maker.vmd_execuable = config['external']['vmd']
    silico.image.vmd.VMD_image_maker.tachyon_executable = config['external']['tachyon']
    silico.file.cube.Fchk_to_cube.cubegen_executable = config['external']['cubegen']
    silico.file.fchk.Chk_to_fchk.cubegen_executable = config['external']['formchk']
    
    # Set angles.
    silico.result.angle.Angle.set_default_angle_units(config['angle_units'])

standard_args = argparse.ArgumentParser(add_help=False)
group = standard_args.add_argument_group("general options", "general options that control various aspects of silico")
group.add_argument("-V", "--verbose", help = "increase verbosity, stack with itself to further increase verbosity (this option overrides log_level)",  action='count', default = None)
group.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'OFF'], help = "the level of messages to print", default = None)
group.add_argument("-S", "--setting", help = "set a config option to a value. Options of this type are parsed as if they were a config file (in yaml format) and are then used to set corresponding options, eg -S \"absorption_graph: {fwhm: 100}\"", nargs = "*", default = [], action = "append")
group.add_argument("--config_files", help = "an additional config file to read from. See the master config file for possible config options. Note that the master config file (at silico/data/config/silico.yaml) and user config file (at ~/.config/silico/silico.yaml) are always read automatically and do not need to be specified here. Multiple files may be given and will be processed in the order specified (the last having highest priority)", nargs = "*", default = [])

def process_standard_args(config, args):
    """
    Add the standard command line arguments to a config object.
    
    :param config: The Silico_options object to add to.
    :param args: An argparse Namespace object.
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


def subprocess_init(*args, **kwargs):
    """
    Init function for subprocess workers.
    This functin disables a signal handler that silico normally monitors as this signal appears to be used by multiprocessing.pool for communication.
    """
    # We unset the silico signal handler for SIGTERM because multiprocessing.pool appears to use this for communication.
    # Perhaps we should unset all silico signal handlers in children?
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
