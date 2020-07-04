import silico.image.vmd
import silico.file.cube
import silico.file.fchk
import silico.result.angle
import silico.config
import silico.logging
import logging
import numpy
import sys
# Hidden import.
import silico.references
from silico.exception.uncatchable import Signal_caught
import signal

#def main_wrapper(logger, inner_func, **kwargs):
def main_wrapper(inner_func, *, arg_parser, arg_to_config = None, logger_name = None, **kwargs):
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
			arg_parser = arg_parser,
			arg_to_config = arg_to_config,
			logger = logger)
	except Exception:
		# Couldn't start.
		logger.error("failed to start program", exc_info = True)
		return -2
		
	try:
		retval = inner_func(args, config, logger, **kwargs)
		
		# Success.
		return retval if retval is not None else 0
	except KeyboardInterrupt:
		logger.error("interrupted by user (^c)", exc_info = True)
		return -1
	except Signal_caught:
		logger.error("stopped by signal", exc_info = True)
		return -2
	except Exception:
		# Something went wrong.
		logger.error("stopped with error", exc_info = True)
		return -2

def init_program(*, arg_parser, arg_to_config = None, logger):
	"""
	Common initialisation routines for all programs in the silico package.
	
	:param arg_parser: An argparse.ArgumentParser object that has already the desired program arguments set up. General, silico-wide options will be added by this function.
	:param logger: The logger to innit. This logger will print to STDERR and will automatically have its log-level set based on the command-line arguments and config file.
	:param arg_to_config: An optional function (or something like a function) that will be called: arg_to_config(args, config) where 'args' will be the argparse namespace object and 'config' will be the silico config object. The intention is this function should set relevant values in the config object from the argparse arguments, if so desired. Optionally, a string can be given, in which case all command-line arguments will be added (without additional modification) to the config object under the key name given by the string. If arg_to_config is None, nothing will be done.
	:return: A tuple of the argparse namespace object, silico config object and logging logger object as (args, config, logger). The logger object can of course also be obtained via logging.getLogger(logger_name).
	"""
	arg_to_config = arg_to_config if arg_to_config is not None else "command_line_args"
	
	# Deal with command line arguments and configuration options.
	# First, add the general silico command line arguments, which are common to all programs (things like verbosity etc) to our command line parser object.
	silico.config.Config.get_standard_args(arg_parser)
	# Process command line arguments.
	args = arg_parser.parse_args()
	
	# Now get a configuration options object (which is a fancy dictionary).
	# This constructor automatically loads in options from the two config files.
	config = silico.config.Config.standard_configs()

	# Add general-silico command line arguments to our config object.
	config.add_standard_args(args)
	
	# If we were given an arg_to_config function call it now so any other command line arguments can be added to the config object.
	if isinstance(arg_to_config, str):
		# Just add all our command line arguments *as-is* to the config object under a new key.
		config.add_config({arg_to_config: vars(args)})
	elif arg_to_config is not None:
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

# def init_logger(logger_name):
# 	"""
# 	Init the program wide logger.
# 	"""
# 	logger = logging.getLogger(logger_name)
# 	
# 	# The console handler, where we'll print most messages.
# 	consoleHandler = logging.StreamHandler(sys.stderr)
# 	# Handle everything.
# 	consoleHandler.setLevel(logging.DEBUG)
# 	# Set its formatter.
# 	consoleHandler.setFormatter(silico.logging.Variable_formatter(logger))
# 	# Add the handler.
# 	logger.addHandler(consoleHandler)
# 	
# 	return logger

def init_from_config(logger, config):
	"""
	Init program wide stuff from a loaded config object.
	
	Silico still has a few program-wide global variables (get off my back; it's in beta) which we set here.
	These globals are:
		- logging object
		- paths to external programs (VMD, cubegen etc).
		- angle units 
	
	Once these globals are refactored away, this function will dissapear.
	"""
	# Set what settings we can.
	# Set our log level.
	config.set_log_level(logger)
	
	# Set external programs.
	silico.image.vmd.VMD_image_maker.vmd_execuable = config['external']['vmd']
	silico.image.vmd.VMD_image_maker.tachyon_executable = config['external']['tachyon']
	silico.file.cube.Cube_maker.cubegen_executable = config['external']['cubegen']
	silico.file.fchk.Fchk_maker.cubegen_executable = config['external']['formchk']
	
	# Set angles.
	silico.result.angle.Angle.set_default_angle_units(config['angle_units'])



