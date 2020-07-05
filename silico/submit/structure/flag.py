from enum import Enum, auto

class Flag(Enum):
	"""
	Class that represents a flag; an empty text file whose name conveys state about a calculation.
	"""
	# Calculation has been formally submitted but no work has been done yet (used by SLURM when we're in the queue).
	PENDING =  auto()
	# Calculation has started. Unlike RUNNING, this flag is permanent.
	STARTED =  auto()
	# Calculation is ongoing. Unlike STARTED, this flag is deleted once the calculation has finished (successfully or otherwise).
	RUNNING =  auto()
	# Calculation has finished successfully.
	SUCCESS = auto()
	# Performing cleanup.
	CLEANUP = auto()
	# Calculation has finished with error.
	ERROR = auto()
	# Performing post analysis on the calculation.
	POST = auto()
	# Optimisation has converged.
	CONVERGED = auto()
	# Optimisation has not converged.
	NOT_CONVERGED = auto()
	# All work has finished (including post analysis, if appropriate).
	DONE = auto()
	