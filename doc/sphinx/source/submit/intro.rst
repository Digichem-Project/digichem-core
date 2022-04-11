Submitting Calculations
=======================

Submitting calculations is the process by which a computational chemistry program is first prepared and then executed to perform a calculation on a set of atomic coordinates.
Silico manages all parts of the calculation submission process and no direct calls to the CC program (Gaussian, Turbomole etc) are required by the user; Silico acts as the sole interface.
Silico can also interface to some job-scheduling and queuing services (currently, only SLURM is supported), in which case there is also no need for the user to to write job submission scripts manually; these will also be handled automatically.

Calculation submission is handled by the ``silico submit`` sub-progam, which is detailed in this section, and a detailed, step-by-step :ref:`walk-through <tutorial>` is also available for those new to Silico.

.. toctree::
    :caption: Contents:
    :maxdepth: 2
	
    silico-submit
    tutorial
    input-files
    structure
    monitoring
    post