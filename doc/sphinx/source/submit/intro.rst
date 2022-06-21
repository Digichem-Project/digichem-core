Introduction
============

Submitting calculations is the process by which a computational chemistry program is first prepared and then executed to perform a calculation on a set of atomic coordinates.
Silico manages all parts of the calculation submission process and no direct calls to the CC program (Gaussian, Turbomole etc) are required by the user; Silico acts as the sole interface.
Silico can also interface to some job-scheduling and queuing services (currently, only SLURM is supported), in which case there is also no need for the user to to write job submission scripts manually; these will also be handled automatically.

The remainder of this section details this submission process, particularly the usage of the :ref:`silico submit <silico_submit>` subprogram. A detailed walkthrough can also be found in the :ref:`submission tutorial <submit_tutorial>` section, which is an ideal starting point for newcomers to Silico.