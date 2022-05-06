.. _file_flags:

Flags
=====

Within each calculation directory, Silico manages a special folder with the name ‘Flags’.
This folder contains a number of empty text files where the name of each file conveys status about the calculation.
These files are created and destroyed at key points in the calculation submission process, so the calculation can be monitored by observing which files are present (or absent) from
the ‘Flags’ folder at any given moment.

Available Flags
---------------

The currently available file flags are as follows:

.. list-table::
    :widths: 20 80
    :header-rows: 1

    * - Name
      - Description
    * - PENDING
      - The calculation has been submitted but has not yet begun; most likely because it is waiting in the queue.
    * - STARTED
      - The calculation has begun. This flag is never deleted, so it is useful for confirming that the calculation at least started, even if it did not finish.
    * - RUNNING
      - The calculation is currently ongoing.
    * - SUCCESS
      - The calculation finished successfully.
    * - CONVERGED
      - The optimisation converged successfully. This flag is only used for optimisation calculations.
    * - NOT_CONVERGED
      - The optimisation did not converge successfully. This flag is only used for optimisation calculations.
    * - CLEANUP
      - The calculation has finished (successfully or otherwise) and Silico is currently cleaning up (saving files etc).
    * - ERROR
      - The calculation has stopped because an error occurred.
    * - POST
      - The calculation has finished and Silico is currently performing post analysis (writing result and report files).
    * - DONE
      - All work (including post-analysis) has been completed; Silico will not make any changes after this flag. It is safe to move, download or delete the calculation folder.
    