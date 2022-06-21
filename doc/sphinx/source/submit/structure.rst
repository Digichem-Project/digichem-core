.. _submit folder structure:

Calculation Folder Structure
============================

Silico will create a series of folders for each coordinate file and each method submitted.
In each case, the topmost folder that is created will be named after the coordinate file, one directory per file.
Within this `molecule directory`, a separate folder will then be created for each calculation performed.
The name of this `calculation directory` will be taken from the method that was used to submit it.
If the same coordinate file is submitted to the same method more than once, a number will be appended to the calculation directory to ensure each calculation is performed in a unique directory.

.. note::
	The name of the ‘molecule’ directory is simply taken from the name of each coordinate file; Silico does not try to automatically name structures.
	
Within each calculation directory, a number of sub folders will be created. They are as follows:

Input
-----

Contains the input file(s) for the calculation, in the format appropriate for the chosen CC program.
The file(s) will be fully prepared for the calculation program, so they can be inspected to determine the specific parameters of the calculation.
If a non-native coordinate format was originally used in the submission (ChemDraw, for example), then this will have been converted appropriately.

.. note::
	The file originally submitted is not stored in the Input directory; only the prepared input file is saved.
	
Output
------

Contains output files written by the calculation program, and only those files (or derivatives thereof).

.. note::
	For Guassian calculations, Silico will automatically convert .chk files to .fchk (before deleting the original .chk) in order to save file space.
	Occasionally, the original .chk file is required for post-analysis, in which case this option can be disabled.
	
Logs
----

Contains log files written by Silico (in plain text).
	
Flags
-----

Contains file flags, text files that convey information about the status of the calculation.
See the section on :ref:`file_flags` for more information about file flags.

Results
-------

Contains text result files that are automatically written by Silico during :ref:`post analysis <result_files>`. This folder will only be created once the calculation has been completed.

Report
------

Contains a PDF report summary of the completed calculation, along with the rendered images used in the report. 
This folder is likewise only created once the calculation has been completed.
