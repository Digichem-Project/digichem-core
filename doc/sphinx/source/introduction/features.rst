.. _features:

Features
========

A general list of the features that Silico offers for computational chemists is described below:

 * Submission to computational programs through a simple and unified interface.
 * Simultaneous submission of multiple molecules/systems.
 * Automatic in series submission of results from completed calculations to subsequent calculations.
 * Automatic conversion of input files (including ChemDraw) to formats appropriate for the selected computational program.
 * Automatic and manual analysis of computation results, including tabulation to comma-separated values (CSV) format.
 * Automatic and manual generation of PDF reports from computation results, including rendered 3D structures, orbital images and graphs.


Separation Of Input Coordinates and Calculation Method
------------------------------------------------------

Traditionally, each computation chemistry (CC) program used to perform calculations (Gaussian, Turbomole, Orca etc.) is associated with a particular (and often unique) input file type.
Generally, these files are made up of two sections which dictate how the calculation is performed:

 * The Atomic Coordinates, specifying the position and element of the system (molecule) under study.
 * The Calculation Commands, specifying the nature of the calculation to perform, such as the method, functional, basis set and associated options.
 
In addition to generally being unique to each CC program, this design pattern can make submission of large numbers of calculations difficult.
This is both because it is tedious and time-consuming to duplicate the calculation commands across multiple files, but also because accidental differences can be introduced across the series of files which might not be detected until after the calculation(s) have completed, which may in extreme cases take days or more.

Silico aims to alleviate these traditional problems by separating the atomic coordinates and calculation commands into separate coordinate and method files respectively,
and by providing a standard interface for both by which all CC programs can be accessed. Additionally, Silico supports the submission of multiple (near infinite, limited only by the resources of the operating system) coordinate files to a single method simultaneously, facilitating large-scale computational screens with ease.


Coordinate Files
____________________

The molecules or systems to study are specified by coordinate files, one molecule per file. Silico supports a wide range of coordinate formats, mostly by interfacing to the obabel file conversion library.\ :cite:p:`Openbabel` This includes, but is not limited to:

 * Specific CC program input files, including Gaussian input files (.gjf, .com etc) and Turbomole coordinate files (.turb). These proprietary formats may be used regardless of the CC program actually being used to perform the calculation, for example a .com file may be used to submit a Turbomole calculation.
 * Independant 3D formats, including .xyz and .cml.
 * 2D drawing formats, including .cdx, .cdxml (ChemDraw) and .mrv (MarvinSketch). These will be automatically converted to an appropriate 3D representation.
 * Calculation output files, including .log files from Gaussian and Turbomole.
 * Crystallographic formats, most notable the .cif format.
 * An independant Silico format, the .si format.
 
.. note::
	Care should be taken when using 2D formats, particularly for complex 3D structures or those with specific steric information (enantiomers, for example).
	The conversion from 2D to 3D employs a rapid molecular-mechanics (MM) optimisation provided by the obabel library\ :cite:p:`Openbabel`. In many cases this will result in a satisfactory starting structure for further optimisation, but occassionally the geometry will become locked in an impossible or high-energy conformation. Similarly, steric information may be destroyed by the optimisation process. In these cases it is recommended to first convert the 2D coordinates to a 3D representation using the convert sub-program and manually inspect the resulting geometry prior to submission.
	
.. note::
	When using a coordinate file that also includes calculation commands (for example, the Gaussian input format), these commands will be ignored. However, charge and multiplicity information, if present, will be respected.


Method Files
________________

The details of the calculation to perform (method, functional, basis set etc) are specified by a number of method files, one calculation per file.
However, Silico comes pre-loaded with a large database of method files which are likely to satisfy the typical user, in which case no method files need to be written at this stage. For more advanced usage, see :ref:`Writing Method Files`\ .