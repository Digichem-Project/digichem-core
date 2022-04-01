Calculation Submission Tutorial
===============================

Perhaps the most powerful aspect of Silico is its functionality for submitting and managing calculations.
This tutorial acts as a walk-through for running through this process and is an excellent start for performing computations with Silico.


Prepare Files For Submission
-------------------------------

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
	
The first step in the submission process is to write the systems that are to be studied to respective coordinate files. The choice of file format is essentially irrelevant, and different file formats can be freely intermixed. Once written, these coordinate files should be transfered to the machine on which the calculation is to be performed (most commonly a remote server or cluster) by whatever method is most convenient (File-Transfer Protocol (FTP), SSH File Transfer Protocol (SFTP) etc).


Method Files
________________

The details of the calculation to perform (method, functional, basis set etc) are specified by a number of method files, one calculation per file.
However, Silico comes pre-loaded with a large database of method files which are likely to satisfy the typical user, in which case no method files need to be written at this stage. For more advanced usage, see :ref:`Writing Method Files`\ .


Connect To The Calculation Server
------------------------------------

In most cases, computational chemistry programs are run on large, distributed server clusters.
If this setup applies to you, at this point you should connect to the server cluster where Silico is installed (for example, by SSH, PuTTy or equivalent).


Run Silico
-------------

Once the required coordinate files (and optionally method files) have been prepared, it is time to run Silico.
This can be done in either a non-interactive or interactive manner (see :ref:`Running Interactively`\ ); for the purposes of this tutorial the interactive interface will primarily be used.
To begin, run the silico submit sub-program followed by a list of the coordinate files to submit. For example, to submit two files named 'Benzene.cdx' and 'Naphthalene.com', run:

.. code-block:: console

	silico sub -I Benzene.cdx Naphthalene.com
	
At this point, an explicit charge and/or multiplicity can be set using the ``-C`` (or ``--charge``) and ``-M`` (or ``--multiplicity``) options. If given, these options will overwrite any charge or multiplicity given in the coordinate files, for all specified coordinate files. For example, to submit all calculations as a radical cation:

.. code-block:: console

	silico sub -I Benzene.cdx Naphthalene.com -C 1 -M 2
	
If any of the file names contain whitespace, or other 'unusual' characters, remember to use quotation marks:

.. code-block:: console

	silico sub -I "Benz ene.cdx" Naphthalene.com

If any of the coordinate files are not in the current directory, the full path should be specified (including directories):

.. code-block:: console

	silico sub -I Aromatic/Benzene.cdx Aromatic/Naphthalene.com
	
.. note::
	Alternatively, you can change the current directory using the ``cd`` command, for example ``cd Aromatic``.
	
.. note::
	You can check which files are in the current directory using the ``ls`` command.

In additional to individual coordinate files, the contents of entire directories can be submitted `via` the unix wildcard character (*):

.. code-block:: console

	silico sub -I Aromatic/*

	
The Interactive Interface
----------------------------

Any of the above commands will run the silico submit sub-program in interactive mode, which will open a window that appears as follows:

.. image:: /_static/submit_interface.png

This interface acts similarly to a graphical user interface (GUI) (although technically it is not).
The various parts of the interface can be navigated by the arrow keys.
Doing so will move the flashing cursor which indicates the part of the interface which is currently selected.


Input Coordinates
_____________________

The upper section of the submission interface displays loaded input coordinates in a table format, along with the relevant molecular
formula, charge and multiplicity. These latter two columns can be edited individually for each system under study.
For example, to change the multiplicity of 'Benzene' in the above example, first move the cursor with the arrow keys to the ``mult:1`` widget for the 'Benzene' row.
Then, the old multiplicity can be removed used the backspace key, and a new multiplicity can be typed.

New input coordinate files can also be loaded at this point using the 'Add new here' widget.
This widget is a button, which can be readily identified by the angle brackets surrounding the text of the button (eg, ``< Button >``).
Buttons can be 'clicked' or 'activated' by first selecting them with the cursor and then pressing 'enter' (or 'space').


The File Browser
++++++++++++++++++++++

Selecting the '< Add new here >' button will open the file browser, which will appear as follows:

.. image:: /_static/file_browser.png

This browser displays a list of files in a 'tree' like format;
each directory (or folder) appears as a 'branch' node (which can be expanded to show its contents) while each file appears as a 'leaf' node (which cannot).
To expand (or 'open') a directory,  navigate up or down with the arrow keys to select it, and then use the 'right arrow' key to expand it.
A directory can similarly be contracted by selecting it and using the 'left arrow' key to hide its contents.






