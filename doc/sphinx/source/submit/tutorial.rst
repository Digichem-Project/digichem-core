.. _tutorial :

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

	$ silico sub -I Benzene.cdx Naphthalene.com
	
At this point, an explicit charge and/or multiplicity can be set using the ``-C`` (or ``--charge``) and ``-M`` (or ``--multiplicity``) options. If given, these options will overwrite any charge or multiplicity given in the coordinate files, for all specified coordinate files. For example, to submit all calculations as a radical cation:

.. code-block:: console

	$ silico sub -I Benzene.cdx Naphthalene.com -C 1 -M 2
	
If any of the file names contain whitespace, or other 'unusual' characters, remember to use quotation marks:

.. code-block:: console

	$ silico sub -I "Benz ene.cdx" Naphthalene.com

If any of the coordinate files are not in the current directory, the full path should be specified (including directories):

.. code-block:: console

	$ silico sub -I Aromatic/Benzene.cdx Aromatic/Naphthalene.com
	
.. note::
	Alternatively, you can change the current directory using the ``cd`` command, for example ``cd Aromatic``.
	
.. note::
	You can check which files are in the current directory using the ``ls`` command.

In additional to individual coordinate files, the contents of entire directories can be submitted `via` the unix wildcard character (*):

.. code-block:: console

	$ silico sub -I Aromatic/*

Finally, you may choose to not specify any input coordinates at this time, in which case they can be loaded later using the interactive interface (see :ref:`interactive coords`)

.. code-block:: console

	$ silico sub -I

	
The Interactive Interface
----------------------------

Any of the above commands will run the silico submit sub-program in interactive mode, which will open a window that appears as follows:

.. image:: /_static/submit_tutorial/submit_interface.png

This interface acts similarly to a graphical user interface (GUI).
The various parts of the interface can be navigated by the arrow keys.
Doing so will move the flashing cursor which indicates the part of the interface which is currently selected.


Input Coordinates
_________________

The upper section of the submission interface displays loaded input coordinates in a table format, along with the relevant molecular
formula, charge and multiplicity. These latter two columns can be edited individually for each system under study.
For example, to change the multiplicity of 'Benzene' in the above example, first move the cursor with the arrow keys to the ``mult:1`` widget for the 'Benzene' row.
Then, the old multiplicity can be removed used the backspace key, and a new multiplicity can be typed.

The three widgets in the right-most column of the coordinate table can be used to control the position of each row.
These widgets are buttons, which can be readily identified by the angle brackets surrounding the text of the button (eg, ``< Button >``).
Buttons can be 'clicked' or 'activated' by first selecting them with the cursor and then pressing 'enter' (or 'space').
In this case, the ``< ↑ >`` and ``< ↓ >`` buttons will move each row up or down one position respectively,
while the ``< r >`` button will delete the given row.


.. _interactive coords:

Adding New Coordinates
++++++++++++++++++++++

New input coordinate files can also be loaded at this point using the ``< Add new here >`` button, which will open the file browser:

.. image:: /_static/submit_tutorial/file_browser.png

This browser displays a list of files in a 'tree' like format;
each directory (or folder) appears as a 'branch' node with a '+' icon (which can be expanded to show its contents) while each file appears as a 'leaf' node (which cannot).
To expand (or 'open') a directory,  navigate up or down with the arrow keys to select it, and then use the 'right arrow' key to expand it.
An expanded directory will show a '-' icon instead of a '+'.
A directory can similarly be contracted by selecting it and using the 'left arrow' key to hide its contents.

To select a coordinate file to load, use the 'space' or 'enter' key to highlight it. If a file is chosen in error, pressing 'space' again will deselect it.
Once the files to be loaded have been selected, navigate to the ``< Confirm >`` button in the bottom right corner and select it.

..	note::
	Instead of using the down arrow key to navigate all the way to the bottom of the page, the 'tab' key can be used to skip directly to the controls at the bottom of the window.
	Similarly, 'shift-tab' (holding shift will pressing tab) will skip back to the browser.

.. image:: /_static/submit_tutorial/file_browser_selected.png

This will load each of the chosen coordinate files.
Once complete, the 'Finished loading coordinates' line be printed, at which point the output window can be closed using the ``< Confirm >`` button:

.. image:: /_static/submit_tutorial/file_browser_output.png


Calculation Methods
___________________

The 'Calculation Methods' section of the submission interface is where the actual calculations to be performed are selected.
In most cases this will be done by selecting a (number of) methods from the built in library.
To do so, 'click' the ``< Browse library >`` button to open the method browser:

.. image:: /_static/submit_tutorial/method_browser.png

Conceptually, each method consists of three parts, which are:

 * The Destination: A logical or physical location where the calculation will be performed, for example a specific SLURM partition.
 * The Program: A CC progam to perform the calculation, for example Gaussian or Turbomole.
 * The Calculation: A specific set of calculation instructions, including, for example, the method, functional and basis set.

Each part of the method is chosen from the browser sequentially. This first item to choose is the destination. These will have been setup by your server administrator, so you should contact them if you are unsure which to choose.
In the above example there is only one destination available, a SLURM partition called 'Single Node SLURM'.
This destination can be expanded in the same manner as the file browser above, by navigating with the arrow keys and expanding the node with the right arrow key.

Doing so will reveal the CC programs that this destination supports, in this example there is a Gaussian 09, a Gaussian 16 and a Turbomole program available:

.. image:: /_static/submit_tutorial/method_browser_program.png

Similarly, expanding a program will reveal the calculations that program supports. For example, the calculations the 'Gaussian 16' program supports are as follows:

.. image:: /_static/submit_tutorial/method_browser_calculation.png

These calculation types are grouped in a hierarchy, where the top-most item describes the general calculation type, for example an 'Optimisation' or calculation of 'Excited States'.
Within each heading the specifics of the calculation can be chosen, for example the below selection is for an optimisation using the popular B3LYP functional and 6-31G(d,p) basis set, in the gas phase:

.. image:: /_static/submit_tutorial/method_browser_selection.png

To choose a given method, select the final item (typically the basis set), highlight it with the 'enter' or 'space' keys and then 'click' the ``< Confirm >`` button.
It will then be added to the method table:

.. image:: /_static/submit_tutorial/method_chosen.png


Method Codes
++++++++++++

You will notice that each of the three items of the method is given a unique code (an integer which is greater than zero). These codes are shown both in the method browser and the method table:

.. image:: /_static/submit_tutorial/method_code_table.png

.. image:: /_static/submit_tutorial/method_code_browser.png

For example, the method chosen above has the method code of `1/2/1489`.
These method codes are unique and stable (they do not change randomly), meaning they can be used as a quick way to refer to a method.
Among other things, this allows a method to be selected by using its code alone by clicking the ``< Add from code >`` button of the method table and entering the relevant code directly:

.. image:: /_static/submit_tutorial/add_by_code.png


Method Queuing
++++++++++++++

Silico allows multiple methods to be queued up to be performed one after another.
This `in-series` calculation queuing works by taking the output geometry of the previous calculation and automatically submitting it to the next calculation.
This is particularly useful for calculations that depend on a certain type of optimised geometry.
For example, the calculation of excited states typically requires a prior optimisation of the geometry which has to be performed as a separate step.
To queue up such a series of calculations, simply add a second method (or as many as are required) after the first. The methods will be processed in the same order as they appear in the table:

.. image:: /_static/submit_tutorial/method_queue.png

.. note::
	Methods can even be queued using different CC programs;
	the output geometry from the previous calculation will automatically be converted to an appropriate input type for the next CC progrm.
	
Submit
------

Once the desired input coordinates and calculation methods have been chosen, the selection can be submitted by selecting the ``< Confirm >`` button.
Information will be shown as each coordinate file is prepared and then submitted.
Once all files have been processed, the 'Successfully submitted x file(s)' line will appear:

.. image:: /_static/submit_tutorial/submission.png

Congratulations, the input coordinates have now been submitted and Silico can be exited (by using the 'esc' key), or further calculations can be queued.
