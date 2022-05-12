.. _tutorial :

Calculation Submission Tutorial
===============================

Perhaps the most powerful aspect of Silico is its functionality for submitting and managing calculations.
This tutorial acts as a walk-through for this process and is an excellent start for performing computations with Silico.


Introduction
------------

When using Silico, there are two types of input file that control how each calculation will be performed.
These are: (1) coordinate files, which specify the elements and starting geometry of each of the molecules or systems under study,
and (2) method files, which control the specifics of each calculation, such as the functional and basis set.

One of the main advantages of using Silico for calculation submission is that multiple coordinate files can be specified at once, in which case all the given coordinate files will be submitted to the same calculation (as specified by the method file) simultaneously.
In addition, multiple method files can also be specified, in which case each coordinate file will first be submitted to the calculation defined by the first method file specified. Once each of these calculations has completed, the resulting atomic geometry will then be submitted to the calculation defined by the second method file, and so on until all method files have been exhausted.

Silico comes pre-loaded with a large database of method files, and so in most cases these do not need to be written by the user.
Hence all that is generally needed to submit a calculation is a number of coordinate file specifying the molecules of interest.

Silico supports a wide variety of coordinate formats, including both 2D and 3D formats, each of which can be written by different programs. Notable entries include:

 * Gaussian and Turbomole output (.log)
 * GaussView input (.com, .gjf, .gjc and .gau)
 * ChemDraw (.cdx and .cdxml)
 * MarvinSketch (.mrv)
 * Independant cheminformatics (.cml, .xyz etc)
 * Crystallographic information file (.cif)
 
Silico uses the OpenBabel library for file conversion. Please see `the OpenBabel documentation <https://open-babel.readthedocs.io/en/latest/FileFormats/Overview.html>`_ for a full list of supported formats.

.. note::
	Care should be taken when using 2D formats, particularly for complex 3D structures or those with specific steric information (enantiomers, for example), as the exact configuration as drawn might not be maintanined in the conversion to 3D coordinates.

In general, the closer to the final atom coordinates the input coordinates are, the fewer optimisation cycles will be required.
As such, it is strongly recommended to favour 3D formats over 2D. If available, a Crystallographic Information File (.cif) or output file from a previous calculation (on the same or similar structure) should be favoured as input, for the same reason.


Prepare Coordinate Files
------------------------

The first step in submitting a calculation is to prepare the coordinate files. These files represent the molecules upon which the calculation (or calculations) will be performed; one molecule per file.
Any program can be used to write these coordinate files (so long as the file can be saved in a format Silico understands), but generally in the EZC group the GaussView program is used.
GaussView can be run on your personal machine `via` `Apps Anywhere <https://appstore.st-andrews.ac.uk/login>`_, or on one of the group computers (where GaussView is already installed).

Upload Coordinate Files
-----------------------

Connect to the Kennedy cluster using WinSCP and upload the prepared coordinate files.
It is recommended that molecules corresponding to different projects be stored in different directories, but this is left to the discretion of the reader.

Run Silico
-------------

It is now time to run Silico, the program which will manage the remainder of the submission process.
To begin, run the ``silico submit`` subprogram followed by a list of the coordinate files to submit. For example, to submit two files named 'Benzene.cdx' and 'Naphthalene.com', run:

.. code-block:: console

	$ silico sub -I Benzene.cdx Naphthalene.com
	
.. note::
	The ``-I`` switch given in the above command informs Silico to run in interactive mode.
	
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

Finally, you may choose to not specify any input coordinates at this time, in which case they can be loaded later using the interactive interface (see :ref:`interactive coords`):

.. code-block:: console

	$ silico sub -I

	
The Interactive Interface
--------------------------

Any of the above commands will run the silico submit subprogram in interactive mode, which will open a window that appears as follows:

.. image:: /_static/submit_tutorial/submit_interface.png
    :width: 80%
    :align: center

This interface acts similarly to a graphical user interface (GUI).
The various parts of the interface can be navigated with the arrow keys.
Doing so will move the flashing cursor which indicates the part of the interface that is currently selected.


Input Coordinates
_________________

The upper section of the submission interface displays loaded input coordinates in a table format, along with the relevant molecular
formula, charge and multiplicity. These latter two columns can be edited individually for each system under study.
For example, to change the multiplicity of 'Benzene' in the above example, first move the cursor with the arrow keys to the ``mult:1`` widget for the 'Benzene' row.
Then, the old multiplicity can be removed used the backspace key, and a new multiplicity can be typed.

The three widgets in the right-most column of the coordinate table can be used to control the position of each row.
These widgets are buttons, which can be readily identified by the angle brackets surrounding the text of the button (eg, ``< Button >``).
Buttons can be 'clicked' or 'activated' by first selecting them with the arrow keys and then pressing 'enter' (or 'space').
In this case, the ``< ↑ >`` and ``< ↓ >`` buttons will move each row up or down one position respectively,
while the ``< r >`` button will delete the given row.


.. _interactive coords:

Adding New Coordinates
++++++++++++++++++++++

New input coordinate files can also be loaded at this point using the ``< Add new here >`` button, which will open the file browser:

.. image:: /_static/submit_tutorial/file_browser.png
    :width: 80%
    :align: center

This browser displays a list of files in a 'tree' like format;
each directory (or folder) appears as a 'branch' node with a '+' icon (which can be expanded to show its contents) while each file appears as a 'leaf' node (which cannot).
To expand (or 'open') a directory,  navigate up or down with the arrow keys to select it, and then use the 'right arrow' key to expand it.
An expanded directory will show a '-' icon instead of a '+'.
A directory can similarly be contracted by selecting it and using the 'left arrow' key to hide its contents.

To select a coordinate file to load, use the 'space' or 'enter' key to highlight it. If a file is chosen in error, pressing 'space' again will deselect it.
Once the files to be loaded have been selected, navigate to the ``< Confirm >`` button in the bottom right corner and select it.

..	tip::
	Instead of using the down arrow key to navigate all the way to the bottom of the page, the 'tab' key can be used to skip directly to the controls at the bottom of the window.
	Similarly, 'shift-tab' (holding shift will pressing tab) will skip back to the browser.

.. image:: /_static/submit_tutorial/file_browser_selected.png
    :width: 80%
    :align: center

This will load each of the chosen coordinate files.
Once complete, the 'Finished loading coordinates' line be printed, at which point the output window can be closed using the ``< Confirm >`` button:

.. image:: /_static/submit_tutorial/file_browser_output.png
    :width: 80%
    :align: center


Calculation Methods
___________________

The 'Calculation Methods' section of the submission interface is where the actual calculations to be performed are selected.
In most cases this will be done by selecting a (number of) methods from the built in library.
To do so, 'click' the ``< Browse library >`` button to open the method browser:

.. image:: /_static/submit_tutorial/method_browser.png
    :width: 80%
    :align: center

Conceptually, each method consists of three parts, which are:

 * The Destination: A logical or physical location where the calculation will be performed, for example a specific SLURM partition.
 * The Program: A CC progam to perform the calculation, for example Gaussian or Turbomole.
 * The Calculation: A specific set of calculation instructions, including, for example, the method, functional and basis set.

Each part of the method is chosen from the browser sequentially. This first item to choose is the destination.
On Kennedy, these destinations represent the different SLURM partitions that can be submitted to, of which only one, the `Single Node SLURM` partition, is available by default.
This single node partition should be chosen *in nearly all cases*. If you believe it does not meet your requirements, discuss with a senior computational group member about the alternative partitions.

To select the `Single Node SLURM` destination, navigate to it with the arrow keys and expand it with the right arrow key.
Doing so will reveal the computational chemistry programs that this destination supports. On Kennedy, three programs are available, which are Gaussian 09, Gaussian 16 and Turbomole:

.. image:: /_static/submit_tutorial/method_browser_program.png
    :width: 80%
    :align: center

Similarly, expanding a program will reveal the calculations that program supports. For example, the calculations the 'Gaussian 16' program supports are as follows:

.. image:: /_static/submit_tutorial/method_browser_calculation.png
    :width: 80%
    :align: center

These calculation types are grouped in a hierarchy, where the top-most item describes the general calculation type, for example an 'Optimisation' or calculation of 'Excited States'.
Within each heading the specifics of the calculation can be chosen, for example the below selection is for an optimisation using the popular B3LYP functional and 6-31G(d,p) basis set, in the gas phase:

.. image:: /_static/submit_tutorial/method_browser_selection.png
    :width: 80%
    :align: center

To choose a given method, select the final item (typically the basis set), highlight it with the 'enter' or 'space' keys and then 'click' the ``< Confirm >`` button.
It will then be added to the method table:

.. image:: /_static/submit_tutorial/method_chosen.png
    :width: 80%
    :align: center


Method Codes
++++++++++++

You will notice that each of the three items of the method is given a unique code (an integer which is greater than zero). These codes are shown both in the method browser and the method table:

.. image:: /_static/submit_tutorial/method_code_table.png
    :width: 80%
    :align: center

.. image:: /_static/submit_tutorial/method_code_browser.png
    :width: 80%
    :align: center

For example, the method chosen above has the method code of `1/2/1489`.
These method codes are unique and stable (they do not change randomly), meaning they can be used as a quick way to refer to a method.
Among other things, this allows a method to be selected by using its code alone by clicking the ``< Add from code >`` button of the method table and entering the relevant code directly:

.. image:: /_static/submit_tutorial/add_by_code.png
    :width: 80%
    :align: center

See :ref:`codes` for a table of common calculation codes.

Method Queuing
++++++++++++++

Silico allows multiple methods to be queued up to be performed one after another.
This `in-series` calculation queuing works by taking the output geometry of the previous calculation and automatically submitting it to the next calculation.
This is particularly useful for calculations that depend on a certain type of optimised geometry.
For example, the calculation of excited states typically requires a prior optimisation of the geometry which has to be performed as a separate step.
To queue up such a series of calculations, simply add a second method (or as many as are required) after the first. The methods will be processed in the same order as they appear in the table:

.. image:: /_static/submit_tutorial/method_queue.png
    :width: 80%
    :align: center

.. note::
	Methods can even be queued using different CC programs;
	the output geometry from the previous calculation will automatically be converted to an appropriate input type for the next CC progrm.
	
Submit
------

Once the desired input coordinates and calculation methods have been chosen, the selection can be submitted by selecting the ``< Confirm >`` button.
Information will be shown as each coordinate file is prepared and then submitted.
Once all files have been processed, the 'Successfully submitted x file(s)' line will appear:

.. image:: /_static/submit_tutorial/submission.png
    :width: 80%
    :align: center

Congratulations, your computations have now been submitted successfully.
You may now quit Silico (by pressing `ESC` or `ctrl-c`), or you may continue to submit further calculations.
