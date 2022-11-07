Introduction
============

About This Document
-------------------

This document is designed to act as a rapid introduction and quick reference guide for the Zysman-Colman group on how to access and use Silico; the computational chemistry management suite.

If you cannot find a solution to your problem within this document, please try referring to the full Silico manual or ask a senior group member for further help.

Images are included for comparison and reference.

Supplementary information and helpful tips are displayed separately, as follows:

.. note::
	This is a tip.
	

Commands that should be typed by the user are displayed in the following format:

.. code-block:: console

	$ echo Hello world
	
The dollar sign character ('$') itself should not be typed by the user; it is simply used to indicate that the text following should be entered as a command.
This is helpful to distinguish between commands (which should be typed) and the resulting output (which will appear without the dollar sign):

.. code-block:: console

	$ echo Hello world
	Hello world

In the above example, the user is being instructed to type: ``echo Hello world``.
The computer, in response, gives the output: ``Hello world``.

Ellipses (...) indicate that the real, full output has been truncated:

.. code-block:: console

	$ cat /etc/fstab
	# /etc/fstab: static file system information.
	#
	...


About Silico
-------------

Computational chemistry, particularly for newcomers, can be an almost impenetrably complex field.
Silico is a software package that attempts to alleviate this problem by automating and/or simplifying as many aspects of the computational chemistry process as possible.

Silico handles the entire calculation submission process (the computational chemistry pipeline) from start to finish; this means that you do not need to write batch files or worry about running computational programs (such as Gaussian or Turbomole) directly, this is all managed for you. As the user, all you need to do to start performing calculations with Silico is to (1) draw the structures you are interested in, and (2) pick the calculations from a list. Silico can even handle submitting multiple calculations in sequence, automatically. Further, Silico parses and processes the raw data from completed calculation results and presents it in formats that are easy to read, understand and distribute to other scientists. This even includes the automation of common post-analysis tasks such as the generation of molecular orbital plots, natural-transition orbitals for excited states, simulated UV-Vis, IR and emission spectra, and much more. For a full list of features and a detailed overview of the Silico project, please refer to the full Silico manual.


A Note on Terminology
---------------------

Computational chemistry programs perform calculations on collections of atoms. These atoms can theoretically describe any structure, such as a compound, polymer or molecule, but generally in the body of this document the term 'molecule' will be used to interchangeably describe the structure of interest (even though in reality, the structure may not truly be a 'molecule').


What's new in this Version (2.0)
----------------------------------

Version 2 is the second major iteration of Silico. This version includes a number of technical improvements, most notable of which is a complete rewrite of how calculation methods are stored internally. This means that all calculations now share a core base of options, regardless of the program the calculation is for. Most calculation options are now much better grouped together too, so it's easier to find the setting you are looking for. Version 2.0 also greatly expands support for Turbomole, both in terms of the number of calculations offered and also the overall handling of the rather complex Turbomole submission mechanism. In particular, support for calculating single-point gradients with ADC(2) and CC2 has been added, meaning difference density plots can now be calculated without a costly ground-state optimisation. Finally, version 2.0 introduces automated testing, so for the first time all major aspects of Silico have been rigorously tested. A full list of the changes introduced in this version are listed below:

 - Calculation options have been standardised across calculation programs.
 - The calculation method is now written to the Input directory for each calculation.
 - SMILES codes have been added to report and result files.
 - The program-independent Silico input file format (.si) has been reworked so geometry is now supported natively in YAML.
 - Geometries can now be extracted from silico.resume.pickle files.
 - Support has been added for SP gradient calculations with both Turbomole and Gaussian.
 - Support has been added for 'fast' DFT approximations (RI-DFT and RIJK-DFT) with Turbomole.
 - Support has been added for MP2, RI-MP3, RI-MP4, RI-CCSD and RI-CCSD(T) with Turbomole.
 - Support has been added for solvent simulation with Turbomole (COSMO).
 - Support has been added for numerical frequency calculations with Turbomole (including for ADC(2) and CC2).
 - Support has been added for excited state optimisations with Turbomole (including for ADC(2) and CC2).
 - An automated testing suite has been introduced.

