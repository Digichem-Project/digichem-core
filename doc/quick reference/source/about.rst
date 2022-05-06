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
Silico is a software package that attempts to alleviate this problem by automating and/or simplifying
as many aspects of computational chemistry as possible, so that the computational process
becomes as close to a black-box as can be achieved.

Silico handles the entire calculation submission process from start to finish; this means that you do not need to write batch files or worry about running computational programs (such as Gaussian or Turbomole) directly, Silico manages all this for you. As the user, all you need to do to start performing calculations is to draw the structures you are interested in and then pick a calculation from a list, that is all. Silico can even handle submitting multiple calculations in sequence, automatically. Further, Silico parses and processes the raw data from completed calculation results and presents it in formats that are easy to read, understand and distribute to other scientists. This even includes the automation of common post-analysis tasks such as the generation of molecular orbital plots, natural-transition orbitals for excited states, simulated UV-Vis, IR and emission spectra, and much more. For a full list of features and a detailed overview of the Silico project, please refer to the full Silico manual.


A Note on Terminology
---------------------

Computational chemistry programs perform calculations on collections of atoms.
These collections of atoms can theoretically describe any structure, ranging in scale and complexity from a single element to a molecule, an ionic compound, a polymer, organometallic complex, protein or even a crystal lattice.
This large scope can make referring to these collections of atoms in an inclusive way difficult, and often the term 'system of interest' is used so as to not exclude any of the above group.
However, this term is somewhat clumsy, and so instead the term 'molecule' will be generally used in the discussion of this document for simplicity's sake.
As such, the term 'molecule' can be exchanged for any of the above terms wherever it is encountered, unless otherwise noted.
