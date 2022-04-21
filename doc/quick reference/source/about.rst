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
A non-exhaustive list of features is as follows:

	* Submission to computational programs through a simple and unified interface.
	* Simultaneous submission of multiple molecules/systems.
	* Automatic in series submission of results from completed calculations to subsequent calculations.
	* Automatic conversion of input files (including ChemDraw) to types appropriate for the selected computational program.
	* Automatic and manual analysis of computation results, including tabulation to CSV format.
	* Automatic and manual generation of PDF reports from computation results, including rendered 3D structures, orbital images and graphs.
	* Various post-analyses including spin-orbit coupling (SOC), natural-transition orbitals (NTOs) and differential-density plots.

