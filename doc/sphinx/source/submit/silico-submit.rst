.. _silico_submit:

Silico Submit
=============================

Calculations are submitted with Silico using the ``silico submit`` command.
To perform a calculation, generally two things need to be provided by the user: (1) a (number of) input coordinate files and (2), a (number) of methods.

Usage
-----

Coordinate files can be given anywhere after the ``silico submit`` command, for example:

.. code-block:: console

	$ silico sub coordinate.file
	
Multiple coordinate files can be specified too; they will all be submitted sequentially to the same method.

.. code-block:: console

	$ silico sub coordinate1.file coordinate2.file coordinate3.file ...
	
.. note::
	Silico does not specify a maximum for the number of coordinate files that can be submitted at any one time.
	However, the operating system and shell will typically impose a maximum character length for a single command and this will act as an upper limit for the number of files that can be submitted simultaneously.
	Fortunately, this limit is typically very large in modern systems and so is unlikely to pose a problem in most use-cases.
	In instances where this limit is reached, it can be overcome by splitting the input coordinate files into separate groups and submitting each separately.
	
Methods, which define the calculations to be performed, can be specified in a number of ways.
Method `files` can be specified using the ``-m`` (or ``--method-files``) option:

.. code-block:: console

	$ silico sub coordinate.file -m method.file
	
While methods from the built in library can be specified using the ``-c`` (or ``--method-codes``) option and the unique code:

.. code-block:: console

	$ silico sub coordinate.file -c 1/1/1
	
Or even by unique name, if these are known:

.. code-block:: console

	$ silico sub coordinate.file -c "Single Node SLURM/Gaussian 16/[Gaussian Optimisation, Gaussian B3LYP (GD3BJ), Gaussian Gas Phase, 'Gaussian 6-31G(d,p)']"

Multiple methods can also be specified, in which case each coordinate file will be automatically submitted to each specified method in sequence (using the output goemetry from the previous method in each case):

.. code-block:: console

	$ silico sub coordinate.file -c 1/1/1 1/1/2 ...
	
.. note::
	Similarly to coordinate files, Silico does not impose a limit on the number of methods that can be queued in this way.
	
Method files and method codes can also be freely intermixed, as desired:

.. code-block:: console

	$ silico sub coordinate.file -c 1/1/1 -m method2.file -c 1/1/2 -m method4.file ...

By default, a separate folder will be created for each specified coordinate file and placed in the current working directory (see :ref:`submit folder structure` for more information).
If this is not desired, an output directory can be specified with the ``-O`` (or ``--output``) option, in which case a folder will be created for each coordinate file in the given directory:

.. code-block:: console

	$ silico sub coordinate.file  -c 1/1/1 -O "Output Folder"
	
Further options are detailed below.


Reference
---------

.. argparse::
   :module: silico.program.main
   :func: get_argparser
   :prog: silico
   :path: submit
   :noepilog:
   :nodescription:

   