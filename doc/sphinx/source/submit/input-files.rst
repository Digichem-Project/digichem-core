Coordinate and Method Files
=========================================

Traditionally, each computation chemistry (CC) program used to perform calculations (Gaussian, Turbomole, Orca etc.) is associated with a particular (and often unique) input file type.
Generally, these files are made up of two sections which dictate how the calculation is performed:

 * The Atomic Coordinates, specifying the positions of the elements that make up the system (molecule) under study.
 * The Calculation Commands, specifying the nature of the calculation to perform, such as the method, functional, basis set and associated options.
 
In addition to generally being unique to each CC program, this design pattern can make submission of large numbers of calculations difficult.
This is both because it is tedious and time-consuming to duplicate the calculation commands across multiple files, but also because accidental differences can be introduced across the series of files which might not be detected until after the calculation(s) have completed, which may in extreme cases take days or more.

Silico aims to alleviate these traditional problems by separating the atomic coordinates and calculation commands into separate coordinate and method files respectively,
and by providing a standard interface for both by which all CC programs can be accessed. Additionally, Silico supports the submission of multiple (near infinite, limited only by the resources of the operating system) coordinate files to a single method simultaneously, facilitating large-scale computational screens with ease.


Coordinate Files
----------------

When submitting a calculation, the molecules or systems to study are specified by coordinate files, one molecule per file.
Silico supports a wide range of coordinate formats, mostly by interfacing to the obabel file conversion library.\ :cite:p:`Openbabel` This includes, but is not limited to:

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


Charge and Multiplicity
_______________________

In addition to elements and their positions, coordinate files in Silico terminology also convey information on the number and occupation of the electrons of the system.
In traditional cheminformatic style, this is represented by the `charge` and `multiplicity`, which are both integers with the following meaning:

 * Charge: The difference in the total number of electrons compared to the total number of protons of the system. Thus a cation with one fewer electrons than protons has a charge of +1, while an anion with one more electrons than protons has a charge of -1.
 * Multiplicity: A measure of the number of unpaired electrons in the system, where multiplcicity, :math:`m = n + 1` (where n is the number of unpaired electrons).
 
.. note::
	It naturally follows that some combinations of charge and multiplicity are impossible, but this depends on the system in question. For example, so-called superoxide O\ :subscript:`2`\ :superscript:`-` has an odd number of electrons (17), and so must have at least one unpaired. Thus a charge of -1 and multiplicity 1 is impossible for O\ :subscript:`2`\ . Silico does not currently check that any given charge and multiplicity combination is valid; but any CC program almost certainly will.
	
However, some common coordinate formats don't support charge or multiplicity information directly (such as the `.xyz` format).
In this case, it is recommended to first convert the given coordinate format to one that does (such as the Silico universal input format).


Universal Input Format (.si)
____________________________

Silico supports a text-based, program-independent input format known as the silico input format (.si).
To create a .si file, use the ``silico convert`` command to convert any file format supported by Silico:

.. code-block:: console

	$ silico con coordinate.file -O coordinate.si
	
Explicit charge and multiplicity information can be specified by the ``-C`` (or ``--charge``) and ``-M`` (or ``--multiplicity``) options:

.. code-block:: console 

	$ silico con coordinate.file -O coordinate.si -C 0 -M 1

The .si format is written in yaml and has the following basic structure::

	name: null
	charge: 0
	multiplicity: 1
	geometry: |-
	  C          -1.73906         3.58846        -1.30468
	  C          -0.74178         3.28843        -2.23496
	  C          -1.96277         2.73749        -0.21917
	  C          -1.18754         1.58343        -0.06306
	  C           0.03510         2.13524        -2.08167
	  C          -0.18716         1.28164        -0.99543
	  H          -2.33756         4.48085        -1.42656
	  H          -2.73647         2.97415         0.49883
	  H          -0.57443         3.95177        -3.07253
	  H           0.80756         1.90409        -2.80343
	  H           0.41535         0.38931        -0.87850
	  H          -1.36380         0.92662         0.77916

These options have the following meaning:

:name: Optional name of the system. If not given (or ``null``), the name of the file will be used instead.	
:charge: Explicit charge of the system. If  not given (or ``null``), a guess will be used (probably of `0`).
:multiplicity: Explicit multiplicity of the system. If  not given (or ``null``), a guess will be used (probably of `1`).
:geometry: The molecular geometry in .xyz format.

Any of these options can be edited as desired (for example, with the ``nano``, ``vi`` or ``emacs`` editors).
This is particularly useful for changing the charge and/or multiplicity of the system, but coordinates and elements can also be changed as necessary.


Methods
------------

The details of the calculation to be performed (method, functional, basis set etc) are specified by `methods`.
Each method, conceptually, contains three parts which together control how the calculation will be performed, which are:

 * The Destination: A logical or physical location where the calculation will be performed, for example a specific SLURM partition.
 * The Program: A CC progam to perform the calculation, for example Gaussian or Turbomole.
 * The Calculation: A specific set of calculation instructions, including, for example, the method, functional and basis set.

Silico contains a large library of such methods built in  (which can be configured by the administrator of the installation), and for most users this internal database will contain more than sufficient methods to choose from.
However, users can also, if they wish, write their own method files.


Method Files
____________

Method files are written in yaml format and contain three basic `keys` (``destination``, ``program`` and ``calculation``), each of which contains information about the three logical parts of the method (the `Destination`, the `Program` or the `Calculation`).
Each of these structures can either contain a custom definition (essentially defining a new method), or refer to part of a method already built into Silico.
This is useful because it allows a method file to use a built in `destination` and `program` definition, which typically depend on the server setup and cannot be changed anyway, while still changing the details of the `calculation` itself.

To refer to a built in method part, specify either the unique code or ID of the relevant part, for example:

.. code-block:: yaml

	destination: SLURM

or:

.. code-block:: yaml

	destination: 1

If the method part is built up from a hierarchy of TAG names, the path can be specified as a list:

.. code-block:: yaml
	
	program: [Gaussian, 16]
	
or:

.. code-block:: yaml

	program:
	    - Gaussian
	    - 16

Any other format will be interpreted as specifying a new method part, in which case sub-keys should be given for the relevant options that are to be set:

.. code-block:: yaml

	calculation:
	    class_name: Gaussian
	    memory: 1GB
	    name: New Calculation
		
Each method file requires all three sections to be set, so a minimal example for a custom Gaussian calculation might look like the following:

.. code-block:: yaml

	calculation:
	    DFT:
	        functional: B3LYP
	    basis_set:
	        internal: 6-31G(d,p)
	    class_name: Gaussian
	    memory: 1GB
	    name: New Calculation
	destination: Single Node SLURM
	program: Gaussian 16
	
There are a great many options that can be set to finely control the specifics of a calculation; see the Method Reference for a full description of the available options.
