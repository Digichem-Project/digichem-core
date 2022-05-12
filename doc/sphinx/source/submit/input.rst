Input Files
===========

When using Silico, there are two types of input file that control how each calculation will be performed.
These are: (1) coordinate files, which specify the elements and starting geometry of each of the molecules or systems under study,
and (2) method files, which control the specifics of each calculation, such as the functional and basis set.

One of the main advantages of using Silico for calculation submission is that multiple coordinate files can be specified at once, in which case all the given coordinate files will be submitted to the same calculation (as specified by the method file) simultaneously.
In addition, multiple method files can also be specified, in which case each coordinate file will first be submitted to the calculation defined by the first method file specified. Once each of these calculations has completed, the resulting atomic geometry will then be submitted to the calculation defined by the second method file, and so on until all method files have been exhausted.

Silico comes pre-loaded with a large database of method files, and so in most cases these do not need to be written by the user.
Hence all that is generally needed to submit a calculation is a number of coordinate file(s) specifying the systems of interest.

Input Coordinates
-----------------

Silico supports a wide variety of coordinate formats, including both 2D and 3D formats, each of which can be written by different programs. Notable entries include:

 * Gaussian and Turbomole output (.log)
 * GaussView input (.com, .gjf, .gjc and .gau)
 * ChemDraw (.cdx and .cdxml)
 * MarvinSketch (.mrv)
 * Independant cheminformatics (.cml, .xyz etc)
 * Crystallographic information file (.cif)
 
Silico uses the OpenBabel library for file conversion. Please see `the OpenBabel documentation <https://open-babel.readthedocs.io/en/latest/FileFormats/Overview.html>`_ for a full list of supported formats.

.. note::
    Care should be taken when using 2D formats, particularly for complex 3D structures or those with specific steric information (enantiomers, for example).
    The conversion from 2D to 3D employs a rapid molecular-mechanics (MM) optimisation provided by the obabel library\ :cite:p:`Openbabel`. In many cases this will result in a satisfactory starting structure for further optimisation, but occassionally the geometry will become locked in an impossible or high-energy conformation. Similarly, steric information may be destroyed by the optimisation process. In these cases it is recommended to first convert the 2D coordinates to a 3D representation using the ``silico convert`` subprogram and manually inspect the resulting geometry prior to submission.
    
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