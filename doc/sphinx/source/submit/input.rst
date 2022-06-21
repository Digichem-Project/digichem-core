Input Files
===========

When submitting calculations with Silico, there are two types of input file that control how each calculation will be performed.
These are: (1) coordinate files, which specify the elements and starting geometry of each of the molecules or systems under study,
and (2) method files, which control the specifics of each calculation, such as the functional and basis set. Of these, input files have to be provided by the user, because Silico cannot know in advance which molecules you will find interesting. However, method files are optional, because Silico contains a large database of methods built-in, and in many cases the methods in this database will be sufficient (particularly for most common DFT problems). However, for greater control over the calculation, or for more unusual or niche calculations, separate method files can also be provided by the user during submission. This section details the format and content of these input files types.

.. _coordinate_files:

Coordinate Files
-----------------

Silico can read input coordinates from a wide range of cheminformatics formats. This includes familiar 3D formats such as Gaussian input files (.com) and .xyz files, as well as 2D formats such as ChemDraw files. When submitting a calculation with Silico, any of the supported coordinate formats can be used interchangeably; even formats that are specific to certain CC programs may be used to perform a calculation with a different CC program. For example, the Gaussian input format (.com, .gjf etc) can be used to perform a calculation with Turbomole. This permits the user to use their preferred program to write their input files for all their calculations without having to worry about incompatibilities. The output files from previously completed calculations (.log etc) can also be used directly as input to new calculations, again regardless of whether the program being submitted to is the same as that which wrote the output file.

The choice of coordinate file format is therefore mostly arbitrary, but with two important considerations:

 #. Two-dimensional formats: Formats that are naturally two-dimensional, such as ChemDraw and similar files, lack information about the z-axis except perhaps for basic stereochemistry (bold and dashed bonds etc). Most molecules, except those that are perfectly planar, have at least some extent in the z-axis and so a conversion to a 3D format must be undertaken to provide a reasonable starting geometry for further optimisation. In Silico, this conversion utilises a rapid molecular-mechanics (MM) optimisation provided by the obabel library\ :cite:p:`Openbabel`. In many cases this will result in a satisfactory structure, but occasionally the geometry will become locked in an impossible or high-energy conformation. Similarly, steric information may be destroyed by the optimisation process. In these cases it is recommended to first convert the 2D coordinates to a 3D representation prior to submission, so the starting geometry can be inspected. This can be achieved using the Silico convert subprogram, for example:
 
  .. code-block:: console

        $ silico con molecule1.cdx -O molecule1.cml
  
  See the convert command for more information. Note also that is generally does not make sense to perform calculations directly on a 2D coordinate format without a prior optimisation at a higher level of theory (for example, with DFT).
 #. Charge and multiplicity: Many coordinate formats contain information only about the elements that make up the molecule and their positions, and do not convey information about the number of electrons. In these cases a neutral, singlet ground state is typically assumed. If this is not suitable, then a format which does contain electron information (normally in the form of charge and multiplicity) should be chosen. See the below sections on :ref:`charge and multiplicity <charge_and_mult>` and the :ref:`silico input format <si_format>` for more information.

.. _coordinate_formats :

Formats
________

 * Gaussian and Turbomole output (.log)
 * GaussView input (.com, .gjf, .gjc and .gau)
 * ChemDraw (.cdx and .cdxml)
 * MarvinSketch (.mrv)
 * Independant cheminformatics (.cml, .xyz etc)
 * Crystallographic information file (.cif)


.. _charge_and_mult :

Charge and Multiplicity
_______________________

In addition to elements and their positions, coordinate files in Silico terminology also convey information on the number and occupation of the electrons of the system.
In traditional cheminformatic style, this is represented by the `charge` and `multiplicity` of the system, which are both integers with the following meaning:

 * Charge: The difference in the total number of electrons compared to the total number of protons of the system. Thus a cation with one fewer electrons than protons has a charge of +1, while an anion with one more electrons than protons has a charge of -1.
 * Multiplicity: A measure of the number of unpaired electrons in the system, where multiplcicity, :math:`m = n + 1` (where n is the number of unpaired electrons).
 
.. note::
    It naturally follows that some combinations of charge and multiplicity are impossible, but this depends on the system in question. For example, so-called superoxide O\ :subscript:`2`\ :superscript:`-` has an odd number of electrons (17), and so must have at least one unpaired. Thus a charge of -1 and multiplicity 1 is impossible for O\ :subscript:`2`\ . Silico does not currently check that any given charge and multiplicity combination is valid; but any CC program almost certainly will.
    
However, some common coordinate formats don't support charge or multiplicity information directly (such as the `.xyz` format).
In this case, it is recommended to first convert the given coordinate format to one that does (such as the Silico universal input format).

.. _si_format :

Universal Input Format (.si)
____________________________

Silico supports a text-based, program-independent input format known as the silico input format (.si).
To create a .si file, use the ``silico convert`` command to convert any file format supported by Silico:

.. code-block:: console

    $ silico con coordinate.file -O coordinate.si
    
Explicit charge and multiplicity information can be specified by the ``-C`` (or ``--charge``) and ``-M`` (or ``--multiplicity``) options:

.. code-block:: console 

    $ silico con coordinate.file -O coordinate.si -C 0 -M 1

The .si format is written in yaml, so a silico input file can also be written from scratch using any text editor.  The .si file has the following basic structure::

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

Method files are written in yaml format and contain three basic `keys` (``destination``, ``program`` and ``calculation``), each of which contains information about the three logical parts of the method (the `Destination`, the `Program` or the `Calculation`):

.. code-block:: yaml

    destination:        # Destination (SLURM partition, storage location etc) information.
    program:            # CC program (Gaussian, Turbomoel etc) information.
    calculation:        # Specific calculation options (functional, method, basis set etc).

Each of these structures can either contain a custom definition (essentially defining a new method), or refer to part of a method already built into Silico.
This is useful because it allows a method file, for example, to use a built in `destination` and `program` definition, which typically depend on the server setup and cannot be changed anyway, while still changing the details of the `calculation` itself.

To refer to a built in method part, specify either the unique code or ID of the relevant part, for example:

.. code-block:: yaml

    destination: SLURM      # Use the built in destination called 'SLURM'.

or:

.. code-block:: yaml

    destination: 1          # Use the built in destination with code of 1.

If the method part is built up from a hierarchy of TAG names, the path can be specified as a list:

.. code-block:: yaml
    
    program: [Gaussian, 16] # Use the built in program with the name "Gaussian" "16".
    
or:

.. code-block:: yaml

    program:
        - Gaussian
        - 16

Any other format will be interpreted as specifying a new method part, in which case sub-keys should be given for the relevant options that are to be set:

.. code-block:: yaml

    calculation:            # Options for a new type of calculation.
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