Introduction to Using Silico
============================


The Sub-Programs
----------------
The various functionality of Silico is split into a number of sub-programs. These can be accessed
by specifying the name of the desired sub-program after the silico command. For example, to use
the submit program:

.. code-block:: console

    $ silico submit

Each sub-program also has a 3-letter and 1-letter short code by which it can be accessed.
There is no functional difference between using the short codes or the full name, so the user is
encouraged to use whichever style they prefer. For the purposes of this document, all examples
will use the 3-letter short code (unless otherwise stated). As example, the submit program can also
be accessed by any of the following commands:

.. code-block:: console

    $ silico sub

or:

.. code-block:: console

    $ silico s

A list of the available sub-programs is given below:

.. list-table::
    :widths: 17 17 17 49
    :header-rows: 1

    * - Command
      - 3-letter
      - 1-letter
      - Description
    * - ``submit``
      - ``sub``
      - ``s``
      - Submit calculations
    * - ``report``
      - ``rep``
      - ``r``
      - Generate reports from completed calculations
    * - ``result``
      - ``res``
      - ``R``
      - Analyse and tabulate calculation results
    * - ``convert``
      - ``con``
      - ``c``
      - Convert different input/output file formats
    * - ``status``
      - ``sta``
      - ``s``
      - Check how busy the queue is
    * - ``interactive``
      - ``int``
      - ``I``
      - Run Silico interactively

.. note ::
    All silico commands, inlcuding the short codes, are cAsE sEnSiTiVe.
    In particular, all full command names as well as the 3-letter short codes should always be in lower case,
    but certain 1-letter short codes can refer to two different programs depending on capitalization,
    such as ``silico r`` (short for ``silico report``) and ``silico R`` (short for ``silico result``).

.. _Running Interactively:

Running Interactively
---------------------
Most of the Silico sub-programs (with the notable exception of ``silico status``) have both a non-interactive text/console interface as well as an interactive interface powered by the urwid library. Both interfaces are accessed entirely through the console and thus neither requires a running graphical user interface (GUI) (such as X-Server) to function. This makes accessing a remote installation of silico trivial, for example *via* the Secure Shell Protocol (SSH) or similar setup. By default, each sub-program will run in the non-interactive mode. This is generally faster and more convenient for experienced users, as well as being  easily incorporated into custom scripts or programs. However, many users, particularly those who are not accustomed or comfortable with computing in a GUI-free environment, will find the interactive interface easier to use and generally more forgiving.

To use the interactive interface for a sub-program, specify the ``-I`` (or ``--interactive``) option after the sub-program name. For example, to use the submit sub-program in an interactive fashion:

.. code-block:: console

	$ silico sub -I
	
.. note::
	As with all aspects of using silico, there is no difference between specifying the 3-letter short code, 1-letter short code or full command name when specifying the ``-I`` option.
	
In addition, the ``silico interactive`` command can be used to launch Silico interactively without specifying a sub-program. This will open the Silico main menu from which all aspects of the program can be accessed interactively.
