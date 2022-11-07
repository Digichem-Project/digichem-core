Introduction to Using Silico
============================


The Sub-Programs
----------------
The various functionality of Silico is split into a number of subprograms. These can be accessed
by specifying the name of the desired subprogram after the silico command. For example, to use
the submit program:

.. code-block:: console

    $ silico submit

Each subprogram also has a 3-letter and 1-letter short code by which it can be accessed.
There is no functional difference between using the short codes or the full name, so the user is
encouraged to use whichever style they prefer. For the purposes of this document, all examples
will use the 3-letter short code (unless otherwise stated). As example, the submit program can also
be accessed by any of the following commands:

.. code-block:: console

    $ silico sub

or:

.. code-block:: console

    $ silico s

A list of the available subprograms is given below:

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
Most of the Silico subprograms (with the notable exception of ``silico status``) support two interfaces. The first is a non-interactive command-line interface, while the second is an interactive, graphical interface powered by the urwid library. By default, each subprogram will run in non-interactive, command-line mode. To use the interactive interface instead, specify the ``-I`` (or ``--interactive``) option after the subprogram name. For example, to use the submit subprogram in an interactive fashion:

.. code-block:: console

    $ silico sub -I
    
In addition, the ``silico interactive`` command can be used to launch Silico interactively without specifying a subprogram. This will open the Silico main menu from which all aspects of the program can be accessed interactively.
