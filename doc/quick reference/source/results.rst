Calculation Results
===================

After the completion of a calculation, Silico will parse the raw calculation data (usually provided in the log file written by the CC program) and processes it in various ways which makes analysis significantly easier. This process includes presenting the data in readable formats, the conversion of results to convenient units (eV, nm, cm\ :sup:`-1` etc), the calculation of simple derived results (HOMO/LUMO energy gap, singlet-triplet splitting energy etc) as well as more involved post-calculations (natural-transition orbitals, orbital density plots, spin-orbit coupling values, simulated spectra etc). The results processed in this way are saved in two separate types of files; these are text-based `result files` and PDF `report files`. The handling of calculation results is described in more detail below.


.. _Text Result Files:

Result Files
------------

The processed calculation results are saved to a number of individual text files in the `Results` directory, an example of which is shown in the following diagram for a calculation performed on the molecule 'Benzene'::

    Results
        Benzene.absorptions.csv
        Benzene.atoms.csv
        Benzene.ES.csv
        Benzene.orbitals.csv
        Benzene.SCF.csv
        Benzene.si
        Benzene.summary.csv
        Benzene.summary.txt
        Benzene.TDM.csv
        Benzene.transitions.csv
        Benzene.UV-Vis.csv
        Benzene.xyz
        
Each file contains data pertaining to a certain part of the calculation results. The .si and .xyz files (`Benzene.si` and `Benzene.xyz` in the above example) hold the output geometry of the calculation in the silico universal input format and xyz format respectively. Either of these files can be specified to ``silico sub`` to submit new calculations based on this geometry, if desired. The .summary.txt file (`Benzene.summary.txt`) contains a summary of the overall calculation results in a simple text format that can be read using any text-processor (notepad, ``less`` etc). The remainder of the result files are in comma-separated values (CSV) format which can be read by spreadsheet processing software (such as Microsoft Excel, LibreOffice Calc or Google Sheets). The data contained in these csv files are described in the following table:

.. csv-table:: Available Result Files
   :file: results.csv
   :widths: 20 80
   :header-rows: 1
   :class: longtable


.. _Report Files:

Report Files
------------

A detailed summary of the calculation results is also saved in portable document format (PDF) to the `Reports` directory. This file contains similar data to the .summary.txt result file, but with greater analysis and detail and in a format that is more easily read on most personal machines. In addition, this graphical format also permits the inclusion of graphs and figures which are generated automatically, such as orbital density plots and simulated spectra. In many cases, all of the results a scientist might be interested in from a calculation can be obtained from the report file alone.

.. note::
    All of the graphs and figures used in the report file are also saved as seperate files in the `Report/image` directory, permitting their inclusion in whatever document the user may require.

