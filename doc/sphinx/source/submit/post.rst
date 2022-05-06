Post Analysis
=============

Following the completion of a calculation, Silico will automatically perform various post analyses of the calculation results to present the data in a more readable format.
This includes the writing of a number of plain-text and comma-separated value (CSV) text files detailing the results of the calculation, as well as a portable-document format (PDF) report for easy reading and distribution of results.


.. _result_files:

Result Files
------------

Text result files are written to the `Results` folder.
This includes a summary file, which can be easily read in a text editor (such as ``less``) to obtain a summary of the most important calculation results:

.. code-block:: console

	$ less Results/MOLECULE.summary
	
Where `MOLECULE` in the above command will be the name of the coordinate file that was actually submitted.
In addition to this summary file, detailed tables of results are also available in CSV format, which can be easily imported into a spreadsheet program for further analysis.
The currently available data files are as follows:

.. list-table::
    :widths: 20 80
    :header-rows: 1

    * - File
      - Description
    * - MOLECULE.atoms.csv
      - Final atom coordinates.
    * - MOLECULE.orbitals.csv
      - Orbital numbers, labels and energies (restricted calculations only).
    * - MOLECULE.alpha.csv
      - Alpha orbital numbers, labels and energies (unrestricted calculations only).
    * - MOLECULE.beta.csv
      - Beta orbital numbers, labels and energies (unrestricted calculations only).
    * - MOLECULE.SCF.csv
      - SCF (self consistent field) energy. For optimisations, includes the energy at each step which can be used to graph the calculation’s convergence. The SCF energy is typically the HF (Hartree–Fock) or DFT energy.
    * - MOLECULE.MP.csv
      - Same as MOLECULE.SCF.csv, but for Møller–Plesset energies. The energies are total (including HF and MP correction) and are from the highest MP level of the calculation.
    * - MOLECULE.CC.csv
      - Same as MOLECULE.SCF.csv, but for coupled-cluster energies. The energies are total (including HF and CC correction, if applicable).
    * - MOLECULE.ES.csv
      - Excited state results (from TD-DFT, TDA etc). Includes energy, symmetry, multiplicity and orbital contributions.
    * - MOLECULE.transitions.csv
      - Excited state transitions (orbital contributions) in a format that is more convenient for directly comparing transitions than MOLECULE.ES.csv.
    * - MOLECULE.TDM.csv
      - Excited state transition dipole moments.
    * - MOLECULE.UV-Vis.csv
      - A simulated UV-Vis spectrum using Gaussian line-broadening in corrected units of wavelength (nm).
    * - MOLECULE.absorptions.csv
      - A simulated UV-Vis spectrum using Gaussian line-broadening in units of energy (eV).
    * - MOLECULE.vibrations.csv
      - Vibrational frequencies and their intensities.
    * - MOLECULE.IR.csv
      - A simulated infra-red (IR) absorption spectrum using Gaussian line-broadening.
    * - MOLECULE.SOC.csv
      - Singlet/triplet spin-orbit coupling values.
    
Finally, the output geometry of the calculation is also saved to the Results folder, in both .xyz and .si format.
These files can be used to easily start a subsequent calculation based on the same geometry, if so desired.