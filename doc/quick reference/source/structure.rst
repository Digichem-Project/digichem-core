Calculation Structure
=====================

When submitting calculations, Silico automatically generates a certain folder structure to ensure individual calculations are stored separately.
This folder structure is described by the following diagram, which shows an example structure for the molecule `Benzene`::

    Benzene
        Combined Reports
            Organic Emitter TDA
        Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p)
            Flags
            Input
            Logs
            Output
            Report
            Result
        Gaussian 16 Optimisation Frequencies PBE1PBE (GD3BJ) Toluene 6-31G(d,p) 02
            Flags
            Input
            Logs
            Output
            Report
            Result
        Gaussian 16 Excited States TDA 10 Singlets 10 Triplets PBE1PBE (GD3BJ) Toluene 6-31G(d,p)
            Flags
            Input
            Logs
            Output
			
The top-most directory, called the `molecule directory`, is named after the molecule on which the calculations have been performed, in this case, `Benzene`. This name is taken directly from the name of the coordinate file that is submitted to Silico without consideration of the actual geometry of the contained molecule. This allows the user greater freedom over the name chosen for the structure, but also implies that molecules with different geometry can be contained within one 'molecule directory' if multiple coordinate files with the same name but different structure are submitted.

Within each molecule directory there are a number of `calculation directories`, which contain data specific to each calculation performed. In the above example there are three calculation directories, two DFT optimisations and one TDA excited states calculation. You will notice that the two optimisation calculations appear to be identical but are still stored in different folders (the second has the numerical suffix `02` added). **This is because Silico always guarantees that each calculation is performed in a separate directory,** so there is no risk of one calculation overwriting data from another.

Within each calculation directory there are a number of sub-folders which contain data specific to that calculation, many of which are self-explanatory. These sub folders are described in detail below:

Flags
    Contains file-flags for the calculation; text files which are created and destroyed by Silico as the calculation progresses. The name of each of these files thus conveys information on the status of the calculation, for example if the calculation is currently running or is complete.
    
Input
    Contains the input file(s) for the calculation. The file(s) will be fully prepared for the calculation program, so they can be inspected to determine the specific parameters of the calculation.
    
Logs
    Contains log files written by Silico, but not by the calculation program itself.
    
Output
    Contains output produced by the calculation program, including log files, but not those written directly by Silico.
    
Report
    Contains a PDF report of the completed calculation, as well as any images or graphs rendered by Silico.
    
Result
    Contains the parsed results of the calculation in various text formats.
    
The `Report` and `Result` folders will only be created once the calculation has completed, and so will not appear while the calculation is running.


Combined Reports
----------------

In addition to calculation directories, each molecule directory may also contain a `Combined Reports` folder. This folder contains compound reports; those generated from the combination of multiple individual calculation results. For instance, the above example contains one combined report in the `Organic Emitter TDA` directory, which contains a combination of results from both the Optimisation and Excited States calculations. These reports are generated automatically when multiple methods are submitted in series.
    