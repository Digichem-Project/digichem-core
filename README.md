# Silico: Computational Chemistry Simplified

A software toolkit that aims to reduce the complexity of computational chemistry, particularly for newcomers,
by automating and/or simplifying the computational pipepline.

## Features

 - Handling of the entire submission process from start to finish. The user does not not need to write batch/submission files or worry about calling the CC programs themselves.
 - Separation of the concepts of ‘starting geometries’ (coordinates) and ‘calculation options’ (methods) into separate logical units.
 - Support for ca. 150 different input coordinate formats. Any of these formats can be used interchangeably with any of the supported CC programs.
 - Introduction of the concept of calculation ‘methods’ which define all aspects of a computation (except any starting geometries). Support for specifying methods as individual files (in program-independent, YAML format) or by selecting from an internal database.
 - Support for submitting multiple input coordinates simultaneously, facilitating large scale computational screens with ease.
 - Support for queueing multiple methods sequentially to support compound jobs in which one calculation depends on the output geometry of another.
 - Automated parsing and analysis of completed calculation results and saving of the resulting data to formats which can be easily read and processed by the user, including plain text and comma-separated values (CSV).
 - Automation of common post-processing tasks, including the generation of 3D structure images, orbital density plots, natural-transition orbitals, difference density plots and simulated spectra.
 - Automated generation of a summary report from completed calculations. This report is saved in portable-document format (PDF) for easy reading and distribution.
 
## Support
Silico currently supports the following:

#### Job Scheduling
 - SLURM

#### Computational Packages
 - Gaussian 09
 - Gaussian 16
 - Turbomole
 
#### Python
 - Python >= 3.6 (and above)
 
#### Operating Systems
 - Linux
