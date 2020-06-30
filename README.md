# Silico (beta)

A software package for managing and automating various aspects of computational chemisty, Silico provides programs and Python 3 libraries that support:
 - Submission to computational programs through a unified interface
 - Simultaneous submission of multiple molecules/systems
 - Automatic in series submission of results from completed calculations to subsequent calculations
 - Automatic conversion of input files to types appropriate for the selected computational program via integration with [Openbabel v3](https://github.com/openbabel/openbabel)
 - Automatic and manual analysis of computation results, including tabulation to csv format of multiple results, via integration with [cclib](https://github.com/cclib/cclib/)
 - Automatic and manual generation of pdf reports from computation results, including rendered 3D structures, orbital images and graphs, via integration with [Weasyprint](https://weasyprint.org/), [VMD](https://www.ks.uiuc.edu/Research/vmd/) and [Matplotlib](https://matplotlib.org/)
 
## Support
Silico currently supports the following:

#### Job Scheduling
 - SLURM

#### Computational Packages
 - Gaussian 09
 - Gaussian 16
 
#### Python
 - Python 3.6 (and above)
 - Python < 3.6 (or missing entirely) through PyInstaller frozen packages (programs only)
 
#### Operating Systems
 - Linux
