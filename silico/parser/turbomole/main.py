# General imports.
from pathlib import Path
import cclib.io

# Silico imports.
from silico.parser.base import Parser
from logging import getLogger
import silico
import itertools
import numpy
from silico.parser.turbomole.orbitals import Turbomole_orbitals


class Turbomole_parser(Parser):
    """
    Top level class for parsing output from Turbomole files.
    """
    
    def __init__(self,
        log_file,
        basis_file = None,
        control_file = None,
        mos_file = None,
        alpha_file = None,
        beta_file = None,
        job_files = (),
        coord_file = None,
        gradient_file = None,
        aoforce_file = None,
        **kwargs):
        """
        Constructor for Turbomole parsers.
        
        :param log_file: Path to the main .log output file.
        :param basis_file: Optional path to the basis file.
        :param control_file: Optional path to the control file.
        :param mos_file: Optional path to the mos file containing molecular orbital info.
        :param alpha_file: Optional path to the alpha file containing alpha molecular orbital info.
        :param beta_file: Optional path to the beta file containing beta molecular orbital info.
        :param job_files: Optional list of files containing optimisation steps (this list should be in order).
        :param coord_file: Optional path to the coord file.
        :param gradient_file: Optional path to the gradient file.
        :param aoforce_file: Optional path to the aoforce file.
        """
        self.basis_file = basis_file
        self.control_file = control_file
        self.mos_file = mos_file
        self.alpha_file = alpha_file
        self.beta_file = beta_file
        self.optimisation_step_file_paths = job_files
        self.coord_file = coord_file
        self.gradient_file = gradient_file
        self.aoforce_file = aoforce_file
        
        # Pass to next constructor.
        super().__init__(
            log_file,
            basis_file,
            control_file,
            mos_file,
            alpha_file,
            beta_file,
            *job_files,
            coord_file,
            gradient_file,
            aoforce_file,
            **kwargs)
        
    @classmethod
    def from_log(self, log_file, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of files from a given log file.
        
        :param log_file: Output file to parse.
        """
        # First, find our parent dir.
        log_file = Path(log_file)
        parent = log_file.parent
                
        # Try and find job files.
        kwargs['job_files'] = []
        
        # These files have names like 'job.0', 'job.1' etc, ending in 'job.last'.
        for number in itertools.count():
            # Get the theoretical file name.
            job_file_path = Path(parent, "job.{}".format(number))
            
            # See if it exists.
            if job_file_path.exists():
                # Add to list.
                kwargs['job_files'].append(job_file_path)
            else:
                # We've found all the numbered files.
                break
            
        # See if the job.last file exists.
        job_last_path = Path(parent, "job.last")
        
        if job_last_path.exists():
            kwargs['job_files'].append(job_last_path)
            
        # Look for other files.
        for maybe_file_name in ("basis", "control", "mos", "alpha", "beta", "coord", "gradient", "aoforce"):
            maybe_file_path = Path(parent, maybe_file_name)
            
            if maybe_file_path.exists():
                # Found it.
                kwargs["{}_file".format(maybe_file_name)] = maybe_file_path
                
        # Continue with normal constructor.
        return self(log_file, **kwargs)
    
    
    def parse(self):
        """
        Extract results from our output files.
        """
        # First, load results from log file.
        self.data = None
        super().parse()
        
        #self.parse_optimisations()
        # cclib puts the final energy first, so we need to reorganise.
        for energy in ["scfenergies", "mpenergies", "ccenergies"]:
            if hasattr(self.data, energy) and len(getattr(self.data, energy)) > 1:
                # Get the first energy.
                last_energy = getattr(self.data, energy)[0]
                
                # Make a new array.
                setattr(self.data, energy, numpy.append(getattr(self.data, energy)[1:], last_energy))
                #getattr(self.data, energy).append(last_energy)
                
        self.parse_metadata()
        
        # Get a list of subparsers.
        subparsers = [Turbomole_orbitals(self)]
        
        # Read main log file.
        with open(self.optimisation_step_file_paths[-1]) as log_file:
            while True:
                try:
                    line = next(log_file)
                except StopIteration:
                    # Nothing more to read.
                    break
            
                for subparser in subparsers:
                    subparser.parse(line, log_file)
        
    def parse_metadata(self):
        """
        Parse additional calculation metadata.
        """
        # Add name.
        if self.name is not None:
            self.data.metadata['name'] = self.name
            
    
    
    def parse_HOMO_LUMO(self):
        """
        Get our HOMO and LUMO.
        """
        # The HOMO and LUMO energies are stored in the main .log file.
        
        
    # This method is not required; cclib can already parse these files.          
    def parse_optimisations(self):
        """
        Parse each of our optimisation step files to get a list of energies.
        """
        energies = {
            'scfenergies': [],
            'mpenergies': [],
            'ccenergies': []
            }
        
        for optimisation_step_file_path in self.optimisation_step_file_paths:
            
            # Output a message (because this is slow).
            getLogger(silico.logger_name).info("Parsing Turbomole step result '{}'".format(optimisation_step_file_path))
            
            with open(optimisation_step_file_path, "rt") as optimisation_step_file:
                ccdata = cclib.io.ccread(optimisation_step_file)
                
                # Add selected attributes to our real data.
                for energy_name, energy in energies.items():
                    if hasattr(ccdata, energy_name):
                        energy.append(getattr(ccdata, energy_name)[-1])
                            
        # If we got any energies, overwrite those in our main result.
        for energy_name, energy in energies.items():
            if len(energy) > 0:
                setattr(self.data, energy_name, numpy.array(energy))
                

    