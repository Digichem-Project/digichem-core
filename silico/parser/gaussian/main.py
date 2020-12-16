# General imports.
from pathlib import Path
import pysoc.io.SOC
from logging import getLogger
from datetime import datetime, timedelta

# Silico imports.
from silico.parser.base import Parser
import silico
from silico.exception.base import Result_unavailable_error


class Gaussian_parser(Parser):
    """
    Top level class for parsing output from Gaussian log files.
    """
    
    # Headers that indicate certain data is about to be printed.
    TDM_HEADER = " Ground to excited state transition electric dipole moments (Au):\n"
    
    def __init__(self, log_file, *, rwf_file = None, **kwargs):
        """
        Constructor for Gaussian output file parsers.
        """
        self.rwf_file_path = rwf_file
        
        super().__init__(log_file, **kwargs)
        
    @classmethod
    def from_log(self, log_file, **kwargs):
        """
        Intelligent constructor that will attempt to guess the location of files from a given log file.
        
        :param log_file: Output file to parse.
        """
        log_file = Path(log_file)
        
        # See if we can find our rwf file.
        if log_file.with_suffix(".rwf").exists():
            kwargs['rwf_file'] = log_file.with_suffix(".rwf")
            
        # Continue with normal constructor.
        return self(log_file, **kwargs)
        
    def parse(self):
        """
        Extract results from our output files.
        """
        # We start by using cclib to get most of the data we need.
        super().parse()
        
        # Output a message (because this is slow).
        getLogger(silico.logger_name).debug("Secondary parsing calculation result '{}'".format(self.description))
                    
        # Date and time.
        self.parse_metadata()
        
        # Next try and get SOC.
        try:
            self.parse_SOC()
        except Exception:
            getLogger(silico.logger_name).debug("Cannot parse SOC from output file '{}'".format(self.log_file_path), exc_info = True)
        
        # Then TDM.
        try:
            self.parse_transition_dipole_moments()
        except Exception:
            getLogger(silico.logger_name).debug("Cannot parse TDM from output file '{}'".format(self.log_file_path), exc_info = True)
        
        # All done.
        return self.data
    
    # Headers for date strings.
    DATE_HEADER = "Normal termination of"
    ELAPSED_TIME_HEADER = "Elapsed time:"
    CPU_TIME_HEADER = "Job cpu time:"
    CPU_HEADER = "Will use up to"
    
    def parse_metadata(self):
        """
        Parse additional calculation metadata.
        """
        # Add name.
        if self.name is not None:
            self.data.metadata['name'] = self.name
        
        # Open our file.
        with open(self.log_file_path, "rt") as log_file:
            # Loop through our lines, looking for a specific line of text.
            for log_file_line in log_file:
                # Although we only need the last ~5 lines from the (possibly huge) log file, we read all the way through because negative seek()ing is tricky.
                # Look for our key string.
                if self.DATE_HEADER in log_file_line:
                    # This line looks like: "Normal termination of Gaussian 16 at Sun Dec  6 19:13:09 2020"
                    date_str = " ".join(log_file_line.split()[-4:])
                    self.data.metadata['date'] = datetime.strptime(date_str, "%b %d %H:%M:%S %Y.").timestamp()
                    
                elif self.ELAPSED_TIME_HEADER in log_file_line:
                    # This line looks like: "Elapsed time:       0 days  2 hours 38 minutes 50.9 seconds."
                    datey = log_file_line.split()[-8:]
                    self.data.metadata['walltime'] = timedelta(days = int(datey[0]), hours = int(datey[2]), minutes = int(datey[4]), seconds = float(datey[6])).total_seconds()
                    
                elif self.CPU_TIME_HEADER in log_file_line:
                    # This line looks like: "Job cpu time:       0 days 20 hours 52 minutes 17.3 seconds."
                    datey = log_file_line.split()[-8:]
                    self.data.metadata['cputime'] = timedelta(days = int(datey[0]), hours = int(datey[2]), minutes = int(datey[4]), seconds = float(datey[6])).total_seconds()
                    
                elif self.CPU_HEADER in log_file_line:
                    # This line looks like: "Will use up to   10 processors via shared memory."
                    self.data.metadata['numcpus'] = int(log_file_line.split()[4])
    
    
    def parse_transition_dipole_moments(self):
        """
        Parse transition dipole moments.
        """
        # Each calculation may calculate TDMs multiple times (for example in an optimisation + excited state).
        transition_dipole_groups = []
        
        # Open our file.
        with open(self.log_file_path, "rt") as log_file:
            # Loop through our lines, looking for a specific line of text.
            #found = False
            for log_file_line in log_file:
                # Look for our key string.
                # It might be better to look for a substring rather than the whole string? Need to benchmark.
                if log_file_line == self.TDM_HEADER:
                    # We found our header.
                    
                    # The next line should be the table header, so we can skip it.
                    log_file.readline()
                    
                    # Now loop through, splitting on white space.
                    transition_dipoles = []
                    for log_file_line in log_file:
                        # Split on whitespace.
                        split_line = log_file_line.split()
                        
                        try:
                            # Add into our list.
                            transition_dipoles.append({
                                'state_level': int(split_line[0]),
                                'origin_coords': (0.0, 0.0, 0.0),
                                'vector_coords': (
                                    self.au_to_debye(float(split_line[1])),
                                    self.au_to_debye(float(split_line[2])),
                                    self.au_to_debye(float(split_line[3]))
                                )
                                })
                        except (ValueError, IndexError):
                            # No more data.
                            break
                    
                    # Add this set of transition dipoles to our big list.
                    transition_dipole_groups.append(transition_dipoles)
        
        # If we have no TDMs, get upset.
        if len(transition_dipole_groups) == 0:
            raise Result_unavailable_error("Transition dipole moment", "There is no transition dipole moment data")

        # We're only interested in the last group of TDMs.
        transition_dipoles = transition_dipole_groups[-1]
        
        # Sort our list of TDMs (probably unnecessary).
        transition_dipoles.sort(key = lambda dipole: dipole['state_level'])
        
        # Now split into a list for cclib.
        # This is a list of the vector coords of each TDM (we assumed the origin is 0,0,0).
        self.data.etmoments = [tdm['vector_coords'] for tdm in transition_dipoles]
            
            
            
    def parse_SOC(self):
        """
        Parse spin-orbit coupling using PySOC.
        """
        # For SOC, we need both .log and .rwf file.
        # We also need etsyms to decide which excited state is which.
        if self.rwf_file_path is None:
            pass
            #raise Result_unavailable_error("Spin-orbit coupling", "There is no .rwf file available")
        elif not hasattr(self.data, "etsyms"):
            raise Result_unavailable_error("Spin-orbit coupling", "There are no excited states available")
        
        
        # Get a PySOC parser.
        soc_calculator = pysoc.io.SOC.Calculator(self.log_file_path, rwf_file_name = self.rwf_file_path)
        soc_calculator.calculate()
        SOC_table = soc_calculator.soc_td.SOC
        
        # We'll split the SOC table given to use by PySOC to better match the format used by cclib.
        socstates = []
        soc = []
        
        for SOC_line in SOC_table:
            # Add states.
            socstates.append([SOC_line.singlet_state, SOC_line.triplet_state])
            
            # Add coupling.
            soc.append([SOC_line.positive_one, SOC_line.zero, SOC_line.negative_one])
                
        # Add to data.
        self.data.socstates = socstates
        self.data.soc = soc
            
            