# General imports.
import pysoc.io.SOC
from datetime import datetime, timedelta

# Silico imports.
from silico.parser.main import Parser
import silico.logging
from silico.exception.base import Result_unavailable_error
import silico.file.types as file_types


class Gaussian_parser(Parser):
    """
    Top level class for parsing output from Gaussian log files.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
            file_types.gaussian_chk_file: "chk_file",
            file_types.gaussian_fchk_file: "fchk_file",
            file_types.gaussian_rwf_file: "rwf_file"
        }
    
    ## Headers that indicate certain data is about to be printed.
    #TDM_HEADER = " Ground to excited state transition electric dipole moments (Au):\n"
    
        
    def parse(self):
        """
        Extract results from our output files.
        """
        # We start by using cclib to get most of the data we need.
        super().parse()
        
        # Output a message (because this is slow).
        silico.logging.get_logger().debug("Secondary parsing calculation result '{}'".format(self.description))
                            
        # Next try and get SOC.
        try:
            self.parse_SOC()
        except Exception:
            silico.logging.get_logger().debug("Cannot parse SOC from output file '{}'".format(self.log_file_path), exc_info = True)
        
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
        super().parse_metadata()
        
        # Assume we used 1 CPU if not otherwise clear (is this a good idea?)
        self.data.metadata['numcpus'] = 1
        
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
            
            
    def parse_SOC(self):
        """
        Parse spin-orbit coupling using PySOC.
        """
        # For SOC, we need both .log and .rwf file.
        # No need to check for these tho; pysoc does that for us.
        # We also need etsyms to decide which excited state is which.
        if not hasattr(self.data, "etsyms"):
            raise Result_unavailable_error("Spin-orbit coupling", "There are no excited states available")
        
        # Get a PySOC parser.
        soc_calculator = pysoc.io.SOC.Calculator(self.log_file_path, rwf_file_name = self.aux_files['rwf_file'])
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
            
            