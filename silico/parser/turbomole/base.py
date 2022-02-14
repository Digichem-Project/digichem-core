# General imports.
from cclib.io.ccio import sort_turbomole_outputs
import re
from datetime import timedelta, datetime

# Silico imports.
from silico.parser.base import Parser


class Turbomole_parser(Parser):
    """
    Top level class for parsing output from Turbomole files.
    """
    DAYS_REGEX = re.compile(r"([0-9.]*) days")
    HOURS_REGEX = re.compile(r"([0-9.]*) hours")
    MINUTES_REGEX = re.compile(r"([0-9.]*) minutes")
    SECONDS_REGEX = re.compile(r"([0-9.]*) seconds")
    
    def duration_to_timedelta(self, duration_str):
        """
        Convert a Turbomole duration string into an equivalent timedelta object.
        """
        time_parts = {'days': 0, 'hours': 0, 'minutes': 0, 'seconds': 0}
        
        for time_part in time_parts:
            # Use regex to look for each part in the string.
            match = getattr(self, time_part.upper() + '_REGEX').search(duration_str)
            if match:
                time_parts[time_part] = float(match.group(1))
                
        # Build a timedelta from our parts.
        duration = timedelta(days = time_parts['days'], hours = time_parts['hours'], minutes = time_parts['minutes'], milliseconds = time_parts['seconds'] * 1000)
        
        # All done.
        return duration
        
    
    def parse_metadata(self):
        """
        Parse additional calculation metadata.
        """
        super().parse_metadata()
                    
        # Look for duration information.
        self.data.metadata['walltime'] = 0.0
        self.data.metadata['cputime'] = 0.0
        
        # For Turbomole, we have multiple files to look through.
        # Fortunately, cclib knows which order to process in.
        for log_file_path in sort_turbomole_outputs(self.log_file_paths):
            with open(log_file_path, "rt") as log_file:
                for line in log_file:
                    
                    # Look for duration.
                    if "total  cpu-time :" in line:
                        self.data.metadata['cputime'] += self.duration_to_timedelta(line).total_seconds()
                    if "total wall-time :" in line:
                        self.data.metadata['walltime'] += self.duration_to_timedelta(line).total_seconds()
                        
                    # And also end date.
                    if ": all done  ****" in line:
                        # Skip 2 lines.
                        line = next(log_file)
                        line = next(log_file)
                        line = next(log_file)
                        try:
                            self.data.metadata['date'] = datetime.strptime(line.strip(), "%Y-%m-%d %H:%M:%S.%f").timestamp()
                        except ValueError:
                            # We couldn't parse.
                            pass
        
        # Delete our durations if they are zero.
        if self.data.metadata['walltime'] == 0.0:
            del(self.data.metadata['walltime'])
            
        if self.data.metadata['cputime'] == 0.0:
            del(self.data.metadata['cputime'])
                    
                

    