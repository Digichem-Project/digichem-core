# General imports.
from pathlib import Path

# Silico imports.
from silico.parser.base import Parser
import itertools


class Turbomole_parser(Parser):
    """
    Top level class for parsing output from Turbomole files.
    """
        
#     @classmethod
#     def find_logs(self, log_file):
#         """
#         Find output (log) files from a given hint.
#         
#         :param log_file: A path to a file to use as a hint to find additional log files. log_file can be a directory, in which case files inside this directory will be found.
#         """
#         log_files, parent = super().find_logs(log_file)
#                 
#         # Try and find job files.
#         # These files have names like 'job.0', 'job.1' etc, ending in 'job.last'.
#         for number in itertools.count():
#             # Get the theoretical file name.
#             job_file_path = Path(parent, "job.{}".format(number))
#             
#             # See if it exists (and isn't the log_file given to us).
#             if job_file_path.exists():
#                 # Add to list.
#                 log_files.append(job_file_path)
#             else:
#                 # We've found all the numbered files.
#                 break
#                     
#         # Look for other files.
#         for maybe_file_name in ("basis", "control", "mos", "alpha", "beta", "coord", "gradient", "aoforce", "job.last"):
#             maybe_file_path = Path(parent, maybe_file_name)
#             
#             if maybe_file_path.exists():
#                 # Found it.
#                 log_files.append(maybe_file_path)
#                  
#         # Done.
#         return log_files, parent
    
        
    def parse_metadata(self):
        """
        Parse additional calculation metadata.
        """
        # Add name.
        if self.name is not None:
            self.data.metadata['name'] = self.name
                

    