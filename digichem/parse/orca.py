from silico.parse.base import Parser
import silico.file.types as file_types


class Orca_parser(Parser):
    """
    Top level class for parsing output from Gaussian log files.
    """
    
    # A dictionary of recognised auxiliary file types.
    INPUT_FILE_TYPES = {
            file_types.orca_gbw_file: "gbw_file",
            file_types.orca_density_file: "density_file",
        }