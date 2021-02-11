
# Silico imports.
from silico.report.turbomole import Turbomole_report
from silico.report.gaussian import Gaussian_report
from silico.report.base.pdf import PDF_report
from silico.result.result import Result_set
from silico.exception import Silico_exception
from silico.parser import get_parser
from silico.result.alignment.base import Alignment

def class_from_result(result):
    """
    Get an appropriate report class based on a result object.
    """
    if result.metadata.package == "Turbomole":
        return Turbomole_report
    elif result.metadata.package == "Gaussian":
        return Gaussian_report
    else:
        return PDF_report
    
def from_files(log_file, *, discover_additional_inputs = True, alignment_class_name = None, options, emission_excited_state = None, **named_input_files):
    """
    Convenience function to load a report object from a given calculation log file.
    
    :param log_file: Path to a log file to read.
    :param discover_additional_inputs: If True, the directory of log_file will be searched for additional input files.
    :param alignment_class_name: Name of the alignment class to use.
    :param options: A silico Config dictionary which contains various options that control the appearance of this report.
    :param named_input_files: Additional input files.
    """
    # First, decide on which alignment class we're using.
    alignment_class = Alignment.from_class_handle(options['alignment'] if alignment_class_name is None else alignment_class_name)
    
    # Load results.
    result = get_parser(log_file).process(alignment_class)
    
    # Load emission results if we're given file names instead of results.
    emission_results = {"emission_excited_state": emission_excited_state}
    for emission in ['vertical_emission_ground_result', 'adiabatic_emission_ground_result', 'emission_excited_result']:
        if emission in named_input_files:
            if named_input_files[emission] is not None and not isinstance(named_input_files[emission], Result_set):
                # This emission 'result' is not a result (assume it is a path); try and load it.
                try:
                    emission_results[emission] = get_parser(named_input_files[emission]).process(alignment_class)
                except Exception:
                    raise Silico_exception("Error loading emission result file '{}'".format(named_input_files[emission]))
            else:
                emission_results[emission] = named_input_files[emission]
                
            # Remove from named_input_files.
            named_input_files.pop(emission)
    
    # Add emission energies to result.
    result.add_emission(
        **emission_results
    )
    
    # Get an appropriate report class.
    report_class = class_from_result(result)
    
    # Remove any named_input_files that are None.
    named_input_files = {named_file:named_input_files[named_file] for named_file in named_input_files if named_input_files[named_file] is not None}
    
    # Now continue in subclass.
    return report_class.from_result(result, log_file_path = log_file if discover_additional_inputs else None, options = options, **named_input_files)
    
