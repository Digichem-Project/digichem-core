# Code for processing results from calculations and generating reports.

# General imports.
from datetime import timedelta, datetime

# Silico imports.
from silico.exception import Result_unavailable_error
from silico.result import Result_object
import math
from silico.misc.time import latest_datetime, total_timedelta
import itertools

class Metadata(Result_object):
    """
    Class for storing calculation metadata.
    """
    
    # A dictionary of known computational orders, where the key is name of the method (upper case) and the value is an integer describing the ordered 'level of theory'.
    METHODS = {
        "HF":       1,
        "DFT":      2,
        "LMP2":     3,
        "DF-MP2":   4,
        "MP2":      5,
        "MP3":      6,
        "MP4":      7,
        "MP5":      8,
        "CCSD":     9,
        "CCSD(T)":  10,
        "CCSD-T":   11
    }
    
    def __init__(
            self,
            name = None,
            log_files = None,
            auxiliary_files = None,
            date = None,
            duration = None,
            package = None,
            package_version = None,
            calculations = None,
            success = None,
            methods = None,
            functional = None,
            basis_set = None,
            charge = None,
            multiplicity = None,
            optimisation_converged = None,
            temperature = None,
            pressure = None,
            orbital_spin_type = None):
        """
        Constructor for result Metadata objects.
        
        :param name: Optional name of this calculation result.
        :param log_files: An optional list of text-based calculation log files from which this result was parsed.
        :param auxiliary_files: An optional dict of auxiliary files associated with this calculation result.
        :param num_calculations: Optional number of individual calculations this metadata represents.
        :param date: Optional date (datetime object) of this calculation result.
        :param duration: Optional duration (timedelta object) of this calculation.
        :param package: Optional string identifying the computational chem program that performed the calculation.
        :param package: Optional string identifying the version of the computational chem program that performed the calculation.
        :param calculations: A list of strings of the different calculations carried out (Opt, Freq, TD, TDA, SP etc).
        :param success: Whether the calculation completed successfully or not.
        :param methods: List of methods (DFT, HF, MPn etc) used in the calculation.
        :param functional: Functional used in the calculation.
        :param basis_set: Basis set used in the calculation.
        :param charge: Charge (positive or negative integer) of the system studied.
        :param multiplicity: The multiplicity of the system.
        :param optimisation_converged: Whether the optimisation converged or not.
        :param temperature: The temperature used in the calculation (not always relevant).
        :param pressure: The pressure used in the calculation (not always relevant).
        :param orbital_spin_type: The types of orbitals that have been calculated, either 'restricted' or 'unrestricted'.
        """
        self.num_calculations = 1
        self.name = name
        self.log_files = log_files if log_files is not None else []
        self.auxiliary_files = auxiliary_files if auxiliary_files is not None else {}
        self.date = date
        self.duration = duration
        self.package = package
        self.package_version = package_version
        self.calculations = calculations if calculations is not None else []
        self.success = success
        self.methods = methods if methods is not None else []
        self.functional = functional
        self.basis_set = basis_set
        self.charge = charge
        self.multiplicity = multiplicity
        self.optimisation_converged = optimisation_converged
        self.temperature = temperature
        self.pressure = pressure
        self.orbital_spin_type = orbital_spin_type
        
    @classmethod
    def merge(self, *multiple_metadatas):
        """
        Merge multiple metadata objects into a single metadata.
        """
        return Merged_metadata.from_metadatas(*multiple_metadatas)
        
    @property
    def identity(self):
        """
        A dictionary of critical attributes that identify a calculation.
        """
        # Get our list of methods, replacing 'DFT' with the functional used if available.
        methods = [method if method != "DFT" and self.functional is not None else self.functional for method in self.methods]
        return {
            "package": self.package_string,
            "calculations": self.calculations,
            "methods": methods,
            "basis": self.basis_set,
            "multiplicity": self.multiplicity,
            "charge": self.charge
        }
    
    @classmethod
    def sorted_methods(self, methods):
        """
        Order a list of methods (HF, DFT, MP2 etc) in terms of 'level of theory', with the lowest level (HF or DFT probably) first and highest (MP4, CCSD etc) last.
        """
        return sorted((method.upper() for method in methods), key = lambda method: self.METHODS[method] if method in self.METHODS else math.inf)
    
    @property
    def package_string(self):
        """
        The comp chem package and version combined into one string.
        """
        package_string = self.package if self.package is not None else ""
        
        # Add version string if we have it.
        package_string += " " if package_string != "" else ""
        package_string += "({})".format(self.package_version)
        
        # Done.
        return package_string 
        
    @property
    def calculations_string(self):
        """
        Get the list of calculation types as a single string, or None if there are no calculations set.
        """
        return ", ".join(self.calculations) if len(self.calculations) != 0 else None
    
    @property
    def methods_string(self):
        """
        Get the list of calculation methods as a single string, or None if there are no methods set.
        """
        return ", ".join(self.methods) if len(self.methods) != 0 else None
            
    @classmethod
    def get_calculation_types_from_cclib(self, ccdata):
        """
        Determine what types of calculation we did from the data provided by cclib.
        
        :param ccdata: The data from cclib.
        :return: A list of strings.
        """
        calcs = []
        # Need to think of a better way to determine what is an optimisation, as frequency calcs contain an optdone for some reason.
        #if hasattr(ccdata, 'optdone'):
        if len(getattr(ccdata, 'scfenergies', [])) > 1:
            calcs.append('Optimisation')
        if hasattr(ccdata, 'vibfreqs'):
            calcs.append('Frequencies')
        if hasattr(ccdata, 'etenergies'):
            calcs.append('Excited States')
        # If our list is empty, assume we did an SP.
        if len(calcs) == 0 and ( hasattr(ccdata, 'scfenergies') or hasattr(ccdata, 'mpenergies') or hasattr(ccdata, 'ccenergies') ):
                calcs = ['Single Point']
        # Return our list.
        return calcs
    
    @classmethod
    def get_methods_from_cclib(self, ccdata):
        """
        Get a list of method types (DFT, MP etc.) from the data provided by cclib.
        
        :param ccdata: The data from cclib.
        :return: A list of strings of the methods used (each method appears once only; the order has no special significance). The list may be empty.
        """
        return self.sorted_methods((set(ccdata.metadata.get('methods', []))))
        #return list(dict.fromkeys(ccdata.metadata.get('methods', [])))
    
    @classmethod
    def get_orbital_spin_type_from_cclib(self, ccdata):
        """
        Determine the orbital type (restricted, unrestricted) from the data provided by cclib.
        
        :return ccdata: The data from cclib.
        :return: A string describing the orbital type ('restricted', 'unrestricted' or 'none').
        """
        #mo_len = len(ccdata.get('moenergies', []))
        mo_len = len(getattr(ccdata, 'moenergies', []))
        if mo_len == 1:
            return "restricted"
        elif mo_len == 2:
            return "unrestricted"
        else:
            return "none"
    
    @classmethod
    def from_parser(self, parser):
        """
        Construct a Metadata object from an output file parser.
        
        :param parser: Output data parser.
        :return: A populated Metadata object.
        """
        try:
            # Do some processing on time objects.
            if parser.data.metadata.get('walltime', None) is not None:
                duration = timedelta(seconds = parser.data.metadata['walltime'])
                
            elif parser.data.metadata.get('cputime', None) is not None and parser.data.metadata.get('numcpus', None) is not None:
                # We can estimate the wall time from the CPU time and number of cpus.
                duration = timedelta(seconds = parser.data.metadata['cputime'] / parser.data.metadata['numcpus'])
                
            else:
                duration = None
                
            if parser.data.metadata.get('date', None) is not None:
                date = datetime.fromtimestamp(parser.data.metadata['date'])
            else:
                date = None
            
            return self(
                name = parser.data.metadata.get('name', None),
                log_files = parser.data.metadata.get('log_files', None),
                auxiliary_files = parser.data.metadata.get('aux_files', None),
                date = date,
                duration = duration,
                package = parser.data.metadata.get('package', None),
                package_version = parser.data.metadata.get('package_version', None),
                calculations = self.get_calculation_types_from_cclib(parser.data),
                success = parser.data.metadata.get('success', None),
                methods = self.get_methods_from_cclib(parser.data),
                functional = parser.data.metadata.get('functional', None),
                basis_set = parser.data.metadata.get('basis_set', None),
                charge = getattr(parser.data, 'charge', None),
                multiplicity = getattr(parser.data, 'mult', None),
                optimisation_converged = getattr(parser.data, 'optdone', None),
                temperature = getattr(parser.data, 'temperature', None),
                pressure = getattr(parser.data, 'pressure', None),
                orbital_spin_type = self.get_orbital_spin_type_from_cclib(parser.data),
            )
        except AttributeError:
            # There is no metadata available, give up.
            raise Result_unavailable_error("Metadata", "no metadata is available")


class Merged_metadata(Metadata):
    """
    A modified metadata class for merged calculation results.
    """
    
    
    def __init__(self, num_calculations, *args, **kwargs):
        """
        :param num_calculations: The number of merged calculations this metadata represents.
        """
        self.num_calculations = num_calculations
        super().__init__(*args, log_files = None, auxiliary_files = None, **kwargs)
        
    @classmethod
    def from_metadatas(self, *multiple_metadatas):
        """
        Create a merged metadata object from multiple metadata objects.
        
        :param multiple_metadatas: A list of metadata objects to merge.
        :returns: A new Merged_metadata object.
        """
        # Our merged metadata.
        merged_metadata = self(num_calculations = len(multiple_metadatas))
        for attr in ("name", "package", "package_version", "functional", "basis_set"):
            setattr(merged_metadata, attr, self.merged_attr(attr, multiple_metadatas))
            
        # We take the latest of the two dates.
        merged_metadata.date = latest_datetime(*[metadata.date for metadata in multiple_metadatas if metadata.date is not None])
        # And the total duration.
        merged_metadata.duration = total_timedelta(*[metadata.duration for metadata in multiple_metadatas if metadata.duration is not None])
        
        # Merge methods and calculations (but keep unique only).
        merged_metadata.calculations = list(dict.fromkeys(itertools.chain(*(metadata.calculations for metadata in multiple_metadatas))))
        merged_metadata.methods = self.sorted_methods(set(itertools.chain(*(metadata.methods for metadata in multiple_metadatas))))
        
        # We are only successful if all calcs are successful.
        merged_metadata.success = all((metadata.success for metadata in multiple_metadatas))
        converged = [metadata.optimisation_converged for metadata in multiple_metadatas]
        merged_metadata.optimisation_converged = False if False in converged else True if True in converged else None
        
        # Only keep these if all the same.
        for attr in ("temperature", "pressure"):
            # Get a unique list of the attrs.
            attr_set = set(getattr(metadata, attr) for metadata in multiple_metadatas)
            
            # Only keep if we have exactly one entry.
            if len(attr_set) == 1:
                setattr(merged_metadata, attr, tuple(attr_set)[0])
                
        # We take the first orbital spin type charge and mult, as this is what our orbital list will actually be.
        merged_metadata.orbital_spin_type = multiple_metadatas[0].orbital_spin_type
        merged_metadata.charge = multiple_metadatas[0].charge
        merged_metadata.multiplicity = multiple_metadatas[0].multiplicity
        
        return merged_metadata
    