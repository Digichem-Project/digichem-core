"""Code for metadata regarding finished calculation results."""

# General imports.
from datetime import timedelta, datetime
import math
import itertools
from deepmerge import conservative_merger
from pathlib import Path
import copy

# Silico imports.
from silico.exception import Result_unavailable_error
from silico.result import Result_object
from silico.misc.time import latest_datetime, total_timedelta, date_to_string,\
    timedelta_to_string
from silico.misc.text import andjoin
import silico


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
            user = None,
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
            orbital_spin_type = None,
            silico_version = None):
        """
        Constructor for result Metadata objects.
        
        :param name: Optional name of this calculation result.
        :param user: The username of the user who parsed this result.
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
        self.user = user
        self.log_files = log_files if log_files is not None else []
        self.auxiliary_files = auxiliary_files if auxiliary_files is not None and len(auxiliary_files) != 0 else {}
        self.date = date
        self.duration = duration
        self.package = package
        self.package_version = package_version
        self.calculations = calculations if calculations is not None else []
        self.success = success
        self.methods = methods if methods is not None else []
        self.functional = functional
        self.basis_set = basis_set
        # TODO: charge and mult should be deprecated here, they are available in ground_state.
        self.charge = charge
        self.multiplicity = multiplicity
        self.optimisation_converged = optimisation_converged
        self.temperature = temperature
        self.pressure = pressure
        self.orbital_spin_type = orbital_spin_type
        # TOOD: Ideally this would be parsed from the calculation output somehow, but this is fine for now.
        self.silico_version = silico.version if silico_version is None else silico_version
    
    # TODO: This is more than a bit clumsy and in general the handling of names should be improved.
    @property
    def molecule_name(self):
        return Path(self.name).name if self.name is not None else None
    
    @property
    def level_of_theory(self):
        """
        A short-hand summary of the methods and basis sets used.
        """
        theories = []
        if len(self.converted_methods) > 0:
            #theories.extend(self.converted_methods)
            theories.append(self.converted_methods[-1])
            
        if self.basis_set is not None:
            theories.append(self.basis_set)
            
        return("/".join(theories))
        
    @classmethod
    def merge(self, *multiple_metadatas):
        """
        Merge multiple metadata objects into a single metadata.
        """
        return Merged_metadata.from_metadatas(*multiple_metadatas)
    
    @property
    def converted_methods(self):
        """
        Similar to the methods attribute but where DFT is replaced with the actual functional used.
        """
        return [self.functional if method == "DFT" and self.functional is not None else method for method in self.methods if method is not None]
    
    @property
    def description(self):
        desc = []
        if self.name is not None:
            desc.append(self.name)
        
        desc.append(self.identity_string)
        
        return ", ".join(desc)
        
    @property
    def identity(self):
        """
        A dictionary of critical attributes that identify a calculation.
        """
        # Get our list of methods, replacing 'DFT' with the functional used if available.
        return {
            "package": self.package,
            "calculations": self.calculations_string,
            "methods": ", ".join(self.converted_methods) if len(self.converted_methods) > 0 else None,
            "basis": self.basis_set,
        }
        
    @property
    def identity_string(self):
        """
        A string that identifies this calculation.
        """
        return " ".join([value for value in self.identity.values() if value is not None])
        
    
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
    
    @property
    def human_calculations_string(self):
        """
        A version of calculations_string that is more pleasant to read.
        """
        calculations = []
        if "Single Point" in self.calculations:
            calculations.append("single point energy")
        
        if "Optimisation" in self.calculations:
            calculations.append("optimised structure")
            
        if "Frequencies" in self.calculations:
            calculations.append("vibrational frequencies")
            
        if "Excited States" in self.calculations:
            calculations.append("excited states")
            
        if len(calculations) == 0:
            # Use a generic term.
            calculations.append("properties")
        
        return andjoin(calculations)
            
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
            if parser.data.metadata.get('wall_time', None) is not None:
                duration = sum(parser.data.metadata['wall_time'], timedelta())
                
            elif parser.data.metadata.get('cpu_time', None) is not None and parser.data.metadata.get('num_cpus', None) is not None:
                # We can estimate the wall time from the CPU time and number of cpus.
                duration = sum(parser.data.metadata['cpu_time'], timedelta()) / parser.data.metadata['num_cpus']
                
            else:
                duration = None
                
            if parser.data.metadata.get('date', None) is not None:
                date = datetime.fromtimestamp(parser.data.metadata['date'])
            else:
                date = None
            
            # TODO: This doesn't seem to make sense; the parser already contains a metadata object...
            return self(
                name = parser.data.metadata.get('name', None),
                user = parser.data.metadata.get('user', None),
                log_files = parser.data.metadata.get('log_files', None),
                auxiliary_files = parser.data.metadata.get('auxiliary_files', None),
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
        
    def dump(self, silico_options):
        """
        Get a representation of this result object in primitive format.
        """
        # Most attributes we can just dump as is.
        attr_dict = {
            "name": self.name,
            "log_files": [str(log_file) for log_file in self.log_files],
            "auxiliary_files": {aux_file_name: str(aux_file) for aux_file_name, aux_file in self.auxiliary_files.items()}
        }
        
        attrs = [
            "charge",
            "multiplicity",
            "user",
            "package",
            "package_version",
            "silico_version",
            "calculations",
            "methods",
            "functional",
            "basis_set",
            "orbital_spin_type",
            "success",
            "optimisation_converged",
        ]
        attr_dict.update({attr: getattr(self, attr) for attr in attrs})
        
        # Add some more complex stuff.
        attr_dict['date'] = {
            "value": self.date.timestamp() if self.date is not None else None,
            "units": "s",
            "string": date_to_string(self.date) if self.date is not None else None
        }
        attr_dict['duration'] = {
            "value": self.duration.total_seconds() if self.duration is not None else None,
            "units": "s",
            "string": timedelta_to_string(self.duration) if self.duration is not None else None
        }
        attr_dict["temperature"] = {
            "value": self.temperature,
            "units": "K"
        }
        attr_dict["pressure"] = {
            "value": self.pressure,
            "units": "atm"
        }
        
        return attr_dict
    
    @classmethod
    def from_dump(self, data, result_set):
        """
        Get an instance of this class from its dumped representation.
        
        :param data: The data to parse.
        :param result_set: The partially constructed result set which is being populated.
        """
        # Assemble our args to pass to the constructor.
        # Most of these can be used from the metadata dict as is.
        kwargs = copy.deepcopy(data)
        
        # For more complex fields, use the data item.
        for attr in ['date', 'duration', 'temperature', "pressure"]:
            kwargs[attr] = data[attr]['value']
        
        kwargs['date'] = datetime.fromtimestamp(kwargs['date']) if kwargs['date'] is not None else None
        kwargs['duration'] = timedelta(seconds = kwargs['duration'])  if kwargs['duration'] is not None else None
        
        return self(**kwargs)


class Merged_metadata(Metadata):
    """
    A modified metadata class for merged calculation results.
    """
    
    
    def __init__(self, num_calculations, *args, **kwargs):
        """
        :param num_calculations: The number of merged calculations this metadata represents.
        """
        super().__init__(*args, log_files = None, auxiliary_files = None, **kwargs)
        self.num_calculations = num_calculations
        
    @classmethod
    def from_metadatas(self, *multiple_metadatas):
        """
        Create a merged metadata object from multiple metadata objects.
        
        :param multiple_metadatas: A list of metadata objects to merge.
        :returns: A new Merged_metadata object.
        """
        # Our merged metadata.
        merged_metadata = self(num_calculations = len(multiple_metadatas))
        for attr in ("name", "user", "package", "package_version", "functional", "basis_set"):
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
        
        # CAREFULLY merge aux files, so that later files do not overwrite earlier ones.
        # This is useful behaviour because it matches how other results are merged, so certain aux files
        # (turbomole density files) will still match their respective results.
        for metadata in multiple_metadatas:
            merged_metadata.auxiliary_files = conservative_merger.merge(merged_metadata.auxiliary_files, metadata.auxiliary_files)
        
        return merged_metadata
    