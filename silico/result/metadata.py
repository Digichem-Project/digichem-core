# Code for processing results from calculations and generating reports.

# General imports.
from datetime import timedelta, datetime

# Silico imports.
from silico.exception import Result_unavailable_error
from silico.result import Result_object

class Metadata(Result_object):
    """
    Class for storing calculation metadata.
    """

    def __init__(
            self,
            name = None,
            date = None,
            duration = None,
            package = None,
            package_version = None,
            calculations = None,
            calc_success = None,
            calc_methods = None,
            calc_functional = None,
            calc_basis_set = None,
            system_charge = None,
            system_multiplicity = None,
            optimisation_converged = None,
            calc_temperature = None,
            calc_pressure = None,
            orbital_spin_type = None):
        """
        Constructor for result Metadata objects.
        
        :param name: Optional name of this calculation result.
        :param date: Optional date (datetime object) of this calculation result.
        :param duration: Optional duration (timedelta object) of this calculation.
        :param package: Optional string identifying the computational chem program that performed the calculation.
        :param package: Optional string identifying the version of the computational chem program that performed the calculation.
        :param calculations: A list of strings of the different calculations carried out (Opt, Freq, TD, TDA, SP etc).
        :param calc_success: Whether the calculation completed successfully or not.
        :param calc_methods: List of methods (DFT, HF, MPn etc) used in the calculation.
        :param calc_functional: Functional used in the calculation.
        :param calc_basis_set: Basis set used in the calculation.
        :param system_charge: Charge (positive or negative integer) of the system studied.
        :param system_multiplicity: The multiplicity of the system.
        :param optimisation_converged: Whether the optimisation converged or not.
        :param calc_temperature: The temperature used in the calculation (not always relevant).
        :param calc_pressure: The pressure used in the calculation (not always relevant).
        :param orbital_spin_type: The types of orbitals that have been calculated, either 'restricted' or 'unrestricted'.
        """
        self.name = name
        self.date = date
        self.duration = duration
        self.package = package
        self.package_version = package_version
        self.calculations = calculations if calculations is not None else []
        self.calc_success = calc_success
        self.calc_methods = calc_methods if calc_methods is not None else []
        self.calc_functional = calc_functional
        self.calc_basis_set = calc_basis_set
        self.system_charge = system_charge
        self.system_multiplicity = system_multiplicity
        self.optimisation_converged = optimisation_converged
        self.calc_temperature = calc_temperature
        self.calc_pressure = calc_pressure
        self.orbital_spin_type = orbital_spin_type
    
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
    def calc_methods_string(self):
        """
        Get the list of calculation methods as a single string, or None if there are no methods set.
        """
        return ", ".join(self.calc_methods) if len(self.calc_methods) != 0 else None
            
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
        return list(set(ccdata.metadata.get('methods', [])))
    
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
    def from_parser(self, parser, name = None, date = None, duration = None):
        """
        Construct a Metadata object from an output file parser.
        
        :param name: Optional name of this calculation result.
        :param date: Optional date (datetime object) of this calculation result.
        :param duration: Optional duration (timedelta object) of this calculation.
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
                date = date,
                duration = duration,
                package = parser.data.metadata.get('package', None),
                package_version = parser.data.metadata.get('package_version', None),
                calculations = self.get_calculation_types_from_cclib(parser.data),
                calc_success = parser.data.metadata.get('success', None),
                calc_methods = self.get_methods_from_cclib(parser.data),
                calc_functional = parser.data.metadata.get('functional', None),
                calc_basis_set = parser.data.metadata.get('basis_set', None),
                system_charge = getattr(parser.data, 'charge', None),
                system_multiplicity = getattr(parser.data, 'mult', None),
                optimisation_converged = getattr(parser.data, 'optdone', None),
                calc_temperature = getattr(parser.data, 'temperature', None),
                calc_pressure = getattr(parser.data, 'pressure', None),
                orbital_spin_type = self.get_orbital_spin_type_from_cclib(parser.data),
            )
        except AttributeError:
            # There is no metadata available, give up.
            raise Result_unavailable_error("Metadata", "no metadata is available")
                        