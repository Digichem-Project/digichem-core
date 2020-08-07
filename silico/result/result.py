# Code for processing results from calculations and generating reports.

# Imports.
from silico.result.alignment.AAA import Adjusted_average_angle
from silico.result.alignment.AA import Average_angle
from silico.result.alignment.FAP import Furthest_atom_pair
from silico.result.alignment import Minimal
from silico.exception import Result_unavailable_error
from silico.result.energy import CC_energy_list, MP_energy_list,\
	SCF_energy_list
from silico.result.atoms import Atom_list
from silico.result.dipole_moment import Dipole_moment
from silico.result.molecular_orbitals import Molecular_orbital_list,\
	Beta_orbital
from silico.result.excited_states import Excited_state_list
from silico.result.vibrations import Vibration_list
from silico.result.base import Result_object
from silico.result.transition_dipole_moment import Transition_dipole_moment
from silico.result.alignment.base import Alignment
from logging import getLogger
from pathlib import Path
import cclib.io
from silico.exception.base import Silico_exception
#import silico
import silico.result.excited_states
from silico.file.cube import Spin_cube_maker
from silico.image.vmd import Spin_density_image_maker
from silico.result.emission import Relaxed_excited_state

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
	def title(self):
		"""
		A string Title describing this result.
		"""
		title = ", ".join(self.calculations)
		if True:
			title += " ({})".format(silico.result.excited_states.Energy_state.multiplicity_number_to_string(self.system_multiplicity).capitalize())
		return title
	
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
	def from_cclib(self, ccdata, *args, **kwargs):
		"""
		Construct a Metadata object from the data provided by cclib.
		
		:param name: Optional name of this calculation result.
		:param date: Optional date (datetime object) of this calculation result.
		:param duration: Optional duration (timedelta object) of this calculation.
		:param ccdata: The data provided by cclib.
		:return: A populated Metadata object.
		"""
		try:
			return self(
				*args,
				package = ccdata.metadata.get('package', None),
				package_version = ccdata.metadata.get('package_version', None),
				calculations = self.get_calculation_types_from_cclib(ccdata),
				calc_success = ccdata.metadata.get('success', None),
				calc_methods = self.get_methods_from_cclib(ccdata),
				calc_functional = ccdata.metadata.get('functional', None),
				calc_basis_set = ccdata.metadata.get('basis_set', None),
				system_charge = getattr(ccdata, 'charge', None),
				system_multiplicity = getattr(ccdata, 'mult', None),
				optimisation_converged = getattr(ccdata, 'optdone', None),
				calc_temperature = getattr(ccdata, 'temperature', None),
				calc_pressure = getattr(ccdata, 'pressure', None),
				orbital_spin_type = self.get_orbital_spin_type_from_cclib(ccdata),
				**kwargs
			)
		except AttributeError:
			# There is no metadata available, give up.
			raise Result_unavailable_error("Metadata", "no metadata is available")
		
		
class Result_set(Result_object):
	"""
	Class that represents a collection of results from a calculation.
	
	This class is a bit heavy and might not be around for long...
	"""
	
	def __init__(
			self,
			gaussian_log_file =  None,
			metadata = None,
			ccenergies = None,
			mpenergies = None,
			scfenergies = None,
			atoms = None,
			alignment = None,
			dipole_moment = None,
			molecular_orbitals = None,
			beta_orbitals = None,
			excited_states = None,
			vertical_emission = None,
			adiabatic_emission = None,
			vibrations = None):
		"""
		Constructor for Result_set objects.
		
		:param gaussian_log_file: The Gaussian log file from which these results were read.
		:param metadata: Optional Metadata result object.
		:param ccenergies: Optional Energy_list object of coupled-cluster energies.
		:param mpenergies: Optional Energy_list object of Moller-Plesset energies.
		:param scfenergies: Optional Energy_list object of self-consistent field energies (SCF is the type of energy printed for normal HF and DFT).
		:param atoms: Optional Atom_list object of atom positions.
		:param dipole_moment: Optional dipole_moment object.
		:param molecular_orbitals: Optional Molecular_orbital_list object.
		:param beta_orbitals: Optional Beta MOs. If this is not None, then molecular_orbitals is assumed to refer to the Alpha MOs.
		:param excited_states: Optional Excited_state_list object.
		:param vertical_emission: A Relaxed_excited_state object representing the vertical emission energy.
		:param adiabatic_emission: A Relaxed_excited_state object representing the adiabatic emission energy.
		:param vibrations: Optional molecular Vibrations object.
		"""
		super().__init__()
		self.gaussian_log_file = Path(gaussian_log_file) if gaussian_log_file is not None else None
		self.metadata = metadata
		self.CC_energies = ccenergies
		self.MP_energies = mpenergies
		self.SCF_energies = scfenergies
		self.dipole_moment = dipole_moment
		self.atoms = atoms
		self.alignment = alignment
		self.molecular_orbitals = molecular_orbitals
		self.beta_orbitals = beta_orbitals
		self.excited_states = excited_states
		self.vibrations = vibrations
		self.vertical_emission = vertical_emission
		self.adiabatic_emission = adiabatic_emission
		
	@property
	def energy(self):
		"""
		The total energy of this calculation.
		
		This convenience property is the energy at the highest level of theory available (CC > MP > SCF).
		
		:raises Result_unavailable_error: If no total energy is available.
		"""
		# Try CC first.
		if len(self.CC_energies) > 0:
			return self.CC_energies.final
		elif len(self.MP_energies) > 0:
			return self.MP_energies.final
		else:
			return self.SCF_energies.final
		
	def set_file_options(self, output_dir, output_name, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		:param output_base: The base directory where all output will be written to.
		:param fchk_file: An optional fchk_file to use to render the MO image. If 'cube_file' is not given, this must be given.
		:param options: A silico Config dictionary (or a dictionary with the same structure at least) of options to set. This should match the format laid out in the silico config file.
		"""
		# Get our cube maker object.
		cube_file = Spin_cube_maker.from_image_options(
			Path(output_dir, "Spin Density", output_name + ".spin.cube"),
			cubegen_type = "Spin",
			orbital = "SCF",
			**kwargs)
		
		# Save our cube file.
		self._files['cube_file'] = cube_file
		
		# And then save our orbital images.
		self._files['spin_image_positive'] = Spin_density_image_maker.from_image_options(
			Path(output_dir, "Spin Density", output_name + ".spin_pos.jpg"),
			cube_file = cube_file,
			spin = "positive",
			**kwargs)
		self._files['spin_image_negative'] = Spin_density_image_maker.from_image_options(
			Path(output_dir, "Spin Density", output_name + ".spin_neg.jpg"),
			cube_file = cube_file,
			spin = "negative",
			**kwargs)
		self._files['spin_image_both'] = Spin_density_image_maker.from_image_options(
			Path(output_dir, "Spin Density", output_name + ".spin_both.jpg"),
			cube_file = cube_file,
			spin = "both",
			**kwargs)
			
	@property
	def spin_image_positive(self):
		return self.get_file('spin_image_positive')
	
	@property
	def spin_image_negative(self):
		return self.get_file('spin_image_negative')
	
	@property
	def spin_image_both(self):
		return self.get_file('spin_image_both')
		
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		# Remove our cube file.
		super().cleanup_intermediate_files('cube_file')	
	
		
	@property
	def transition_dipole_moment(self):
		"""
		The S1 dipole moment, commonly referred to as THE transition dipole moment (although this name is ambiguous).
		
		None is returned if the S1 dipole moment is not available.		
		"""
		try:
			S1 = self.excited_states.get_state("S(1)")
			return S1.transition_dipole_moment
		except Result_unavailable_error:
			# No S1 available.
			return None

		
	@classmethod
	def from_cclib(
			self,
			ccdata,
			*args,
			gaussian_log_file = None,
			name = None,
			date = None,
			duration = None,
			alignment_class_name = None,
			vertical_emission_ground_result = None,
			adiabatic_emission_ground_result = None,
			emission_excited_result = None,
			emission_excited_state = None,
			**kwargs):
		"""
		Construct a Result_set object from the data provided by cclib.
		
		Missing results will be populated as empty lists (if a list-like object) or as None if singletons.
		
		:param ccdata: Result object as provided by cclib.
		:param name: Optional name of this calculation result.
		:param date: Optional date (datetime object) of this calculation result.
		:param duration: Optional duration (timedelta object) of this calculation.
		:param alignment_class_name: The name/handle of an alignment class to use for molecule alignment. Defaults to Kebab+ (the average angle method).
		:param vertical_emission_ground_result: An optional additional Result_set object which will be used as the ground state for calculation of vertical emission energy.
		:param adiabatic_emission_ground_result: An optional additional Result_set object which will be used as the ground state for calculation of adiabatic emision energy.
		:param emission_excited_result: An optional additional Result_set object which will be used as the excited state for calculation of vertical and/or adiabatic emission energy. If emission_ground_result is give, then emission_excited_result is not optional.
		:param emission_excited_state: Optionally either an Excited_state object or a string describing one ('S(1)' etc) for use in calculating vertical and/or adiabatic emission energy.
		:return: A populated Result_set object.
		"""
		alignment_class_name = "MIN" if alignment_class_name is None else alignment_class_name
		
		#try:
		# First get our list of MOs (because we need them for excited states too.
		molecular_orbitals = Molecular_orbital_list.from_cclib(ccdata)
		beta_orbitals = Molecular_orbital_list.from_cclib(ccdata, cls = Beta_orbital)
		
		# Metadata.
		metadata = Metadata.from_cclib(ccdata, name = name, date = date, duration = duration)
		
		# Our alignment orientation data.
		alignment =	Alignment.from_class_handle(alignment_class_name).from_cclib(ccdata, charge = metadata.system_charge)
		atoms = Atom_list.from_cclib(ccdata, charge = metadata.system_charge)
		
		# Get transition dipoles.
		if gaussian_log_file is not None:
			transition_dipole_moments = Transition_dipole_moment.list_from_log(gaussian_log_file, atoms = alignment)
		else:
			transition_dipole_moments = []
		
		
		# Get our object.
		result = self(
			*args,
			gaussian_log_file = gaussian_log_file,
			metadata = 		metadata,
			ccenergies = 	CC_energy_list.from_cclib(ccdata),
			mpenergies = 	MP_energy_list.from_cclib(ccdata),
			scfenergies = 	SCF_energy_list.from_cclib(ccdata),
			atoms = 		atoms,
			alignment = 	alignment,
			dipole_moment = Dipole_moment.from_cclib(ccdata, atoms = alignment),
			molecular_orbitals = 	molecular_orbitals,
			beta_orbitals = 		beta_orbitals,
			excited_states = 		Excited_state_list.from_cclib(ccdata, molecular_orbitals, beta_orbitals, transition_dipole_moments),
			vibrations = 	Vibration_list.from_cclib(ccdata),
			**kwargs
			)
		
		# Set emission energy if we have it.
	
		# Set our emission energies.
		# For vertical, we can also use the adiabatic excited energy
		try:
			self._set_emission(result, 'vertical', vertical_emission_ground_result, emission_excited_result, emission_excited_state)
		except Exception:
			getLogger(silico.logger_name).warning("Could not load vertical emission energy", exc_info = True)
		try:
			self._set_emission(result, 'adiabatic', adiabatic_emission_ground_result, emission_excited_result, emission_excited_state)
		except Exception:
			getLogger(silico.logger_name).warning("Could not load adiabatic emission energy", exc_info = True)

		# All done.
		return result
		#except Exception:
		#	raise Silico_exception("Unable to parse calculation '{}'".format(name))
		
	@classmethod
	def _set_emission(self, result, transition_type, emission_ground_result = None, emission_excited_result = None, emission_excited_state = None):
		"""
		"""
		if transition_type != "vertical" and transition_type != "adiabatic":
			raise ValueError("Unknown transition_type '{}'".format(transition_type))
		
		# Set our emission energy.
		if emission_excited_result is not None:
			setattr(result, transition_type + "_emission", Relaxed_excited_state.from_results(
				result,
				ground_state_result = emission_ground_result,
				excited_state_result = emission_excited_result,
				transition_type = transition_type,
				excited_state = emission_excited_state))
		elif emission_ground_result is not None:
			raise Silico_exception("Cannot calculate emission energy; no excited state given for ground state '{}'".format(emission_ground_result.metadata.name))
		
	@classmethod
	def from_calculation_file(self,
			calc_file_path,
			*args,
			gaussian_log_file = None,
			name = None,
			alignment_class_name = None,
			**kwargs):
		"""
		Construct a Result_set object from a calculation file.
		
		:param calc_file_path: Path to a calculation file to load.
		"""
		# Load emission results.
		for emission in ['vertical_emission_ground_result', 'adiabatic_emission_ground_result', 'emission_excited_result']:
			if kwargs.get(emission, None) is not None and not isinstance(kwargs.get(emission, None), self):
				# This emission 'result' is not a result (assume it is a path); try and load it.
				try:
					kwargs[emission] = Result_set.from_calculation_file(kwargs[emission], alignment_class_name = alignment_class_name)
				except Exception:
					raise Silico_exception("Error loading emission result file '{}'".format(kwargs[emission]))
		
		try:		
			# Output a message (because this is slow).
			getLogger(silico.logger_name).info("Reading result file '{}'".format(calc_file_path))
			
			# Read and parse the file with cclib.
			with open(calc_file_path, "rt") as calc_file:
				ccdata = cclib.io.ccread(calc_file)
				
			# Now call our next 'constructor'.
			return self.from_cclib(
				ccdata,
				*args,
				name = name if name is not None else str(calc_file_path),
				alignment_class_name = alignment_class_name,
				gaussian_log_file = calc_file_path if gaussian_log_file is None else gaussian_log_file,
				**kwargs)
		except Exception:
			raise Silico_exception("Unable to load calculation result file '{}'".format(calc_file_path))
		