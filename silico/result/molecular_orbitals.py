# Classes for representing MOs.
from itertools import filterfalse
from silico.result import Result_container
from silico.result.base import Result_object
from silico.exception import Result_unavailable_error
from silico.image.vmd import Orbital_image_maker, Combined_orbital_image_maker
from silico.file.cube import Cube_maker
from pathlib import Path
from silico.image.orbitals import Orbital_diagram_maker

class Molecular_orbital_list(Result_container):
	"""
	Class for representing a group of MOs.
	"""
	
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
	@property
	def HOMO_energy(self):
		"""
		Get the energy of the highest occupied orbital in this list.
		
		Use get_orbital(HOMO_difference = 0) to retrieve the HOMO as an object.
		
		:raises Result_unavailable_error: If the HOMO is not available.
		:return The orbital energy (in eV).
		"""
		return self.get_orbital(HOMO_difference = 0).energy
	
	@property
	def LUMO_energy(self):
		"""
		Get the energy of the lowest unoccupied orbital in this list.
		
		Use get_orbital(HOMO_difference = 1) to retrieve the LUMO as an object.
		
		:raises Result_unavailable_error: If the LUMO is not available.
		:return The orbital energy (in eV).
		"""
		return self.get_orbital(HOMO_difference = 1).energy
		
	@property
	def HOMO_LUMO_energy(self):
		"""
		Get ΔE HOMO-LUMO; the energy difference between the HOMO and LUMO.
		
		Depending on how the orbitals are occupied, it might not make sense to calculate the HOMO-LUMO energy.
		
		:raises Result_unavailable_error: If the HOMO or LUMO is not available.
		:return: The HOMO-LUMO energy gap (in  eV).
		"""
		# Get out orbitals.
		HOMO = self.get_orbital(HOMO_difference = 0)
		LUMO = self.get_orbital(HOMO_difference = 1)
		# Return the difference.
		return float(LUMO) - float(HOMO)

	@property
	def spin_type(self):
		"""
		Get the spin type (alpha, beta etc.) of the orbitals in this list.
		
		:raises Result_unavailable_error: If there are no orbitals in this list.
		:return: The spin type, one of either 'alpha' or 'beta' for unrestricted calcs, 'none' for restricted calcs, or 'mixed' if multiple spin types are present in this list.
		"""
		# Unique spin types in our list.
		spin_types = list(set([orbital.spin_type for orbital in self]))
		
		if len(spin_types) == 1:
			return spin_types[0]
		elif len(spin_types) > 1:
			return "mixed"
		else:
			raise Result_unavailable_error("Orbital Spin", "There are no orbitals")
	
	def get_orbital(self, criteria = None, *, label = None, HOMO_difference = None, level = None):
		"""
		Retrieve an orbital based on some property.
		
		Only one of the criteria should be specified.
		
		:raises Result_unavailable_error: If the requested MO could not be found.
		:param criteria: A string describing the orbital to get. The meaning of criteria is determined automatically based on its content. If criteria begins with '+' or '-' then it is used as a HOMO_difference. Otherwise, if criteria is a valid integer, then it is used as a level. Otherwise, criteria is assumed to be an orbital label.
		:param label: The label of the orbital to get.
		:param HOMO_difference: The distance from the HOMO of the orbital to get.
		:param level: The level of the orbital to get.
		:return: The Molecular_orbital object.
		"""
		# Use the search() method to do our work for us.
		try:
			return self.search(criteria, label = label, HOMO_difference = HOMO_difference, level = level)[0]
		except IndexError:
			#raise ValueError("Unable to find orbital '{}'".format(label))
			raise Result_unavailable_error("Orbital", "could not find orbital with criteria = '{}', label = '{}', HOMO_difference = '{}', level = '{}'".format(criteria, label, HOMO_difference, level))
	
	def search(self, criteria = None, *, label = None, HOMO_difference = None, level = None):
		"""
		Attempt to retrieve a number of orbitals based on some property.
		
		Only one of the criteria should be specified.
		
		:param criteria: A string describing the orbitals to get. The meaning of criteria is determined automatically based on its content. If criteria begins with '+' or '-' then it is used as a HOMO_difference. Otherwise, if criteria is a valid integer, then it is used as a level. Otherwise, criteria is assumed to be an orbital label.
		:param label: The label of the orbitals to get.
		:param HOMO_difference: The distance from the HOMO of the orbitals to get.
		:param level: The level of the orbitals to get.
		:return: A (possibly empty) Molecular_orbital_list object.
		"""
		# If we've been given a generic criteria, decide what it is actually talking about.		
		if criteria is not None:
			try:				
				# If our string start with a sing (+ or -) then it's a HOMO_difference.
				if criteria[:1] == "+" or criteria[:1] == "-":
					HOMO_difference = int(criteria)
				# If its and integer (or a string that looks like one), then its a level.
				elif criteria.isdigit() or isinstance(criteria, int):
					level = int(criteria)
				# Otherwise we assume it a label.
				else:
					label = criteria
				
			except Exception:
				# We couldn't parse criteria, get upset.
				raise ValueError("Unable to understand given search criteria '{}'".format(criteria))
		
		# Now get our filter func.
		if label is not None:
			filter_func = lambda mo: mo.label != label
		elif HOMO_difference is not None:
			filter_func = lambda mo: mo.HOMO_difference != HOMO_difference
		elif level is not None:
			filter_func = lambda mo: mo.level != level
		else:
			raise ValueError("Missing criteria to search by; specify one of 'criteria', 'label', 'HOMO_difference' or 'level'")
		
		# Now search.
		return type(self)(filterfalse(filter_func, self))
	
	def ordered(self):
		"""
		Return a copy of this list of MOs that is ordered in terms of energy and removes duplicate MOs.
		"""
		ordered_list = type(self)(set(self))
		ordered_list.sort(key = lambda mo: mo.level)
		return ordered_list
		
	@classmethod
	def from_cclib(self, ccdata, cls = None):
		"""
		Construct a Molecular_orbital_list object from the data provided by cclib.
		
		:param ccdata: Result object as provided by cclib.
		:param cls: Optional class of objects to populate this list with, should inherit from Molecular_orbital. Defaults to Molecular_orbital if only one set of orbitals are available, or Alpha_orbital if both alpha and beta are available (in which case you should call Molecular_orbital_list.from_cclib() again with cls = Beta_orbital to get beta as well).
		:returns: The new Molecular_orbital_list object. The list will be empty if no MO data is available.
		"""
		try:
			# Set our default class if we've not been given one.
			if cls is None:
				# Check to see if we have only 'alpha' or beta as well.
				if len(ccdata.moenergies) == 1:
					cls = Molecular_orbital
				else:
					cls = Alpha_orbital
			# Get our list.
			return self(cls.list_from_cclib(ccdata))
		except AttributeError:
			return self()
		
	def set_file_options(self, output_dir, output_name, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		This method will also call set_file_options on the molecular orbitals contained in this list.
		
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		"""
		# First set options for our list (because we want to steal some of their cube files later).
		for mo in self:
			mo.set_file_options(output_dir, output_name, **kwargs)
		
		try:
			# Now get our spin type.
			if self.spin_type == "alpha":
				spin_type = "alpha_"
			elif self.spin_type == "beta":
				spin_type = "beta_"
			else:
				spin_type = ""
				
			# Save our orbital diagram.
			self._files['energy_diagram'] = Orbital_diagram_maker.from_image_options(
				Path(output_dir, "Orbital Diagram", output_name + ".{}orbitals.png".format(spin_type)),
				molecular_orbitals = self,
				**kwargs
			)
			# A version of the diagram with only the HOMO/LUMO
			self._files['HOMO_LUMO_energy_diagram'] = Orbital_diagram_maker.from_image_options(
				Path(output_dir, "Orbital Diagram", output_name + ".{}HOMO_LUMO.png".format(spin_type)),
				molecular_orbitals = type(self)(orbital for orbital in self if orbital.HOMO_difference == 0 or orbital.HOMO_difference == 1),
				**kwargs
			)
			
			# Also get our HOMO/LUMO combined image.
			# Get our FMOs.
			HOMO = self.get_orbital(HOMO_difference = 0)
			LUMO = self.get_orbital(HOMO_difference = 1)
			
			self._files['HOMO_LUMO_image'] = Combined_orbital_image_maker.from_image_options(
				Path(output_dir, "HOMO LUMO", output_name + ".{}HOMO_LUMO.jpg".format(spin_type)),
				HOMO_cube_file = HOMO.get_file('cube_file'),
				LUMO_cube_file = LUMO.get_file('cube_file'),
				**kwargs
			)
		except Result_unavailable_error:
			# We couldn't find our HOMO/LUMO.
			self._files.pop('HOMO_LUMO_image', None)
	
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		# No files from ourself to cleanup, but could be lots in our children.
		for orbital in self:
			orbital.cleanup_intermediate_files()
		
	@property
	def energy_diagram(self):
		return self.get_file('energy_diagram')
	
	@property
	def HOMO_LUMO_image(self):
		return self.get_file('HOMO_LUMO_image')
	
	def find_common_level(self, *other_lists, HOMO_difference):
		"""
		Find either:
			The orbital with the lowest level that has no less than the given negative HOMO_difference.
		or:
			The orbital with the highest level that has no more than the given positive HOMO_difference.
		Across one or more orbital lists.
		
		The method is useful for determining which orbitals to traverse between two limits from the HOMO/LUMO gap.
		
		:raises Result_unavailable_error: If all of the given orbital_lists (including this one) are empty.
		:param *other_lists: Optional lists to search. If none are given, then only this orbital_list is search.
		:param HOMO_difference: The distance from the HOMO to search for. Negative values indicate HOMO-n, positive values indicate LUMO+(n-1). The LUMO should be at +1 by definition.
		:return: The level (as an integer) of the matching orbital.
		"""
		# Our list of orbitals that match our criteria.
		found_orbitals = []
		
		# Our list of orbital_list objects to look through.
		orbital_lists = list(other_lists)
		orbital_lists.append(self)
		
		# Cant think of a better way to do this...
		if HOMO_difference <= 0:
			search_func = lambda orbital: orbital.HOMO_difference < HOMO_difference
		else:
			search_func = lambda orbital: orbital.HOMO_difference > HOMO_difference
		
		# Loop through each list and search.
		for orbital_list in orbital_lists:
			# Now search each list for orbitals that match our criteria.
			matching_orbitals = list(filterfalse(search_func, orbital_list))
			
			# We can just add all the orbitals we find because we'll only look at the lowest/highest anyway.
			found_orbitals.extend(matching_orbitals)
			
		# Get a list of orbital levels that match our criteria.
		orbital_levels = [orbital.level for orbital in found_orbitals]
		
		# Now either return the smallest or largest orbital level, depending on what we were asked for.
		try:
			if HOMO_difference <= 0:
				# The lowest orbital below HOMO.
				return min(orbital_levels)
			else:
				# The highest orbital above LUMO.
				return max(orbital_levels)
		except ValueError:
			# Min/Max couldn't find anything, this should only happen if all orbital_lists are completely empty.
			raise Result_unavailable_error("Common orbital level", "there are no orbitals")
	

class Molecular_orbital(Result_object):
	"""
	Class representing a molecular orbital.
	"""
	
	# True MOs don't have a spin.
	spin_type = "none"
	
	def __init__(self, level, HOMO_difference, symmetry, energy):
		"""
		Constructor for MOs.
		
		:param level: The 'level' of this MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
		:param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
		:param symmetry: The symmetry of this MO.
		:param energy: The energy of this MO (in eV).
		"""
		super().__init__()
		self.level = level
		self.HOMO_difference = HOMO_difference
		self.symmetry = symmetry
		self.energy = energy
	
	def __float__(self):
		return self.energy
	
	@property
	def HOMO_level(self):
		"""
		The level of the HOMO in the collection of orbitals of which this orbital is a member.
		"""
		return self.level - self.HOMO_difference
	
	@property
	def LUMO_level(self):
		"""
		The level of the LUMO in the collection of orbitals of which this orbital is a member.
		"""
		return self.HOMO_level +1		
	
	@property
	def label(self):
		"""
		A label describing this MO in terms of its proximity to the HOMO and LUMO.
		
		:return: A string label, of the form either HOMO-n or LUMO+n.
		"""
		# The label we return depends on how close to the HOMO we are.
		if self.level == self.HOMO_level:
			# We are the HOMO.
			label = "HOMO"
		elif self.level < self.HOMO_level:
			# We are below the HOMO (and presumably occupied).
			label = "HOMO{}".format(self.level - self.HOMO_level)
		elif self.level == self.LUMO_level:
			# We are the LUMO.
			label = "LUMO"
		else:
			# We are above the LUMO (and presumably unoccupied).
			label = "LUMO{0:+}".format(self.level - self.LUMO_level)
			
		return label
	
	def __eq__(self, other):
		"""
		Equality operator between MOs.
		"""
		return self.label == other.label
	
	def __hash__(self):
		"""
		Hash operator.
		"""
		return hash(tuple(self.label))
		
	
	def set_file_options(self, output_dir, output_name, *, cube_file = None, **kwargs):
		"""
		Set the options that will be used to create images from this object.
		
		:param output_dir: A pathlib Path object to the directory within which our files should be created.
		:param output_name: A string that will be used as the start of the file name of the files we create.
		:param output_base: The base directory where all output will be written to.
		:param fchk_file: An optional fchk_file to use to render the MO image. If 'cube_file' is not given, this must be given.
		:param cube_file: An optional cube file to use to render the MO image. If None, an appropriate cube file will be created automatically.
		:param options: A silico Config dictionary (or a dictionary with the same structure at least) of options to set. This should match the format laid out in the silico config file.
		"""
		# If we haven't been given a cube file (which is pretty likely), then create one from the fchk file we should have been given.
		if cube_file is None:
			# Decide what type of orbital we need.
			if self.spin_type == "alpha":
				cubegen_type = "AMO"
			elif self.spin_type == "beta":
				cubegen_type = "BMO"
			else:
				cubegen_type = "MO"
			
			# Get our cube maker object.
			cube_file = Cube_maker.from_image_options(
				Path(output_dir, self.label, output_name + ".{}.cube".format(self.label)),
				cubegen_type = cubegen_type,
				orbital = self.level,
				**kwargs)
		
		# Save our cube file.
		self._files['cube_file'] = cube_file
		
		# And then save our orbital image.
		self._files['orbital_image'] = Orbital_image_maker.from_image_options(
			Path(output_dir, self.label, output_name + ".{}.jpg".format(self.label)),
			cube_file = cube_file,
			**kwargs)
		
	def cleanup_intermediate_files(self):
		"""
		Remove any intermediate files that may have been created by this object.
		"""
		# Remove our cube file.
		super().cleanup_intermediate_files('cube_file')
	
	@property
	def orbital_image(self):
		"""
		Get an Orbital_image_maker object that can create an orbital density image of this orbital.
		"""
		return self.get_file("orbital_image")
	
	@property
	def cube_file(self):
		"""
		Get a Cube_maker object that can create a cube file of this orbital.
		"""
		return self.get_file("cube_file")

	@classmethod
	def from_cclib(self, index, symmetry, energy, HOMO_index):
		"""
		Create a Molecular_orbital object from the data provided by cclib.
		
		:param index: The index of this MO.
		:param symmetry: The symmetry of this MO.
		:param energy: The energy of this MO (in eV).
		:param HOMO_index: The index of the HOMO in this MO set.
		"""
		return self(index +1, index - HOMO_index, symmetry, energy)
	
	# The index used to access data from cclib (which always has two lists, one for alpha one for beta).
	ccdata_index = 0
	
	@classmethod
	def list_from_cclib(self, ccdata):
		"""
		Create a list of Molecular_orbital objects from the data provided by cclib.
		
		:param ccdata: Result object as provided by cclib.
		:return: A list of Molecular_orbital objects. The list will be empty if no MO is available.
		"""
		try:
			return [self.from_cclib(index, symmetry, energy, ccdata.homos[self.ccdata_index]) for index, (symmetry, energy) in enumerate(zip(ccdata.mosyms[self.ccdata_index], ccdata.moenergies[self.ccdata_index]))]
		except (AttributeError, IndexError):
			return []
	
class Unrestricted_orbital(Molecular_orbital):
	"""
	Top-level class for unrestricted orbitals.
	"""
	
	def __init__(self, level, HOMO_difference, symmetry, energy, spin_type):
		"""
		Constructor for MOs.
		
		:param level: The 'level' of the MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
		:param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
		:param symmetry: The symmetry of the MO.
		:param energy: The energy of the MO (in eV).
		:param spin_type: The spin of this spin-orbital (either alpha or beta).
		"""
		# Call parent first.
		super().__init__(level, HOMO_difference, symmetry, energy)
		self.spin_type = spin_type
	
	@property
	def label(self):
		# Get the base of the label first.
		label = super().label
		# Append our spin type.
		label = "{} ({})".format(label, self.spin_type)
		return label

class Alpha_orbital(Unrestricted_orbital):
	"""
	An alpha spin orbital (these types of orbitals are only singly occupied, electrons are spin-up).
	"""
	
	def __init__(self, level, HOMO_difference, symmetry, energy):
		"""
		Constructor for alpha MOs.
		
		:param level: The 'level' of the MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
		:param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
		:param symmetry: The symmetry of the MO.
		:param energy: The energy of the MO (in eV).
		"""
		super().__init__(level, HOMO_difference, symmetry, energy, "alpha")
		
class Beta_orbital(Unrestricted_orbital):
	"""
	A beta spin orbital (these types of orbitals are only singly occupied, electrons are spin-down).
	"""
	
	# Beta orbitals use the other list in cclib.
	ccdata_index = 1
	
	def __init__(self, level, HOMO_difference, symmetry, energy):
		"""
		Constructor for beta MOs.
		
		:param level: The 'level' of the MO (essentially an index), where the lowest MO has a level of 1, increasing by 1 for each higher orbital.
		:param HOMO_difference: The distance of this MO from the HOMO. A negative value means this orbital is HOMO-n. A positive value means this orbital is HOMO+n (or LUMO+(n-1). A value of 0 means this orbital is the HOMO. A value of +1 means this orbital is the LUMO.
		:param symmetry: The symmetry of the MO.
		:param energy: The energy of the MO (in eV).
		"""
		super().__init__(level, HOMO_difference, symmetry, energy, "beta")
	