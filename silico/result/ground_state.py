# Contains classes for representing the ground state of our system.
# This file is closely related to the excited_states.py file.
# These definitions may change or move.

from silico.result.excited_states import Energy_state
from silico.exception.base import Result_unavailable_error

class Ground_state(Energy_state):
	"""
	Class for representing the ground state of a system/molecule.
	"""
	
	def __init__(self, charge, multiplicity, energy):
		"""
		Constructor for ground state objects.
		
		:param charge: The charge of the ground states as a signed integer.
		:param multiplicity: The multiplicity as an integer of float.
		:param energy: The absolute energy of the GS in eV.
		"""
		# The 'level' refers to the ordering of energy states, the ground state is always the lowest so have level 0.
		# Similarly, the multiplicity_level refers to the ordering of energy states that share the same multiplicity. Regardless of what multiplicity our GS has, we will be the lowest by definition.
		# Energy is our energy above the GS, so that will always be 0.
		super().__init__(level = 0, multiplicity = multiplicity, multiplicity_level = 0, energy = 0)
		self.charge = charge
		# Also save our absolute energy.
		self.absolute_energy = energy
		
	@classmethod
	def from_results(self, results):
		"""
		Construct a Ground_state object from a set of results.
		
		As different energy types can be present in a single calculation, the following (arbitrary?) precedence is used:
		1) Coupled-cluster energy.
		2) The highest level Moller-Plesset energy.
		3) SCF energy.
		
		:param results: A result_set object.
		:return: The populated Ground_state object. Note this object is NOT added to the given Result_set object.
		"""
		if len(results.CC_energies) > 0:
			energy = results.CC_energies.final
		elif len(results.MP_energies) > 0:
			energy = results.MP_energies.final
		elif len(results.SCF_energies) > 0:
			energy = results.SCF_energies.final
		else:
			# We have no energies that we can understand. Give up.
			raise Result_unavailable_error("Ground_state", "No ground state energies that we understand are available")
		# Return our constructed object.
		return self(results.metadata.system_charge, results.metadata.system_multiplicity, energy)
		
		