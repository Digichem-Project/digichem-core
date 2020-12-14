from mako.lookup import TemplateLookup
from silico.extract import Result_extractor
from silico.extract import Result_extractor_group
import silico.extract
#from silico.extract.summary import *

class Text_summary_group_extractor(Result_extractor_group):
	"""
	Class for extracting particular sections from a number of Result_set objects in text format.
	"""
	
	@classmethod
	def get_default_extractors(self, **kwargs):
		"""
		Get a list of default extractor objects that can be used to convert a result set to dict format.
		"""
		return [
			Metadata_text_summary_extractor(**kwargs),
			Vibrations_text_summary_extractor(**kwargs),
			Geometry_text_summary_extractor(**kwargs),
			SCF_text_summary_extractor(**kwargs),
			MP_text_summary_extractor(**kwargs),
			CC_text_summary_extractor(**kwargs),
			Orbitals_text_summary_extractor(**kwargs),
			Beta_text_summary_extractor(**kwargs),
			PDM_text_summary_extractor(**kwargs),
			TDM_text_summary_extractor(**kwargs),
			Excited_states_text_summary_extractor(**kwargs),
			SOC_text_summary_extractor(**kwargs)
		]
	
	def join_results(self, extracted_results):
		"""
		Method called to combine a list of extracted results from multiple result sets.
		"""
		return "\n------------------------------------------------\n\n".join([result for result in extracted_results if result != ""])
	
	def extract(self, results):
		return super().extract(results, display_name = len(results) > 1)
	
	
	@classmethod
	def recursive_subclasses(self):
		"""
		Recursively get all the subclasses of this class.
		
		Note that text extractors don't actually extend from this group extractor.
		"""
		return Text_summary_extractor.recursive_subclasses()

class Text_summary_extractor(Result_extractor):
	"""
	Top-level class for all Text_summary_extractor classes.
	"""
	
	# The template file this extractor should use (relative to the '/extract/text/' template directory). Inheriting classes should modify this.
	TEMPLATE = None
	# A dictionary of keyword arguments that should be passed to the template. Inheriting classes should modify this.
	TEMPLATE_ARGS = {}
	
	def _get_template_args(self, result):
		"""
		Get the arguments (as a dict) that are to be passed to this extractor's template.

		This default implementation returns the value of self.TEMPLATE_ARGS, assuming that each value in TEMPLATE_ARGS is the name of an attribute in a Result_set.
		If your class needs more complex behaviour, write your own implementation.
		"""
		return {key:getattr(result, self.TEMPLATE_ARGS[key]) for key in self.TEMPLATE_ARGS}
			
	def _extract(self, result, display_name = False, **template_args):
		"""
		Convert a result set contained into text.
		
		:param display_name: Whether or not to include the file name.
		:param **template_args: Optional keyword arguments to pass to the template. If 0 are given, then self.TEMPLATE_ARGS are used.
		:return: The extracted results as a string.
		"""
		# Get our template.
		return TemplateLookup(directories = str(silico.default_template_directory())).get_template("/extract/text/{}".format(self.TEMPLATE)).render_unicode(
			result_name = result.metadata.name if display_name else "",
			**(self._get_template_args(result) if len(template_args) == 0 else template_args)
		)

		
class PDM_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for the PDM.
	"""
	CLASS_HANDLE = silico.extract.PDM_CLASS_HANDLE
	TEMPLATE = "dipole.mako"
	TEMPLATE_ARGS ={'dipole_moment': 'dipole_moment'}
	
class TDM_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for the TDM.
	"""
	CLASS_HANDLE = silico.extract.TDM_CLASS_HANDLE
	ALLOW_CRITERIA = True
	TEMPLATE = "dipole.mako"
	TEMPLATE_ARGS ={'dipole_moment': 'transition_dipole_moment'}
	
	def _extract_with_criteria(self, excited_state_id, *, result, **kwargs):
		"""
		Extract a specific TDM.
		
		:param excited_state_id: A string or int describing the excited state to get the dipole of. See Excited_state.get_state().
		"""
		# We've been asked for a specific dipole.
		dipole = result.excited_states.get_state(excited_state_id).transition_dipole_moment
		return super()._extract(result, dipole_moment = dipole, **kwargs)

	def _extract(self, result, *args, **kwargs):
		"""
		Extract the standard TDM.
		"""
		return super()._extract(result, *args, **kwargs)
	
class SCF_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for the SCF energy.
	"""
	CLASS_HANDLE = silico.extract.SCF_CLASS_HANDLE
	TEMPLATE = "energy.mako"
	TEMPLATE_ARGS ={'energy': 'SCF_energies'}
	
class MP_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for the MP (Moller-Plesset) energy.
	"""
	CLASS_HANDLE = silico.extract.MP_CLASS_HANDLE
	TEMPLATE = "energy.mako"
	TEMPLATE_ARGS ={'energy': 'MP_energies'}
	
class CC_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for the CC (coupled cluster) energy.
	"""
	CLASS_HANDLE = silico.extract.CC_CLASS_HANDLE
	TEMPLATE = "energy.mako"
	TEMPLATE_ARGS ={'energy': 'CC_energies'}
	
class Excited_states_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for excited states.
	"""
	CLASS_HANDLE = silico.extract.EXCITED_STATE_CLASS_HANDLE
	TEMPLATE = "excited_states.mako"
	TEMPLATE_ARGS ={'excited_states': 'excited_states'}
	
class SOC_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for spin orbit coupling.
	"""
	CLASS_HANDLE = silico.extract.SOC_CLASS_HANDLE
	TEMPLATE = "spin_orbit_coupling.mako"
	TEMPLATE_ARGS ={'SOC_list': 'spin_orbit_coupling'}
		
class Geometry_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for geometry/atom data.
	"""
	CLASS_HANDLE = silico.extract.GEOM_CLASS_HANDLE
	TEMPLATE = "geometry.mako"
	TEMPLATE_ARGS ={'alignment': 'alignment'}
	
class Metadata_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for metadata.
	"""
	CLASS_HANDLE = silico.extract.METADATA_CLASS_HANDLE
	TEMPLATE = "metadata.mako"
	TEMPLATE_ARGS ={'metadata': 'metadata'}
	
class Orbitals_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for (alpha) orbitals.
	"""
	CLASS_HANDLE = silico.extract.ORBITALS_CLASS_HANDLE
	TEMPLATE = "orbitals.mako"
	TEMPLATE_ARGS ={'orbitals': 'molecular_orbitals'}
	
class Beta_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for (beta) orbitals.
	"""
	CLASS_HANDLE = silico.extract.BETA_CLASS_HANDLE
	TEMPLATE = "orbitals.mako"
	TEMPLATE_ARGS ={'orbitals': 'beta_orbitals'}
	
class Vibrations_text_summary_extractor(Text_summary_extractor):
	"""
	Text extractor for vibrational frequencies..
	"""
	CLASS_HANDLE = silico.extract.VIBRATIONS_CLASS_HANDLE
	TEMPLATE = "vibrations.mako"
	TEMPLATE_ARGS ={'vibrations': 'vibrations'}
	
	