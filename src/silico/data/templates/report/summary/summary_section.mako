## -*- coding: utf-8 -*-

<%!
	from silico.exception import Result_unavailable_error
%>

<%page args="result" />

<%
	# Get our transition_dipole_moment if we can.
	try:
		S1 = result.excited_states.get_state("S(1)")
		transition_dipole_moment = S1.transition_dipole_moment
	except Result_unavailable_error:
		# No S1 available.
		transition_dipole_moment = None
%>

<div class="section section--summary">
	<h2 class="section__header">Summary of Results</h2>
	<div class="section__body section__body--summary">
	%if result.metadata is not None:
		<%include file="/metadata/metadata_results.mako" args="metadata = result.metadata"/>
	%endif
	%if len(result.SCF_energies) > 0:
		<%include file="/energy/energy_results.mako" args="energies = result.SCF_energies"/>
	%endif
	%if len(result.MP_energies) > 0:
		<%include file="/energy/energy_results.mako" args="energies = result.MP_energies"/>
	%endif
	%if len(result.CC_energies) > 0:
		<%include file="/energy/energy_results.mako" args="energies = result.CC_energies"/>
	%endif
	%if len(result.alignment) > 0:
		<%include file="/geometry/geometry_results.mako" args="alignment = result.alignment"/>
	%endif
	%if len(result.molecular_orbitals) > 0:
		<%include file="/orbitals/HOMO_LUMO_results.mako" args="molecular_orbitals = result.molecular_orbitals"/>
	%endif
	%if len(result.beta_orbitals) > 0:
		<%include file="/orbitals/HOMO_LUMO_results.mako" args="molecular_orbitals = result.beta_orbitals"/>
	%endif
	%if result.dipole_moment is not None:
		<%include file="/dipole_moment/dipole_moment_results.mako" args="dipole_moment = result.dipole_moment"/>
	%endif
	%if transition_dipole_moment is not None:
		<%include file="/dipole_moment/dipole_moment_results.mako" args="dipole_moment = transition_dipole_moment"/>
	%endif
	%if len(result.vibrations) > 0:
		<%include file="/vibrations/vibrations_results.mako" args="vibrations = result.vibrations"/>
	%endif
	%if result.vertical_emission is not None:
		<%include file="/emission/emission_results.mako" args="relaxed_excited_state = result.vertical_emission"/><span></span>
	%endif
	%if result.adiabatic_emission is not None:
		<%include file="/emission/emission_results.mako" args="relaxed_excited_state = result.adiabatic_emission"/><span></span>
	%endif
	%if len(result.excited_states) > 0:
		<%include file="/excited_states/excited_states_results.mako" args="excited_states = result.excited_states"/><span></span>
	%endif
	</div>
</div>