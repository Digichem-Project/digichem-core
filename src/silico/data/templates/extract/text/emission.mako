## -*- coding: utf-8 -*-
##
<%!
	from silico.result.excited_states import Energy_state
	from silico.exception import Result_unavailable_error
%>\
##
<%page args="relaxed_excited_state, result_name = ''"/>\
##
##
<%
	if relaxed_excited_state is None:
		raise Result_unavailable_error("relaxed excited state", "there is no emission energy of the requested type")
		
	title = relaxed_excited_state.transition_type.capitalize() + " " + relaxed_excited_state.state_symbol + " " + "Emission Energy"
%>\
##
<%include file="title.mako" args="title=title, result_name=result_name"/>
##
## Transition data.
Excited energy /eV: ${"{:0.2f}".format(relaxed_excited_state.excited_energy)}
Excited multiplicity: ${Energy_state.multiplicity_number_to_string(relaxed_excited_state.excited_multiplicity).capitalize()}
Ground energy /eV: ${"{:0.2f}".format(relaxed_excited_state.ground_energy)}
Ground multiplicity: ${Energy_state.multiplicity_number_to_string(relaxed_excited_state.ground_multiplicity).capitalize()}
Emission type: ${relaxed_excited_state.emission_type.capitalize()}
##
## Excited state data
${relaxed_excited_state.state_symbol} energy /eV: ${"{:0.2f}".format(relaxed_excited_state.energy)}
${relaxed_excited_state.state_symbol} wavelength /nm: ${"{:0.0f}".format(relaxed_excited_state.wavelength)}
${relaxed_excited_state.state_symbol} colour: ${relaxed_excited_state.color}
%if relaxed_excited_state.oscillator_strength is not None:
${relaxed_excited_state.state_symbol} oscillator strength: ${"{:0.2f}".format(relaxed_excited_state.oscillator_strength)}
%endif
