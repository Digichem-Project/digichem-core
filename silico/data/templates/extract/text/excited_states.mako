## -*- coding: utf-8 -*-
##
<%!
    from silico.result.excited_states import Energy_state
    from silico.exception import Result_unavailable_error
%>\
##
<%page args="excited_states, result_name = ''"/>\
##
##
<%
    if len(excited_states) == 0:
        raise Result_unavailable_error("excited states", "there are no excited states")
    
    # A dictionary where each key is a multiplicity.
    grouped_states = excited_states.group()
%>\
##
<%include file="title.mako" args="title='Excited States', result_name=result_name"/>
##
## dEst (if available).
%if excited_states.safe_get('singlet_triplet_energy') is not None:
ΔEst /eV: ${"{:0.2f}".format(excited_states.singlet_triplet_energy)}
%endif
##
## Now print something interesting about lowest state of each mult.
##
%for group in grouped_states:
No. ${Energy_state.multiplicity_number_to_string(group)}s: ${len(grouped_states[group])}
##
<%
    state = grouped_states[group][0]
%>\
##
${state.state_symbol} energy /eV: ${"{:0.2f}".format(state.energy)}
${state.state_symbol} wavelength /nm: ${"{:0.0f}".format(state.wavelength)}
${state.state_symbol} colour: ${state.color}
${state.state_symbol} CIE X: ${"{:0.2f}".format(state.CIE_xy[0])}
${state.state_symbol} CIE Y: ${"{:0.2f}".format(state.CIE_xy[1])}
${state.state_symbol} oscillator strength: ${"{:0.2f}".format(state.oscillator_strength)}
%endfor
