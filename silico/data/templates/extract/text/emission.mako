## -*- coding: utf-8 -*-
##
<%!
    from silico.result.excited_states import Energy_state
    from silico.exception import Result_unavailable_error
%>\
##
<%page args="emission, result_name = ''"/>\
##
##
<%
    if emission is None:
        raise Result_unavailable_error("relaxed excited state", "there is no emission energy of the requested type")
        
    title = emission.transition_type.capitalize() + " " + emission.state_symbol + " " + "Emission Energy"
    
    emission_rate = emission.safe_get('emission_rate')
%>\
##
<%include file="title.mako" args="title=title, result_name=result_name"/>
##
## Transition data.
Excited energy /eV: ${"{:0.2f}".format(emission.excited_energy)}
Excited multiplicity: ${Energy_state.multiplicity_number_to_string(emission.excited_multiplicity).capitalize()}
Ground energy /eV: ${"{:0.2f}".format(emission.ground_energy)}
Ground multiplicity: ${Energy_state.multiplicity_number_to_string(emission.ground_multiplicity).capitalize()}
Emission type: ${emission.emission_type.capitalize()}
##
## Excited state data
${emission.state_symbol} energy /eV: ${"{:0.2f}".format(emission.energy)}
${emission.state_symbol} wavelength /nm: ${"{:0.0f}".format(emission.wavelength)}
${emission.state_symbol} colour: ${emission.color}
%if emission.oscillator_strength is not None:
${emission.state_symbol} oscillator strength: ${"{:0.2f}".format(emission.oscillator_strength)}
%endif
%if emission_rate is not None:
${emission.state_symbol} rate /s-1: ${"{:.2e}".format(emission_rate)}
%endif
