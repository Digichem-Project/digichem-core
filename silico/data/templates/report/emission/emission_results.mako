## -*- coding: utf-8 -*-

<%!
    from silico.result.excited_states import Energy_state
%>

<%page args="relaxed_excited_state"/>

<div class="resultsContainer resultsContainer--excitedStates">
    ##<div class="reportHeader reportHeader--minor reportHeader--results">Geometrically Relaxed Excited State</div>
    <div class="reportHeader reportHeader--minor reportHeader--results reportHeader--excitedStates">${relaxed_excited_state.transition_type.capitalize()} Emission Energy</div>
    <table class="results results--et re">
        ##<tr>
        ##    <td class="results__name">Transition type:</td>
        ##    <td class="results__value">${relaxed_excited_state.transition_type.capitalize()}</td>
        ##</tr>
        <tr>
            <td class="results__name">Excited energy:</td>
            <td class="results__value">${"{:0.2f}".format(relaxed_excited_state.excited_energy)} eV</td>
        </tr>
        <tr>
            <td class="results__name">Excited multiplicity:</td>
            <td class="results__value">${Energy_state.multiplicity_number_to_string(relaxed_excited_state.excited_multiplicity).capitalize()}</td>
        </tr>
        <tr>
            <td class="results__name">Ground energy:</td>
            <td class="results__value">${"{:0.2f}".format(relaxed_excited_state.ground_energy)} eV</td>
        </tr><tr>
            <td class="results__name">Ground multiplicity:</td>
            <td class="results__value">${Energy_state.multiplicity_number_to_string(relaxed_excited_state.ground_multiplicity).capitalize()}</td>
        </tr>
##         <tr>
##             <td class="results__name">Transition energy:</td>
##             <td class="results__value">${"{:0.2f}".format(relaxed_excited_state.energy)} eV</td>
##         </tr>
##         <tr>
##             <td class="results__name">Emission wavelength:</td>
##             <td class="results__value">${"{:0.0f}".format(relaxed_excited_state.wavelength)} nm</td>
##         </tr>
##         <tr>
##             <td class="results__name">Emission colour:</td>
##             <td class="results__value">${relaxed_excited_state.color}</td>
##         </tr>
        <tr>
            <td class="results__name">Emission type:</td>
            <td class="results__value">${relaxed_excited_state.emission_type.capitalize()}</td>
        </tr>
        <%include file="/excited_states/state_result_rows.mako" args="state=relaxed_excited_state"/>
    </table>
</div>
