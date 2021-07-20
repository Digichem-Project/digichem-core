## -*- coding: utf-8 -*-

<%!
    from silico.result.excited_states import Energy_state
%>

<%page args="emission"/>

<div class="resultsContainer resultsContainer--excitedStates">
    ##<div class="reportHeader reportHeader--minor reportHeader--results">Geometrically Relaxed Excited State</div>
    <div class="reportHeader reportHeader--minor reportHeader--results reportHeader--excitedStates">${emission.transition_type.capitalize()} Emission Energy</div>
    <table class="results results--et re">
        ##<tr>
        ##    <td class="results__name">Transition type:</td>
        ##    <td class="results__value">${emission.transition_type.capitalize()}</td>
        ##</tr>
        <tr>
            <td class="results__name">Excited energy:</td>
            <td class="results__value">${"{:0.2f}".format(emission.excited_energy)} eV</td>
        </tr>
        <tr>
            <td class="results__name">Excited multiplicity:</td>
            <td class="results__value">${Energy_state.multiplicity_number_to_string(emission.excited_multiplicity).capitalize()}</td>
        </tr>
        <tr>
            <td class="results__name">Ground energy:</td>
            <td class="results__value">${"{:0.2f}".format(emission.ground_energy)} eV</td>
        </tr><tr>
            <td class="results__name">Ground multiplicity:</td>
            <td class="results__value">${Energy_state.multiplicity_number_to_string(emission.ground_multiplicity).capitalize()}</td>
        </tr>
##         <tr>
##             <td class="results__name">Transition energy:</td>
##             <td class="results__value">${"{:0.2f}".format(emission.energy)} eV</td>
##         </tr>
##         <tr>
##             <td class="results__name">Emission wavelength:</td>
##             <td class="results__value">${"{:0.0f}".format(emission.wavelength)} nm</td>
##         </tr>
##         <tr>
##             <td class="results__name">Emission colour:</td>
##             <td class="results__value">${emission.color}</td>
##         </tr>
        <tr>
            <td class="results__name">Emission type:</td>
            <td class="results__value">${emission.emission_type.capitalize()}</td>
        </tr>
        <%include file="/excited_states/state_result_rows.mako" args="state=emission"/>
    </table>
</div>
