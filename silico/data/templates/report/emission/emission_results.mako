## -*- coding: utf-8 -*-

<%!
    from silico.result.excited_states import Energy_state
%>

<%
	emission_rate = emission.safe_get('emission_rate')
%>

<%page args="emission"/>

<div class="resultsContainer resultsContainer--excitedStates">
    <div class="reportHeader reportHeader--minor reportHeader--results reportHeader--excitedStates">${emission.transition_type.capitalize()} ${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> Emission</div>
    <table class="results results--et re">
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
        <tr>
            <td class="results__name">Emission type:</td>
            <td class="results__value">${emission.emission_type.capitalize()}</td>
        </tr>
        <%include file="/excited_states/state_result_rows.mako" args="state=emission"/>
        %if emission_rate is not None:
        <tr>
			<td class="results__name">${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> rate /s<sup>-1</sup>:</td>
			<td class="results__value">${"{:.2e}".format(emission_rate)}</td>
		</tr>
		%endif
    </table>
</div>
