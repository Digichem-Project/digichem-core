## -*- coding: utf-8 -*-

<%!
    from silico.result.excited_states import Energy_state
%>

<%
	emission_rate = emission.safe_get('emission_rate')
%>

<%page args="emission"/>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "{} {} summary".format(emission.state_symbol, emission.transition_type))}:</div>
		Summary of the ${emission.transition_type} emission from the ${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> state.
	</div>
    <table class="resultsTable__table">
        <tr>
            <td class="resultsTable__title resultsTable__cell">Excited energy</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(emission.excited_energy)} eV</td>
        </tr>
        <tr>
            <td class="resultsTable__title resultsTable__cell">Excited multiplicity</td>
            <td class="resultsTable__value resultsTable__cell">${Energy_state.multiplicity_number_to_string(emission.excited_multiplicity).capitalize()}</td>
        </tr>
        <tr>
            <td class="resultsTable__title resultsTable__cell">Ground energy</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(emission.ground_energy)} eV</td>
        </tr><tr>
            <td class="resultsTable__title resultsTable__cell">Ground multiplicity</td>
            <td class="resultsTable__value resultsTable__cell">${Energy_state.multiplicity_number_to_string(emission.ground_multiplicity).capitalize()}</td>
        </tr>
        <tr>
            <td class="resultsTable__title resultsTable__cell">Emission type</td>
            <td class="resultsTable__value resultsTable__cell">${emission.emission_type.capitalize()}</td>
        </tr>
        <tr>
		    <td class="resultsTable__title resultsTable__cell">${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> energy</td>
		    <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(emission.energy)} eV</td>
		</tr>
		<tr class="resultsTable__row">
        		<td class="resultsTable__title resultsTable__cell">${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> wavelength (colour, CIE)</td>
        		<td class="resultsTable__value resultsTable__cell">
        			${"{:.0f}".format(emission.wavelength)} nm (${emission.color} <%include file="/excited_states/color.mako" args="colorRGB = emission.rgb"/>, ${"({:.2f}, {:.2f})".format(*emission.CIE_xy)})
        		</td>
        	</tr>
		<tr>
		    <td class="resultsTable__title resultsTable__cell">${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> oscillator strength</td>
		    <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(emission.oscillator_strength) if emission.safe_get('oscillator_strength') is not None else "N/A"}</td>
		</tr>
        <tr>
			<td class="resultsTable__title resultsTable__cell">${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> rate</td>
			%if emission_rate is not None:
				<td class="resultsTable__value resultsTable__cell">${"{:.2e}".format(emission_rate)} /s<sup>-1</sup></td>
			%else:
				<td class="resultsTable__value resultsTable__cell">N/A</td>
			%endif
		</tr>
    </table>
</div>
