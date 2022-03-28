## -*- coding: utf-8 -*-

<%page args="excited_states"/>

<%!
	from silico.misc.text import text_float, andjoin
	from silico.result.excited_states import Energy_state
	import inflect
%>

<%
    inflector = inflect.engine()

    # Get our excited states grouped by multiplicity.
    mult_pairs, mults = report.result.excited_states.group_pairs()
%>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "excited states summary")}:</div>
		Summary of the calculated excited states.
	</div>
    <table class="resultsTable__table">
        %for mult, states in mults.items():
        	<%
        		state = states[0]
        	%>
        	<tr class="resultsTable__row">
        		<td class="resultsTable__title resultsTable__cell">
        			No. calculated
        			%if mult is not None and mult < 5:
        				${inflector.plural(Energy_state.multiplicity_number_to_string(mult))}
        			%else:
        				states with multiplicity ${mult}
        			%endif
        		</td>
        		<td class="resultsTable__value resultsTable__cell">${len(states)}</td>
        	</tr>
        	<tr class="resultsTable__row">
        		<td class="resultsTable__title resultsTable__cell">E<sub>${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub></sub></td>
        		<td class="resultsTable__value resultsTable__cell">${"{:.2f}".format(state.energy)} eV</td>
        	</tr>
        	<tr class="resultsTable__row">
        		<td class="resultsTable__title resultsTable__cell">${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub> wavelength (colour, CIE)</td>
        		<td class="resultsTable__value resultsTable__cell">
        			${"{:.0f}".format(state.wavelength)} nm (${state.color} <%include file="/excited_states/color.mako" args="colorRGB = state.rgb"/>, ${"({}, {})".format(*state.CIE_xy)})
        		</td>
        	</tr>
        	<tr class="resultsTable__row">
        		<td class="resultsTable__title resultsTable__cell">${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub> oscillator strength</td>
        		<td class="resultsTable__value resultsTable__cell">${text_float(state.oscillator_strength)}</td>
        	</tr>
        %endfor
		%for mult_pair in mult_pairs:
			<%
			state1 = mults[mult_pair[0]][0]
			state2 = mults[mult_pair[1]][0]
			%>
			<tr class="resultsTable__row">
	        		<td class="resultsTable__title resultsTable__cell">ΔE<sub>${state1.multiplicity_symbol}${state2.multiplicity_symbol}</sub></td>
	        		<td class="resultsTable__value resultsTable__cell">${text_float(state1.energy - state2.energy)} eV</td>
	        	</tr>
		%endfor
		<tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Simulated Absorption Peaks</td>
            <td class="resultsTable__value resultsTable__cell">
            %if 'simulated_absorption_graph' in report.images and report.images['simulated_absorption_graph'].safe_get_file() is not None:
            	<%
            		selected_peaks = report.images['simulated_absorption_graph'].selected_peaks(0, 5)
            	%>
            	${andjoin(selected_peaks)}
            	%if len(report.images['simulated_absorption_graph'].peaks) > len(selected_peaks):
            	...
            	%endif
            	nm
            %else:
            	N/A
            %endif
            </td>
        </tr>
    </table>
</div>
