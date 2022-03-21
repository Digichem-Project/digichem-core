## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	from silico.misc.text import text_float
%>

<div class="resultsTable">
	<div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "excited_states")}:</div> Energies and other properties of the calculated excited states.
		
	</div>
	<table class="resultsTable__table">
		##############
	    ## Headers. ##
	    ##############
	    <tr class="resultTable__row">
	    	<td class="resultsTable__title resultsTable__cell">Number</td>
	    	<td class="resultsTable__title resultsTable__cell">Symbol</td>
	    	<td class="resultsTable__title resultsTable__cell">Symmetry</td>
	    	<td class="resultsTable__title resultsTable__cell">Energy /eV</td>
	    	<td class="resultsTable__title resultsTable__cell">Wavelength /nm</td>
	    	<td class="resultsTable__title resultsTable__cell">Colour (CIE x,y)</td>
	    	<td class="resultsTable__title resultsTable__cell">Oscillator Strength</td>
	    	<td class="resultsTable__title resultsTable__cell">Transitions (Probability)</td>
	    </tr>
	    
	    ###########
	    ## Data. ##
	    ###########
	    %for excited_state in report.result.excited_states:
	    <tr class="resultTable__row">
	    	<td class="resultsTable__value resultsTable__cell">${excited_state.level}</td>
	    	<td class="resultsTable__value resultsTable__cell">${excited_state.multiplicity_symbol}<sub>${excited_state.multiplicity_level}</sub></td>
	    	<td class="resultsTable__value resultsTable__cell">${excited_state.symmetry}</td>
	    	<td class="resultsTable__value resultsTable__cell">${"{:.4f}".format(excited_state.energy)}</td>
	    	<td class="resultsTable__value resultsTable__cell">${"{:.2f}".format(excited_state.wavelength)}</td>
	    	<td class="resultsTable__value resultsTable__cell">${excited_state.color} <%include file="/excited_states/color.mako" args="colorRGB=excited_state.rgb"/> (${"{:.2f}, {:.2f}".format(excited_state.CIE_xy[0], excited_state.CIE_xy[1])})</td>
	    	<td class="resultsTable__value resultsTable__cell">${text_float(excited_state.oscillator_strength, 4)}</td>
	    	<td class="resultsTable__value resultsTable__cell resultsTable__cell--excitedStateTransitions">
	    		%for transition in excited_state.transitions:
                <div class="resultsTable__value resultsTable__value--excitedStateTransition">
                    ${transition.starting_mo.label} → ${transition.ending_mo.label} (${"{:0.2f}".format(transition.probability)})
                </div>
                %endfor
	    	</td>
	    </tr>
	    %endfor
	</table>
</div>