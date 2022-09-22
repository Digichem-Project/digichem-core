## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	from silico.misc.text import text_float
%>

<div class="resultsTable">
	<div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "SOC")}:</div>
		Calculated SOC values between singlet and triplet states.
		[a]: SOC between the singlet state and triplet sub-state with quantum number +1.
		[b]: The same with the triplet sub-state with quantum number 0.
		[c]: The same with the triplet sub-state with quantum number +1.
		[d]: Root sum square of the SOC between the singlet state and all three triplet sub-states.
		[e]: The first order mixing coefficient (H<sub>SO</sub>/ΔE<sub>ST</sub>) between the singlet and triplet state.
		
	</div>
	<table class="resultsTable__table">
		##############
	    ## Headers. ##
	    ##############
	    <tr class="resultTable__row">
	    	<td class="resultsTable__title resultsTable__cell">Singlet State</td>
	    	<td class="resultsTable__title resultsTable__cell">Triplet State</td>
	    	<td class="resultsTable__title resultsTable__cell">H<sub>SO</sub> +1<sup>[a]</sup><br>/cm<sup>-1</sup></td>
	    	<td class="resultsTable__title resultsTable__cell">H<sub>SO</sub> 0<sup>[b]</sup><br>/cm<sup>-1</sup></td>
	    	<td class="resultsTable__title resultsTable__cell">H<sub>SO</sub> -1<sup>[c]</sup><br>/cm<sup>-1</sup></td>
	    	<td class="resultsTable__title resultsTable__cell">H<sub>SO</sub> Root Sum Square<sup>[d]</sup><br>/cm<sup>-1</sup></td>
	    	<td class="resultsTable__title resultsTable__cell">H<sub>SO</sub> Root Sum Square<sup>[d]</sup><br>/eV</td>
	    	<td class="resultsTable__title resultsTable__cell">λ<sup>[e]</sup></td>
	    </tr>
	    
	    ###########
	    ## Data. ##
	    ###########
	    %for soc in report.result.soc:
	    <tr class="resultTable__row">
	    	
	    	<td class="resultsTable__value resultsTable__cell">${soc.singlet_state.multiplicity_symbol}<sub>${soc.singlet_state.multiplicity_level}</sub></td>
            <td class="resultsTable__value resultsTable__cell">${soc.triplet_state.multiplicity_symbol}<sub>${soc.triplet_state.multiplicity_level}</sub></td>
            <td class="resultsTable__value resultsTable__cell">${text_float(soc.positive_one, 4)}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(soc.zero, 4)}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(soc.negative_one, 4)}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(soc.wavenumbers, 4)}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(soc.energy, 4)}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(soc.mixing_coefficient, 4)}</td>   
	    </tr>
	    %endfor
	</table>
</div>