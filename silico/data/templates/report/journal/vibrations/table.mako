## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	from silico.misc.text import text_float
%>

<div class="resultsTable">
	<div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "vibrations")}:</div> Energies of the calculated vibrational frequencies.
		
	</div>
	<table class="resultsTable__table">
		##############
	    ## Headers. ##
	    ##############
	    <tr class="resultTable__row">
	    	<td class="resultsTable__title resultsTable__cell">Number</td>
	    	<td class="resultsTable__title resultsTable__cell">Symmetry</td>
	    	<td class="resultsTable__title resultsTable__cell">Frequency /cm<sup>-1</sup></td>
	    	<td class="resultsTable__title resultsTable__cell">Intensity /km mol<sup>-1</sup></td>
	    </tr>
	    
	    ###########
	    ## Data. ##
	    ###########
	    %for vibration in report.result.vibrations:
	    <tr class="resultTable__row">
	    	<td class="resultsTable__value resultsTable__cell">${vibration.level}</td>
	    	<td class="resultsTable__value resultsTable__cell">${vibration.symmetry}</td>
	    	<td class="resultsTable__value resultsTable__cell">${"{:.4f}".format(vibration.frequency)}</td>
	    	<td class="resultsTable__value resultsTable__cell">${text_float(vibration.intensity, 4)}</td>
	    </tr>
	    %endfor
	</table>
</div>