## -*- coding: utf-8 -*-

<%page args="report"/>

<div class="resultsTable">
	<div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "atoms")}:</div>
		Coordinates of the atoms of the system under study, as aligned to the cartesian axes by the ${report.result.alignment.human_method_type} method.
	</div>
	<table class="resultsTable__table">
		##############
	    ## Headers. ##
	    ##############
	    <tr class="resultTable__row">
	    	<td class="resultsTable__title resultsTable__cell">Element</td>
	    	<td class="resultsTable__title resultsTable__cell">X Coord /Å</td>
	    	<td class="resultsTable__title resultsTable__cell">Y Coord /Å</td>
	    	<td class="resultsTable__title resultsTable__cell">Z Coord /Å</td>
	    </tr>
	    
	    ###########
	    ## Data. ##
	    ###########
	    %for atom in report.result.alignment:
	    <tr class="resultTable__row">
	    	<td class="resultsTable__value resultsTable__cell">${atom.element.symbol}</td>
	    	<td class="resultsTable__value resultsTable__cell">${"{:0.7f}".format(atom.coords[0])}</td>
	    	<td class="resultsTable__value resultsTable__cell">${"{:0.7f}".format(atom.coords[1])}</td>
	    	<td class="resultsTable__value resultsTable__cell">${"{:0.7f}".format(atom.coords[2])}</td>
	    </tr>
	    %endfor
	</table>
</div>