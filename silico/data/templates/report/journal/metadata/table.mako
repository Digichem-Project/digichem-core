<%page args="report"/>

<%!
	from silico import misc
	from silico.result.excited_states import Energy_state
	import itertools
%>

<%
	# An easier list of metadatas to process.
	# If this is a merged result, report.result.metadatas will not include the merged metadata, this one will.
	if len(report.result.metadatas) == 1:
		metadatas = [report.result.metadata]
	else:
		metadatas = list(itertools.chain([report.result.metadata], report.result.metadatas)) 
%>

<div class="resultsTable resultsTable--complete">
	<div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "metadata")}:</div>
		Summary of overall calculation metadata. [a]: The date and time at which the calculation was completed. [b]: Total combined duration in real-time (wall-time) for all components of the calculation. [c]: Temperature used for thermochemistry analysis. [b]: Pressure used for thermochemistry analysis.
	</div>
	<table class="resultsTable__table">
		##############
	    ## Headers. ##
	    ##############
	    <tr class="resultTable__row">
	    	%if len(metadatas) > 1:
	    	<td class="resultsTable__title resultsTable__cell">Calculation no.</td>
	    	%endif
	    	<td class="resultsTable__title resultsTable__cell">Date<sup>[a]</sup></td>
	    	<td class="resultsTable__title resultsTable__cell">Duration<sup>[b]</sup></td>
			<td class="resultsTable__title resultsTable__cell">Success (Converged)</td>
	        <td class="resultsTable__title resultsTable__cell">Computational package</td>
		    <td class="resultsTable__title resultsTable__cell">Level of theory</div>
	        <td class="resultsTable__title resultsTable__cell">Calculations</td>
	        <td class="resultsTable__title resultsTable__cell">Wavefunction</td>
	        <td class="resultsTable__title resultsTable__cell">Multiplicity</td>
	        <td class="resultsTable__title resultsTable__cell">T<sup>[c]</sup> <span class="iblock">/ K</span></td>
	        <td class="resultsTable__title resultsTable__cell">P<sup>[d]</sup> <span class="iblock">/ atm</span></td>
	        ##<td class="resultsTable__title resultsTable__cell">No. merged calculations</td>
	    </tr>
	    
	    ###########
	    ## Data. ##
	    ###########
	    %for index, metadata in enumerate(metadatas):
	    <tr class="resultTable__row">
	    	%if len(metadatas) > 1:
	    	<td class="resultsTable__value resultsTable__cell">
	    		${index if index != 0 else "Combined"}
	    	</td>
	    	%endif
	    	<td class="resultsTable__value resultsTable__cell">
	    	%if metadata.date is not None:
	        	${misc.date_to_string(metadata.date)}
	        %else:
	        	N/A
	    	%endif
	    	</td>
	    	<td class="resultsTable__value resultsTable__cell">
	    	%if metadata.duration is not None:
		        ${misc.timedelta_to_string(metadata.duration)}
		    %else:
		    	N/A
		    %endif
	    	</td>
	    	<td class="resultsTable__value resultsTable__cell">
		    	%if metadata.success:
			    	<span class="resultsTable__value--good">True</span>
			    %elif not metadata.success:
			    	<span class="resultsTable__value--bad">False</span>
			    %else:
			    	<span class="resultsTable__value--bad">N/A</span>
			    %endif
			    %if metadata.optimisation_converged is None:
			    	<span class="">(N/A)</span>
			    %elif metadata.optimisation_converged:
			        (<span class="resultsTable__value--good">True</span>)
			    %else:
			        (<span class="resultsTable__value--bad">False</span>)
			    %endif
		    </td>
		    <td class="resultsTable__value resultsTable__cell">
		    %if metadata.package_string != "":
		    	${metadata.package_string}
		    %else:
		    	N/A
		    %endif
		    </td>
	    	<td class="resultsTable__value resultsTable__cell">${metadata.level_of_theory}</div>
	    	<td class="resultsTable__value resultsTable__cell">
		    %if len(metadata.calculations) > 0:
		        ${metadata.calculations_string}
		    %else:
		    	N/A
		    %endif
		    </td>
	        <td class="resultsTable__value resultsTable__cell">
		    %if metadata.orbital_spin_type is not None:
		        ${metadata.orbital_spin_type}
		    %else:
		    	N/A
		    %endif
		    </td>
	    	<td class="resultsTable__value resultsTable__cell">
	    	%if metadata.multiplicity is not None:
	    		${metadata.multiplicity} (${Energy_state.multiplicity_number_to_string(metadata.multiplicity)})
	    	%else:
	    		N/A
			%endif
			</td>
	    	<td class="resultsTable__value resultsTable__cell">
	    	%if metadata.temperature is not None:
		        ${metadata.temperature}
		    %else:
		    	N/A
		    %endif
	    	</td>
	    	<td class="resultsTable__value resultsTable__cell">
		    %if metadata.pressure is not None:
		        ${metadata.pressure}
		    %else:
		    	N/A
		    %endif
		    </td>
		    ##<td class="resultsTable__value resultsTable__cell">${metadata.num_calculations}</td>
	    </tr>
	    %endfor
	</table>
</div>