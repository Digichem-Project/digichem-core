## -*- coding: utf-8 -*-

<%page args="excited_states"/>

<%
	# Get our excited states grouped by multiplicity.
	grouped_states = excited_states.group()
	
	num_singlets = excited_states.num_singlets
	num_triplets = excited_states.num_triplets
	num_other = len(excited_states) - (num_singlets + num_triplets)
%>

<div class="resultsContainer resultsContainer--excitedStates">
	<div class="reportHeader reportHeader--minor reportHeader--results reportHeader--excitedStates">Excited States</div>
	<table class="results results--et">
		<%include file="/excited_states/dest_result_rows.mako" args="excited_states = excited_states" />
		%if num_singlets > 0:
		<tr>
			<td class="results__name">No. of singlets:</td>
			<td class="results__value">${num_singlets}</td>
		</tr>
		%endif
		%if num_triplets > 0:
		<tr>
			<td class="results__name">No. of triplets:</td>
			<td class="results__value">${num_triplets}</td>
		</tr>
		%endif
		%if num_other > 0:
		<tr>
			<td class="results__name">No. other states:</td>
			<td class="results__value">${num_other}</td>
		</tr>
		%endif
	</table>
</div>
