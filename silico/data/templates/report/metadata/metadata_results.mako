## -*- coding: utf-8 -*-

<%!
	from silico.result.excited_states import Energy_state
	from silico import misc
%>

<%page args="metadata" />


<div class="resultsContainer">
	<div class="reportHeader reportHeader--minor reportHeader--results">Metadata</div>
	<table class="results">
		%if metadata.date is not None:
		<tr>
			<td class="results__name">Date:</td>
			<td class="results__value">${misc.date_to_string(metadata.date)}</td>
		</tr>
		%endif
		%if metadata.duration is not None:
		<tr>
			<td class="results__name">Duration:</td>
			<td class="results__value">${misc.timedelta_to_string(metadata.duration)}</td>
		</tr>
		%endif
		<tr>
			<td class="results__name">Success:</td>
			%if metadata.calc_success:
			<td class="results__value results__value--good">True</td>
			%elif not metadata.calc_success:
			<td class="results__value results__value--bad">False</td>
			%else:
			<td class="results__value results__value--bad">Unknown</td>
			%endif
		</tr>
		%if metadata.optimisation_converged is not None:
		<tr>
			<td class="results__name">Converged:</td>
			%if metadata.optimisation_converged:
			<td class="results__value results__value--good">True</td>
			%else:
			<td class="results__value results__value--bad">False</td>
			%endif
		</tr>
		%endif
		%if metadata.package_string != "":
		<tr>
			<td class="results__name">Computational package:</td>
			<td class="results__value">${metadata.package_string}</td>
		</tr>
		%endif
		%if len(metadata.calc_methods) > 0:
		<tr>
			<td class="results__name">Methods:</td>
			<td class="results__value">${", ".join(metadata.calc_methods)}</td>
		</tr>
		%endif
		%if metadata.calc_functional is not None:
		<tr>
			<td class="results__name">Functional:</td>
			<td class="results__value">${metadata.calc_functional}</td>
		</tr>
		%endif
		%if metadata.calc_basis_set is not None:
		<tr>
			<td class="results__name">Basis set:</td>
			<td class="results__value">${metadata.calc_basis_set}</td>
		</tr>
		%endif
		%if metadata.orbital_spin_type is not None:
		<tr>
			<td class="results__name">Orbital spin:</td>
			<td class="results__value">${metadata.orbital_spin_type}</td>
		</tr>
		%endif
		%if metadata.system_multiplicity is not None:
		<tr>
			<td class="results__name">Multiplicity:</td>
			<td class="results__value">${metadata.system_multiplicity} (${Energy_state.multiplicity_number_to_string(metadata.system_multiplicity)})</td>
		</tr>
		%endif
		%if metadata.calc_temperature is not None:
		<tr>
			<td class="results__name">Calc temperature:</td>
			<td class="results__value">${metadata.calc_temperature} K</td>
		</tr>
		%endif
		%if metadata.calc_pressure is not None:
		<tr>
			<td class="results__name">Calc pressure:</td>
			<td class="results__value">${metadata.calc_pressure} atm</td>
		</tr>
		%endif
	</table>
</div>