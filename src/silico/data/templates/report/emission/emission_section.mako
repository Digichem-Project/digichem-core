## -*- coding: utf-8 -*-

<%page args="relaxed_excited_state"/>

<div class="section">
	<h2 class="section__header">${relaxed_excited_state.transition_type.capitalize()} ${"{}<sub>{}</sub>".format(relaxed_excited_state.multiplicity_symbol, relaxed_excited_state.multiplicity_level)} Emission Energy</h2>
	<div class="section__body">
		<%include file="/excited_states/excited_states_diagram.mako" args="excited_states = relaxed_excited_state"/>
		<%include file="emission_results.mako" args="relaxed_excited_state = relaxed_excited_state"/>
	</div>
</div>
%if relaxed_excited_state.simulated_emission_graph is not None:
<div class="section">
	<h2 class="section__header">Emission</h2>
	<div class="section__body">
		<div class="image--absorptionGraph">
			<div class="image__aligner image__aligner--absorptionGraph">
				<img class="image__img image__img--absorptionGraph" src="${relaxed_excited_state.simulated_emission_graph.relative_path()}">
			</div>
			<div class="image__caption">Emission spectrum (simulated FWHM: ${relaxed_excited_state.simulated_emission_graph.fwhm} nm)</div>
		</div>
	</div>
</div>
%endif