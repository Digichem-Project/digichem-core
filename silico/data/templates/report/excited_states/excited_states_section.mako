## -*- coding: utf-8 -*-

<%page args="excited_states"/>

<div class="section section--fullPage">
	<h2 class="section__header">Excited States</h2>
	<div class="section__body">
		##<div class="image__aligner image__aligner--excitedStatesDiagram">
		##	<img class="image__img image__img--excitedStatesDiagram" src="${excited_states.excited_states_diagram.relative_path()}">
		##</div>
		<%include file="excited_states_diagram.mako" args="excited_states = excited_states"/>
		<%include file="excited_states_results.mako" args="excited_states = excited_states"/>
	</div>
</div>
%if excited_states.simulated_absorption_graph.relative_path() is not None:
<div class="section">
	<h2 class="section__header">Absorptions</h2>
	<div class="section__body">
		<div class="image--absorptionGraph">
			<div class="image__aligner image__aligner--absorptionGraph">
				<img class="image__img image__img--absorptionGraph" src="${excited_states.simulated_absorption_graph.relative_path()}">
			</div>
			<div class="image__caption">Absorption spectrum (simulated Gaussian functions with FWHM: ${excited_states.simulated_absorption_graph.fwhm} eV)</div>
		</div>
	</div>
</div>
%endif
<div class="section section--fullPage">
	<h2 class="section__header">Table of Excited States</h2>
	<div class="section__body section__body--table">
	<%include file="excited_states_table.mako" args="excited_states = excited_states"/>
	</div>
</div>