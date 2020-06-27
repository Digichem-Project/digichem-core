## -*- coding: utf-8 -*-

<%page args="energies"/>

<div class="section">
	<h2 class="section__header">${energies.energy_type} Energies</h2>
	<div class="section__body">
		<div class="imageBlock imageBlock--wide">
			<div class="image__aligner image__aligner--wide">
				<img class="image__img image__img--wide image__img--energyGraph" src="${energies.convergence_graph.relative_path()}">
			</div>
		</div>
		<%include file="/energy/energy_results.mako" args="energies = energies"/>
	</div>
</div>