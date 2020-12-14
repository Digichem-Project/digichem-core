## -*- coding: utf-8 -*-
<%!
	import math
%>

<%page args="spin_orbit_coupling"/>

<div class="section">
	<h2 class="section__header">Table of Spin-Orbit Coupling</h2>
	<div class="section__body section__body--table">
		<div class="tableBorder">
			<table class="table">
				<tr class="table__row table__row--header">
					<th class="table__header">Singlet</th>
					<th class="table__header">Triplet</th>
					<th class="table__header">SOC +1 <br>/cm<sup>-1</sup></th>
					<th class="table__header">SOC 0 <br>/cm<sup>-1</sup></th>
					<th class="table__header">SOC -1 <br>/cm<sup>-1</sup></th>
					<th class="table__header">SOC Root<br>Sum Square /cm<sup>-1</sup></th>
					<th class="table__header">H<sub>SO</sub> /eV</th>
					<th class="table__header">ΔE /eV</th>
					<th class="table__header">First Order<br>Mixing Coefficient</th>
				</tr>
				%for soc in spin_orbit_coupling:
				<tr class="table__row">
					<td class="table__header">${soc.singlet_state.multiplicity_symbol}<sub>${soc.singlet_state.multiplicity_level}</sub></td>
					<td class="table__header">${soc.triplet_state.multiplicity_symbol}<sub>${soc.triplet_state.multiplicity_level}</sub></td>
					<td class="table__header">${"{:0.4f}".format(soc.positive_one)}</td>
					<td class="table__header">${"{:0.4f}".format(soc.zero)}</td>
					<td class="table__header">${"{:0.4f}".format(soc.negative_one)}</td>
					<td class="table__header">${"{:0.4f}".format(soc.wavenumbers)}</td>
					<td class="table__header">${"{:0.4f}".format(soc.energy)}</td>
					<td class="table__header">${"{:0.4f}".format(soc.splitting_energy)}</td>
					<td class="table__header">${"{:0.4f}".format(soc.mixing_coefficient)}</td>			
				</tr>
				%endfor
			</table>
		</div>
	</div>
</div>