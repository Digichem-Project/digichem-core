## -*- coding: utf-8 -*-

<%page args="dipole_moment"/>

<%
	# First work out our title, which changes slightly depending on whether this is the ground or excited state dipole.
	if dipole_moment.dipole_type == "permanent":
		# This is the ground state dipole.
		dipole_title = 'Permanent Dipole Moment'
	else:
		# This is an excited states dipole.
		dipole_title = 'Transition ({}<sub>{}</sub>) Dipole Moment'.format(dipole_moment.excited_state.multiplicity_symbol, dipole_moment.excited_state.multiplicity_level) 
%>

<div class="section">
	<h2 class="section__header">${dipole_title}</h2>
	<div class="section__body">
		%if dipole_image is not None:
		<div class="imageBlock imageBlock--multi">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${dipole_moment.dipole_image.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${dipole_moment.dipole_image.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${dipole_moment.dipole_image.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${dipole_moment.dipole_image.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">Aligned structure (dipole moment in red)</div>
		</div>
		%endif
		<div class="resultsContainer">
			<div class="reportHeader reportHeader--minor reportHeader--results">Dipole Moment</div>
			<table class="results">
				<tr>
					<td class="results__name">Origin X:</td>
					<td class="results__value">${"{:0.2f}".format(dipole_moment.origin_coords[0])} D</td>
				</tr>
				<tr>
					<td class="results__name">Origin Y:</td>
					<td class="results__value">${"{:0.2f}".format(dipole_moment.origin_coords[1])} D</td>
				</tr>
				<tr>
					<td class="results__name">Origin Z:</td>
					<td class="results__value">${"{:0.2f}".format(dipole_moment.origin_coords[2])} D</td>
				</tr>
				<tr>
					<td class="results__name">Vector X:</td>
					<td class="results__value">${"{:0.2f}".format(dipole_moment.vector_coords[0])} D</td>
				</tr>
				<tr>
					<td class="results__name">Vector Y:</td>
					<td class="results__value">${"{:0.2f}".format(dipole_moment.vector_coords[1])} D</td>
				</tr>
				<tr>
					<td class="results__name">Vector Z:</td>
					<td class="results__value">${"{:0.2f}".format(dipole_moment.vector_coords[2])} D</td>
				</tr>
				<%include file="dipole_moment_main_results.mako" args="dipole_moment = dipole_moment"/>
			</table>
		</div>
	</div>
</div>