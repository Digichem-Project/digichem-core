## -*- coding: utf-8 -*-

<%!
	from silico.exception import Result_unavailable_error 
%>

<%page args="molecular_orbitals"/>

<%
	# Render our alternative orbital energy diagrams (we don't actually use these in the template (yet?), just write to file.
	molecular_orbitals.get_file('HOMO_LUMO_energy_diagram').relative_path()


	# First work out the type of orbitals we are dealing with (alpha, beta, or restricted).
	orbital_type = molecular_orbitals.spin_type
	
	# Our title.
	title = "HOMO & LUMO"
	
	# Add on some more to our title if we have unrestricted orbitals.
	if orbital_type != "none":
		title += " ({})".format(orbital_type.capitalize())
		
	# Get our FMOs.
	try:
		HOMO = molecular_orbitals.get_orbital(HOMO_difference = 0)
	except Result_unavailable_error:
		HOMO = None
		
	try:
		LUMO = molecular_orbitals.get_orbital(HOMO_difference = 1)
	except Result_unavailable_error:
		LUMO = None
		
	try:
		combined_image = molecular_orbitals.HOMO_LUMO_image
	except Result_unavailable_error:
		combined_image = None
%>

<div class="section">
	<h2 class="section__header">${title}</h2>
	<div class="section__body section__body--orbital">
		%if HOMO is not None:
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${HOMO.orbital_image.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${HOMO.orbital_image.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${HOMO.orbital_image.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${HOMO.orbital_image.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">HOMO density (isovalue: ${HOMO.orbital_image.isovalue})</div>
		</div>
		%endif
		%if LUMO is not None:
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${LUMO.orbital_image.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${LUMO.orbital_image.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${LUMO.orbital_image.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${LUMO.orbital_image.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">LUMO density (isovalue: ${LUMO.orbital_image.isovalue})</div>
		</div>
		%endif
		%if combined_image is not None:
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${combined_image.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${combined_image.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${combined_image.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${combined_image.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">HOMO (${combined_image.primary_colour}) & LUMO (${combined_image.secondary_colour}) density<br>(isovalue: ${combined_image.isovalue})</div>
		</div>
		%endif
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image__aligner image__aligner--orbital_diagram">
				<img class="image_img image__img--orbital_diagram" src="${molecular_orbitals.energy_diagram.relative_path()}">
			</div>
		</div>
	</div>
</div>