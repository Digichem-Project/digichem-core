## -*- coding: utf-8 -*-

## Template to display spin density.

<%page args="result" />

<%
	# Decide which spin value (up/alpha or down/beta) corresponds to the 'hole' and 'electron'.
	# We do this by comparing the alpha and beta HOMO values, the higher one is the electron.
	alpha_energy = result.molecular_orbitals.HOMO_energy
	beta_energy = result.beta_orbitals.HOMO_energy
	
	if alpha_energy > beta_energy:
		# Alpha/up/positive is electron
		alpha_designation = "electron"
		beta_designation = "hole"
	else:
		# Beta/down/negative is electron.
		alpha_designation = "hole"
		beta_designation = "electron"
%>


## We still use flexbox for some layout here, which breaks pagination in weasyprint. For now, we fix by splitting our list into fours (which is the number of MOs per page).

## This solution was taken from https://stackoverflow.com/questions/6614891/turning-a-list-into-nested-lists-in-python/6615011
<div class="section">
	<h2 class="section__header">Spin Density</h2>
	<div class="section__body section__body--orbital">
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_positive.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_positive.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_positive.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_positive.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">Positive spin density (${alpha_designation}) (isovalue: ${result.spin_image_positive.isovalue})</div>
		</div>
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_negative.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_negative.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_negative.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_negative.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">Negative spin density (${beta_designation}) (isovalue: ${result.spin_image_negative.isovalue})</div>
		</div>
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_both.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_both.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_both.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${result.spin_image_both.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">Positive (${alpha_designation}) (${result.spin_image_both.primary_colour}) & negative (${beta_designation}) (${result.spin_image_both.secondary_colour}) spin density (isovalue: ${result.spin_image_both.isovalue})</div>
		</div>
	</div>
</div>