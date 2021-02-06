## -*- coding: utf-8 -*-

## Template to display images of generic orbitals (ie, not the FMOs)

<%page args="molecular_orbitals, report" />

## We still use flexbox for some layout here, which breaks pagination in weasyprint. For now, we fix by splitting our list into fours (which is the number of MOs per page).

## This solution was taken from https://stackoverflow.com/questions/6614891/turning-a-list-into-nested-lists-in-python/6615011
%for mo_group in [molecular_orbitals[i:i+4] for i in range(0, len(molecular_orbitals), 4)]:
<div class="section section--fullPage">
	<h2 class="section__header">${", ".join([orbital.label for orbital in mo_group])}</h2>
	<div class="section__body section__body--orbital">
		%for orbital in mo_group:
		<div class="imageBlock imageBlock--multi imageBlock--orbital">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${report.relative_image(orbital.label, 'x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${report.relative_image(orbital.label, 'x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${report.relative_image(orbital.label, 'x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${report.relative_image(orbital.label, 'x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">${orbital.label} density (isovalue: ${report.images[orbital.label].isovalue})</div>
		</div>
		%endfor
	</div>
</div>
%endfor