## -*- coding: utf-8 -*-

<%page args="alignment"/>

<div class="section">
	<h2 class="section__header">Geometry</h2>
	<div class="section__body">
		<div class="imageBlock imageBlock--multi">
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${alignment.structure_image.relative_path('x0y0z0')}">
				</div>
				<div class="image__caption">X/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${alignment.structure_image.relative_path('x90y0z0')}">
				</div>
				<div class="image__caption">X/Z plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${alignment.structure_image.relative_path('x0y90z0')}">
				</div>
				<div class="image__caption">Z/Y plane</div>
			</div>
			<div class="image">
				<div class="image__aligner">
					<img class="image__img" src="${alignment.structure_image.relative_path('x45y45z45')}">
				</div>
				<div class="image__caption">45° to axes</div>
			</div>
			<div class="imageBlock__caption">Aligned structure</div>
		</div>
		<%include file="geometry_results.mako" args="alignment = alignment"/>
	</div>
</div>