## -*- coding: utf-8 -*-

<%page args="report, density_image_name"/>

<%
	density_image = report.images[density_image_name]
%>

<div class="section">
    <h2 class="section__header">${density_image.type} Density</h2>
    <div class="section__body">
        <div class="imageBlock imageBlock--multi">
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(density_image_name, 'x0y0z0')}">
                </div>
                <div class="image__caption">X/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(density_image_name, 'x90y0z0')}">
                </div>
                <div class="image__caption">X/Z plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(density_image_name, 'x0y90z0')}">
                </div>
                <div class="image__caption">Z/Y plane</div>
            </div>
            <div class="image">
                <div class="image__aligner">
                    <img class="image__img" src="${report.relative_image(density_image_name, 'x45y45z45')}">
                </div>
                <div class="image__caption">45° to axes</div>
            </div>
            <div class="imageBlock__caption">${density_image.type} density (isovalue: ${density_image.isovalue})</div>
        </div>
    </div>
</div>