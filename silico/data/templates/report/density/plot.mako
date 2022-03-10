## -*- coding: utf-8 -*-
##
## Template to display a collection of four images showing some sort of density plot (eg, orbitals, differential density, NTO etc).
##
<%page args="image_name, caption, report" />


<div class="imageBlock imageBlock--multi imageBlock--orbital">
    <div class="image">
        <div class="image__aligner">
            <img class="image__img" src="${report.relative_image(image_name, 'x0y0z0')}">
        </div>
        <div class="image__caption">X/Y plane</div>
    </div>
    <div class="image">
        <div class="image__aligner">
            <img class="image__img" src="${report.relative_image(image_name, 'x90y0z0')}">
        </div>
        <div class="image__caption">X/Z plane</div>
    </div>
    <div class="image">
        <div class="image__aligner">
            <img class="image__img" src="${report.relative_image(image_name, 'x0y90z0')}">
        </div>
        <div class="image__caption">Z/Y plane</div>
    </div>
    <div class="image">
        <div class="image__aligner">
            <img class="image__img" src="${report.relative_image(image_name, 'x45y45z45')}">
        </div>
        <div class="image__caption">45° to axes</div>
    </div>
    <div class="imageBlock__caption">${caption}</div>
</div>