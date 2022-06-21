## -*- coding: utf-8 -*-
##
## Template to display a collection of four images showing some sort of density plot (eg, orbitals, difference density, NTO etc).
##
<%page args="image_name, caption, report" />


<div class="resultImage resultImage--multi">
    <div class="resultImage__aligner">
    	<div class="result_Image__caption--subImage">A:</div>
        <img class="resultImage__image resultImage__image--subImage" src="${report.relative_image(image_name, 'x0y0z0')}">
    </div><!--
    --><div class="resultImage__aligner">
    	<div class="result_Image__caption--subImage">B:</div>
        <img class="resultImage__image resultImage__image--subImage" src="${report.relative_image(image_name, 'x90y0z0')}">
    </div><!--
    --><div class="resultImage__aligner">
    	<div class="result_Image__caption--subImage">C:</div>
        <img class="resultImage__image resultImage__image--subImage" src="${report.relative_image(image_name, 'x0y90z0')}">
    </div><!--
    --><div class="resultImage__aligner">
    	<div class="result_Image__caption--subImage">D:</div>
        <img class="resultImage__image resultImage__image--subImage" src="${report.relative_image(image_name, 'x45y45z45')}">
    </div>
    <div class="resultImage__caption resultImage__caption--multi"><div class="caption">Figure ${report.captions("figure", image_name)}:</div> ${caption}. A: In the X/Y plane, B: In the X/Z plane, C: In the Z/Y plane, D: 45° to the axes.</div>
</div>