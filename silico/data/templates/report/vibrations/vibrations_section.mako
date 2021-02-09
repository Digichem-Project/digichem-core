<%page args="vibrations, report" />

%if 'simulated_IR_graph' in report.images:
<div class="section">
    <h2 class="section__header">Vibrations</h2>
    <div class="section__body">
        <div class="imageBlock imageBlock--wide">
            <div class="image__aligner image__aligner--wide">
                <img class="image__img image__img--wide" src="${report.relative_image('simulated_IR_graph')}">
            </div>
            <div class="imageBlock__caption">IR spectrum (simulated Gaussian functions with FWHM: ${report.images['simulated_IR_graph'].fwhm} cm<sup>-1</sup>)</div>
        </div>
    </div>
</div>
%endif