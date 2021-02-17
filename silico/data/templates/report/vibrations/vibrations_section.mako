<%page args="vibrations, report" />

%if 'simulated_IR_graph' in report.images:
<%
	peaks = sorted(list(set([int(peak) for peak in report.images['simulated_IR_graph'].peaks])))
	peaks = ["{}".format(peak) for peak in peaks]
%>
<div class="section">
    <h2 class="section__header">Vibrations</h2>
    <div class="section__body">
        <div class="imageBlock imageBlock--wide">
            <div class="image__aligner image__aligner--wide">
                <img class="image__img image__img--wide" src="${report.relative_image('simulated_IR_graph')}">
            </div>
            <div class="imageBlock__caption">
                IR spectrum (simulated Gaussian functions with FWHM: ${report.images['simulated_IR_graph'].fwhm} cm<sup>-1</sup>)<br>
                Peaks /cm<sup>-1</sup>: ${", ".join(peaks)}.
            </div>
        </div>
    </div>
</div>
%endif