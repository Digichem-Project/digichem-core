## -*- coding: utf-8 -*-

<%page args="relaxed_excited_state, report, energies_image_name, graph_image_name"/>

<div class="section">
    <h2 class="section__header">${relaxed_excited_state.transition_type.capitalize()} ${"{}<sub>{}</sub>".format(relaxed_excited_state.multiplicity_symbol, relaxed_excited_state.multiplicity_level)} Emission Energy</h2>
    <div class="section__body">
        <%include file="/excited_states/excited_states_diagram.mako" args="excited_states = relaxed_excited_state, report = report, energies_image_name = energies_image_name"/>
        <%include file="emission_results.mako" args="relaxed_excited_state = relaxed_excited_state"/>
    </div>
</div>
%if graph_image_name in report.images and report.images[graph_image_name].safe_get_file() is not None:
<%
	peaks = sorted(list(set([int(peak) for peak in report.images[graph_image_name].peaks])))
	peaks = ["{}".format(peak) for peak in peaks]
%>
<div class="section">
    <h2 class="section__header">Emission</h2>
    <div class="section__body">
        <div class="image--absorptionGraph">
            <div class="image__aligner image__aligner--absorptionGraph">
                <img class="image__img image__img--absorptionGraph" src="${report.relative_image(graph_image_name)}">
            </div>
            <div class="imageBlock__caption">
            	Emission spectrum (simulated Gaussian functions with FWHM: ${report.images[graph_image_name].fwhm} eV)<br>
            	Peaks /nm: ${", ".join(peaks)}.
            </div>
        </div>
    </div>
</div>
%endif