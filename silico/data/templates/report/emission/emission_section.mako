## -*- coding: utf-8 -*-

<%page args="emission, report"/>
<%
	# Decide on the images to use.
	energies_image_name = '{}_{}_emission_energies'.format(emission.transition_type, emission.state_symbol)
	graph_image_name = 'simulated_{}_{}_emission_graph'.format(emission.transition_type, emission.state_symbol)
%>

<div class="section">
    <h2 class="section__header">${emission.transition_type.capitalize()} ${"{}<sub>{}</sub>".format(emission.multiplicity_symbol, emission.multiplicity_level)} Emission Energy</h2>
    <div class="section__body">
        <%include file="/excited_states/excited_states_diagram.mako" args="excited_states = emission, report = report, energies_image_name = energies_image_name"/>
        <%include file="emission_results.mako" args="emission = emission"/>
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