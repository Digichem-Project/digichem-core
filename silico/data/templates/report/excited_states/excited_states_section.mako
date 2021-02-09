## -*- coding: utf-8 -*-

<%page args="excited_states, report"/>

<%
ES_threshold = 100
%>

<div class="section section--fullPage">
    <h2 class="section__header">Excited States</h2>
    <div class="section__body">
        <%include file="excited_states_diagram.mako" args="excited_states = excited_states, report = report, energies_image_name = 'excited_state_energies'"/>
        <%include file="excited_states_results.mako" args="excited_states = excited_states"/>
    </div>
</div>
%if 'simulated_absorption_graph' in report.images:
<div class="section">
    <h2 class="section__header">Absorptions</h2>
    <div class="section__body">
        <div class="image--absorptionGraph">
            <div class="image__aligner image__aligner--absorptionGraph">
                <img class="image__img image__img--absorptionGraph" src="${report.relative_image('simulated_absorption_graph')}">
            </div>
            <div class="image__caption">
            Absorption spectrum (simulated Gaussian functions with FWHM: ${report.images['simulated_absorption_graph'].fwhm} eV).<br>
            %if excited_states[-1].wavelength > ES_threshold:
            <span class="warning_msg">
            Note: high energy absorption peaks are not simulated.<br>For a complete absorption spectrum, use more excited states.
            </span>
            %endif
            </div>
        </div>
    </div>
</div>
%endif
<div class="section section--fullPage">
    <h2 class="section__header">Table of Excited States</h2>
    <div class="section__body section__body--table">
    <%include file="excited_states_table.mako" args="excited_states = excited_states"/>
    </div>
</div>