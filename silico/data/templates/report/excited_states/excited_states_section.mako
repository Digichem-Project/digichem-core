## -*- coding: utf-8 -*-

<%page args="excited_states"/>

<%
# The number of ES below which we'll issue a warning.
#ES_threshold = 50
ES_threshold = 100
# Excited states with oscillator strength > 0.
#bright_states = [state for state in excited_states if state.oscillator_strength is not None and state.oscillator_strength > 0]
%>

<div class="section section--fullPage">
	<h2 class="section__header">Excited States</h2>
	<div class="section__body">
		##<div class="image__aligner image__aligner--excitedStatesDiagram">
		##	<img class="image__img image__img--excitedStatesDiagram" src="${excited_states.excited_states_diagram.relative_path()}">
		##</div>
		<%include file="excited_states_diagram.mako" args="excited_states = excited_states"/>
		<%include file="excited_states_results.mako" args="excited_states = excited_states"/>
	</div>
</div>
%if excited_states.simulated_absorption_graph.relative_path() is not None:
<div class="section">
	<h2 class="section__header">Absorptions</h2>
	<div class="section__body">
		<div class="image--absorptionGraph">
			<div class="image__aligner image__aligner--absorptionGraph">
				<img class="image__img image__img--absorptionGraph" src="${excited_states.simulated_absorption_graph.relative_path()}">
			</div>
			<div class="image__caption">
			Absorption spectrum (simulated Gaussian functions with FWHM: ${excited_states.simulated_absorption_graph.fwhm} eV).<br>
			%if excited_states[-1].wavelength > ES_threshold:
			<span class="warning_msg">
			Note: high energy absorption peaks are not simulated.<br>For a complete absorption spectrum, use more excited states.
			</span>
			%endif
## 			%if len(bright_states) < ES_threshold:
## 				Generated from ${len(bright_states)} bright (f > 0) excited states, <span class="warning_msg">include more to improve accuracy.</span>
## 			%else:
## 				Generated from ${len(bright_states)} bright (f > 0) excited states.
## 			%endif
			</div>
		</div>
## 		%if len(bright_states) < ES_threshold:
## 		<div class="warning_msg">
## 		Warning: Absorption spectrum simulated with low number of excited states (<${ES_threshold})
## 		</div>
## 		%endif
	</div>
</div>
%endif
<div class="section section--fullPage">
	<h2 class="section__header">Table of Excited States</h2>
	<div class="section__body section__body--table">
	<%include file="excited_states_table.mako" args="excited_states = excited_states"/>
	</div>
</div>