## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	from silico.misc.text import text_integer, andjoin
	import inflect
%>

<%
	inflector = inflect.engine()
%>

<div class="content">
	<h5>Vibrational Frequencies</h5>
	The energies of a total of ${text_integer(len(report.result.vibrations))} vibrational transitions were calculated and vibrational absorption peaks were simulated using a gaussian function with full-width at half maximum (FHWM) of ${report.images['simulated_IR_graph'].fwhm} cm<sup>-1</sup>.
	From this analysis the <div class="result"><div class="result__title">${text_integer(len(report.images['simulated_IR_graph'].selected_peaks(0, 5)))} most intense vibrational peaks</div> were found at <div class="result__value">${andjoin(report.images['simulated_IR_graph'].selected_peaks(0, 5))} cm<sup>-1</sup></div></div>.
	The full simulated vibrational frequency spectrum is shown in figure ${report.captions("figure", "simulated_IR_graph")}.
	Finally there ${inflector.plural("was", len(report.result.vibrations.negative_frequencies))} <div class="result"><div class="result__value">${text_integer(len(report.result.vibrations.negative_frequencies))}</div> <div class="result__title">calculated negative frequencies</div></div>.
</div>
%if len(report.result.vibrations) > 0 and report.images['simulated_IR_graph'].safe_get_file() is not None:
<div class="resultImage resultImage--large resultImage--graph">
	<img class="resultImage__image" src="${report.relative_image('simulated_IR_graph')}">
	<div class="resultImage__caption"><div class="caption">Figure ${report.captions("figure", 'simulated_IR_graph')}:</div> Graph of simulated vibrational spectrum. Calculated vibrational frequencies are shown as vetical black bars while simulated peaks with a gaussian function with FHWM: ${report.images['simulated_IR_graph'].fwhm} cm<sup>-1</sup> are shown as a blue line. Peaks can be found at: ${andjoin(report.images['simulated_IR_graph'].selected_peaks())} cm<sup>-1</sup>.</div>
</div>
%endif