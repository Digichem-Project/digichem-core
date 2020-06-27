## -*- coding: utf-8 -*-

<%!
	from logging import getLogger
	import silico
%>

<%page args="excited_states"/>

<%
	# Weasyprint hack for broken max-width.
	# Our desired max dimensions.
	max_width = 420
	max_height = 450
	try:
		dimensions = excited_states.excited_states_diagram.get_constrained_size(max_width, max_height)
	except Exception:
		getLogger(silico.logger_name).error("Could not load image", exc_info = True)
		dimensions = None
	# Now we open our diagram with 
%>

%if dimensions is not None:
<div class="image__aligner image__aligner--excitedStatesDiagram">
	<img style="width: ${dimensions[0]}px; height: ${dimensions[1]}px" class="image__img image__img--excitedStatesDiagram" src="${excited_states.excited_states_diagram.relative_path()}">
</div>
%endif