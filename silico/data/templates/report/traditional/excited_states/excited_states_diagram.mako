## -*- coding: utf-8 -*-

<%!
    import silico.logging
%>

<%page args="excited_states, report, energies_image_name"/>

<%
    # Weasyprint hack for broken max-width.
    # Our desired max dimensions.
    max_width = 420
    max_height = 450
    try:
        dimensions = report.images[energies_image_name].get_constrained_size(max_width, max_height)
    except Exception:
        silico.logging.get_logger().error("Could not load image", exc_info = True)
        dimensions = None
    # Now we open our diagram with 
%>

%if dimensions is not None:
<div class="image__aligner image__aligner--excitedStatesDiagram">
    <img style="width: ${dimensions[0]}px; height: ${dimensions[1]}px" class="image__img image__img--excitedStatesDiagram" src="${report.relative_image(energies_image_name)}">
</div>
%endif