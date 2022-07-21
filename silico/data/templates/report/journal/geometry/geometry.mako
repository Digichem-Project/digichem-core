## -*- coding: utf-8 -*-

<%page args="report"/>

<%
	alignment = report.result.alignment
%>

<div class="content">
	<h5>Geometry</h5>
	The <div class="result"><div class="result__title">empirical formula</div> of the studied system was
	<div class="result__value"><%include file="/geometry/formula.mako", args="atoms = alignment"/></div></div>,
	corresponding to a <div class="result"><div class="result__title">molecular mass</div> of  <div class="result__value">${"{:0.2f}".format(alignment.molar_mass)} gmol<sup>-1</sup></div></div>${"." if alignment.safe_get('mass') is None else ""}
	%if alignment.safe_get('mass') is not None:
    and an <div class="result"><div class="result__title">exact mass</div>, considering only specific atomic isotopes, of <div class="result__value">${"{:0.2f}".format(alignment.mass)} gmol<sup>-1</sup></div></div>.
    %endif
    The molecular geometry was aligned to the cartesian (X, Y and Z) axes by the <div class="result"><div class="result__title">${alignment.human_method_type} (${alignment.method_type})</div></div> method${"," if "structure" in report.images else "."}
    %if "structure" in report.images:
    and the resulting atomic position are displayed in figure ${report.captions("figure", "structure")}.
    %endif
    Using this method, the <div class="result"><div class="result__title">extent of the molecular system</div> in the X, Y and Z axes (L<sub>X</sub>, L<sub>Y</sub> and L<sub>Z</sub>, corresponding to the molecular width, length and height respectively) was determined to be <div class="result__value">${"{:0.2f}".format(alignment.X_length)}, ${"{:0.2f}".format(alignment.Y_length)} and ${"{:0.2f}".format(alignment.Z_length)} Å</div> respectively</div>.
    These extensions give rise to a <div class="result"><div class="result__title">molecular linearity ratio</div> (1-(L<sub>Y</sub>/L<sub>X</sub>)) and <div class="result__title">planarity ratio</div> (1-(L<sub>X</sub>/L<sub>Y</sub>)) of <div class="result__value">${"{:0.2f}".format(alignment.get_linear_ratio())}</div> and <div class="result__value">${"{:0.2f}".format(alignment.get_planar_ratio())}</div></div> respectively.
    %if "structure" in report.images:
    <%include file="/geometry/image.mako" args="image_name = 'structure', caption = 'The molecular structure, aligned using the {} ({}) method'.format(alignment.human_method_type,  alignment.method_type), report = report" />
    %endif
</div>