## -*- coding: utf-8 -*-

<%page args="dipole_moment, report"/>

<%namespace name="dipole_titles" file="/dipole_moment/title.mako"/>

<%!
    from silico.misc.text import text_float
%>

<%
    image_name = 'dipole_moment'
%>

<div class="content">
	<h5>${dipole_titles.dipole_title(dipole_moment)}</h5>
	<%include file="/dipole_moment/analysis.mako" args="dipole_moment = dipole_moment, units = 'D', report = report"/>
	%if image_name in report.images and dipole_moment.magnitude != 0:
    A plot of the permanent dipole moment is shown in figure ${report.captions("figure", image_name)}.
	<%include file="/geometry/image.mako" args="image_name = image_name, caption = 'The {} ({} arrow) plotted against the aligned molecular geometry with a scale of 1 Å = {:.1f} D'.format(capture(dipole_titles.dipole_name, dipole_moment), report.images[image_name].electric_arrow_colour, 1 / report.images[image_name].scaling), report = report" />
	%endif
</div>