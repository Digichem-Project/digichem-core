## -*- coding: utf-8 -*-

<%page args="dipole_moment, report, image_name"/>

<%!
    from silico.misc.text import text_float
%>
<%namespace name="dipole_titles" file="/dipole_moment/title.mako"/>

<div class="content">
	<h5>${dipole_titles.dipole_title(dipole_moment)}</h5>
	%if dipole_moment.dipole_type == "permanent":
		The calculated <div class="result"><div class="result__title">permanent dipole moment (PDM)</div>
	%else:
		The calculated <div class="result"><div class="result__title">transition dipole moment (TDM)</div> between the ground state and the ${dipole_moment.excited_state.multiplicity_symbol}<sub>${dipole_moment.excited_state.multiplicity_level}</sub> excited state
	%endif
	%if dipole_moment.magnitude != 0:
		was <div class="result__value">${text_float(dipole_moment.magnitude)} D</div></div>, with a vector (X,Y,Z) of (${"{:0.2f}".format(dipole_moment.vector_coords[0])}, ${"{:0.2f}".format(dipole_moment.vector_coords[1])}, ${"{:0.2f}".format(dipole_moment.vector_coords[2])}) D.
		
		%if dipole_moment.origin_coords != (0,0,0):
		The origin point of the dipole vector was (${"{:0.2f}".format(dipole_moment.origin_coords[0])}, ${"{:0.2f}".format(dipole_moment.origin_coords[1])}, ${"{:0.2f}".format(dipole_moment.origin_coords[2])}).
		%endif
		%if dipole_moment.X_axis_angle is not None:
		The angle between the dipole moment vector and the X-axis was ${"{:0.2f} {}".format(dipole_moment.X_axis_angle.angle, dipole_moment.X_axis_angle.pretty_units)}${"." if dipole_moment.XY_plane_angle is None else ","}
		%endif
		%if dipole_moment.X_axis_angle is not None and dipole_moment.XY_plane_angle:
		while the
		%elif dipole_moment.XY_plane_angle is not None:
		The
		%endif
		%if dipole_moment.XY_plane_angle is not None:
		angle between the dipole moment and the XY-plane was ${"{:0.2f} {}".format(dipole_moment.XY_plane_angle.angle, dipole_moment.XY_plane_angle.pretty_units)}.
		%endif
		%if image_name in report.images:
		<%include file="/geometry/image.mako" args="image_name = image_name, caption = 'The {} (red arrow) plotted against the aligned molecular geometry'.format(dipole_titles.dipole_name(dipole_moment)), report = report" />
		%endif
	%else:
		was <div class="result__value">exactly 0 D</div></div>.
	%endif
</div>