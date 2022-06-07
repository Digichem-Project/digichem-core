## -*- coding: utf-8 -*-

<%page args="dipole_moment, units, report"/>

<%namespace name="dipole_titles" file="/dipole_moment/title.mako"/>

<%!
    from silico.misc.text import text_float
%>

The calculated
<div class="result"><div class="result__title">
    %if dipole_moment.dipole_type == "transition" and dipole_moment.electromagnetic_type == "electric":
       electric
    %elif dipole_moment.dipole_type == "transition":
       magnetic
    %endif
    ${dipole_titles.dipole_name(dipole_moment)}
</div>
%if dipole_moment.magnitude != 0:
    was <div class="result__value">${text_float(dipole_moment.magnitude)} ${units}</div></div>, with a vector (x,y,z) of ${"{:0.2f}".format(dipole_moment.vector_coords[0])}, ${"{:0.2f}".format(dipole_moment.vector_coords[1])}, ${"{:0.2f}".format(dipole_moment.vector_coords[2])} ${units}.
    
    %if dipole_moment.origin_coords != (0,0,0):
        The origin point of the dipole vector was (${"{:0.2f}".format(dipole_moment.origin_coords[0])}, ${"{:0.2f}".format(dipole_moment.origin_coords[1])}, ${"{:0.2f}".format(dipole_moment.origin_coords[2])}).
    %endif
    
    %if dipole_moment.X_axis_angle is not None:
        The angle between the dipole moment vector and the x-axis was ${"{:0.2f} {}".format(dipole_moment.X_axis_angle.angle, dipole_moment.X_axis_angle.pretty_units)}${"." if dipole_moment.XY_plane_angle is None else ","}
    %endif
    
    %if dipole_moment.X_axis_angle is not None and dipole_moment.XY_plane_angle:
        while the
    %elif dipole_moment.XY_plane_angle is not None:
        The
    %endif
    
    %if dipole_moment.XY_plane_angle is not None:
        angle between the dipole moment and the xy-plane was ${"{:0.2f} {}".format(dipole_moment.XY_plane_angle.angle, dipole_moment.XY_plane_angle.pretty_units)}.
    %endif
%else:
    was <div class="result__value">exactly 0 ${units}</div></div>.
%endif