## -*- coding: utf-8 -*-
##
<%!
    from silico.exception import Result_unavailable_error
%>\
##
<%page args="dipole_moment, result_name = ''"/>\
##
<%
    # Get upset if we've been given nothing.
    if dipole_moment is None:
        raise Result_unavailable_error("dipole moment", "there is no dipole of the requested type")
    
    
    if dipole_moment.dipole_type == "permanent":
        title = "Permanent Dipole Moment"
    else:
        title = "{} Transition Dipole Moment".format(dipole_moment.excited_state.state_symbol)
%>\
##
##
<%include file="title.mako" args="title=title, result_name=result_name"/>
##
##
Total /D: ${"{:0.2f}".format(dipole_moment.magnitude)}
Vector x /D: ${"{:0.2f}".format(dipole_moment.vector_coords[0])}
Vector y /D: ${"{:0.2f}".format(dipole_moment.vector_coords[1])}
Vector z /D: ${"{:0.2f}".format(dipole_moment.vector_coords[2])}
x axis angle /${dipole_moment.X_axis_angle.pretty_units}: ${"{:0.2f}".format(dipole_moment.X_axis_angle.angle)}
xy plane angle /${dipole_moment.XY_plane_angle.pretty_units}: ${"{:0.2f}".format(dipole_moment.XY_plane_angle.angle)}
