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
    
    
    title = "{} Transition Dipole Moment".format(dipole_moment.excited_state.state_symbol)
    
    # The two types of TDM.
    electric = dipole_moment.electric
    magnetic = dipole_moment.magnetic
%>\
##
##
<%include file="title.mako" args="title=title, result_name=result_name"/>
##
##
## Electric
μ /D: ${"{:0.2f}".format(electric.magnitude) if electric is not None else "N/A"}
μ x /D: ${"{:0.2f}".format(electric.vector_coords[0]) if electric is not None else "N/A"}
μ y /D: ${"{:0.2f}".format(electric.vector_coords[1]) if electric is not None else "N/A"}
μ z /D: ${"{:0.2f}".format(electric.vector_coords[2]) if electric is not None else "N/A"}
θ(μ,x) /${electric.X_axis_angle.pretty_units if electric is not None else ""}: ${"{:0.2f}".format(electric.X_axis_angle.angle) if electric is not None else "N/A"}
θ(μ,xy) /${electric.XY_plane_angle.pretty_units if electric is not None else ""}: ${"{:0.2f}".format(electric.XY_plane_angle.angle) if electric is not None else "N/A"}
## Magnetic
m /au: ${"{:0.2f}".format(magnetic.magnitude) if magnetic is not None else "N/A"}
m x /au: ${"{:0.2f}".format(magnetic.vector_coords[0]) if magnetic is not None else "N/A"}
m y /au: ${"{:0.2f}".format(magnetic.vector_coords[1]) if magnetic is not None else "N/A"}
m z /au: ${"{:0.2f}".format(magnetic.vector_coords[2]) if magnetic is not None else "N/A"}
θ(m,x) /${magnetic.X_axis_angle.pretty_units if magnetic is not None else ""}: ${"{:0.2f}".format(magnetic.X_axis_angle.angle) if magnetic is not None else "N/A"}
θ(m,xy) /${magnetic.XY_plane_angle.pretty_units if magnetic is not None else ""}: ${"{:0.2f}".format(magnetic.XY_plane_angle.angle) if magnetic is not None else "N/A"}
## Comparison
μ /esu⋅cm: ${"{:0.2e}".format(electric.gaussian_cgs) if electric is not None else "N/A"}
m /erg⋅G-1: ${"{:0.2e}".format(magnetic.gaussian_cgs) if magnetic is not None else "N/A"}
%if electric is not None and magnetic is not None:
θ(μ,m) /${dipole_moment.angle().pretty_units}: ${"{:0.2e}".format(dipole_moment.angle().angle)}
cos(θ(μ,m)): ${"{:0.2e}".format(dipole_moment.cos_angle())}
g(lum): ${"{:0.3e}".format(dipole_moment.g_value)}
%else:
θ(μ,m) /: N/A
cos(θ(μ,m)): N/A
g(lum): N/A
%endif