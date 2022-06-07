## -*- coding: utf-8 -*-

<%page args="dipole_moment, report"/>

<%namespace name="dipole_titles" file="/dipole_moment/title.mako"/>

<%!
    from silico.misc.text import text_float
%>

<%
    image_name = '{}_dipole'.format(dipole_moment.excited_state.state_symbol)
%>

<div class="content">
	<h5>${dipole_titles.dipole_title(dipole_moment)}</h5>
##  The content of this section changes depending on whether we have both TEDM and TMDM or just one.
	%if dipole_moment.electric is not None and dipole_moment.magnetic is not None:
##      Probably the most common case; both dipoles available.
    	The calculated <div class="result"><div class="result__title">electric (TEDM, μ)</div> and <div class="result__title">magnetic (TMDM, m)</div> transition dipole moments between the ground state and the ${dipole_moment.excited_state.multiplicity_symbol}<sub>${dipole_moment.excited_state.multiplicity_level}</sub> excited state were <div class="result__value">${text_float(dipole_moment.electric.magnitude)} D</div> and <div class="result__value">${text_float(dipole_moment.magnetic.magnitude)} au</div></div> respectively. The corresponding vector components (x,y,z) were μ = ${"{:0.2f}".format(dipole_moment.electric.vector_coords[0])}, ${"{:0.2f}".format(dipole_moment.electric.vector_coords[1])}, ${"{:0.2f}".format(dipole_moment.electric.vector_coords[2])} D and m = ${"{:0.2f}".format(dipole_moment.magnetic.vector_coords[0])}, ${"{:0.2f}".format(dipole_moment.magnetic.vector_coords[1])}, ${"{:0.2f}".format(dipole_moment.magnetic.vector_coords[2])} au.
##      Dipole Vs Geometry angle info.
##      The only conceivable reason why angle info would be missing is because there are no atoms, in which case it will be missing for both, so this check is fine.
        %if dipole_moment.electric.X_axis_angle is not None and dipole_moment.electric.X_axis_angle is not None:
            In comparison to the molecular geometry, the angle between each dipole moment and the longest axis of the molecule (the x-axis) was θ<sub>μ,x</sub> = ${"{:0.2f} {}".format(dipole_moment.electric.X_axis_angle.angle, dipole_moment.electric.X_axis_angle.pretty_units)} and θ<sub>m,x</sub> = ${"{:0.2f} {}".format(dipole_moment.magnetic.X_axis_angle.angle, dipole_moment.magnetic.X_axis_angle.pretty_units)}, while the angle between each dipole moment and the xy-plane was θ<sub>μ,xy</sub> = ${"{:0.2f} {}".format(dipole_moment.electric.XY_plane_angle.angle, dipole_moment.electric.XY_plane_angle.pretty_units)} and θ<sub>m,xy</sub> = ${"{:0.2f} {}".format(dipole_moment.magnetic.XY_plane_angle.angle, dipole_moment.magnetic.XY_plane_angle.pretty_units)}.
    	%endif
##      Dipole Vs Dipole angle info.
        In Gaussian-CGS units, in which the magnetic and electric transition dipole moments can be directly compared, the magnitude of each dipole moment was μ = ${"{:0.2e}".format(dipole_moment.electric.gaussian_cgs)} esu⋅cm and m = ${"{:0.2e}".format(dipole_moment.magnetic.gaussian_cgs)} erg⋅G<sup>-1</sup>, while the <div class="result"><div class="result__title">angle between the two dipole moments</div> was θ<sub>μ,m</sub> = <div class="result__value">${"{:0.2f} {}".format(dipole_moment.angle().angle, dipole_moment.angle().pretty_units)}</div></div>. Correspondingly, the cosine of the angle was cos(θ<sub>μ,m</sub>) = ${"{:0.2f}".format(dipole_moment.cos_angle())}, and the <div class="result"><div class="result__title">dissymmetry factor</div> of the excited state transition was g<sub>lum</sub> = <div class="result__value">${text_float(dipole_moment.g_value, 3)}</div></div>.
        %if image_name in report.images:
            A plot of the electric and magnetic transition dipole moments is shown in figure ${report.captions("figure", image_name)}.
        %endif
    %elif dipole_moment.electric is not None:
##      Only the electric available.
        <%include file="/dipole_moment/analysis.mako", args="dipole_moment = dipole_moment.electric, units = 'D', report = report"/>
        %if image_name in report.images:
            A plot of the electric transition dipole moment is shown in figure ${report.captions("figure", image_name)}.
        %endif
    %else:
##      Only the magnetic available.
        <%include file="/dipole_moment/analysis.mako", args="dipole_moment = dipole_moment.magnetic, units = 'au', report = report"/>
        %if image_name in report.images:
            A plot of the magnetic transition dipole moment is shown in figure ${report.captions("figure", image_name)}.
        %endif
    %endif
##  Image.
	%if image_name in report.images:
        <%include file="/geometry/image.mako" args="image_name = image_name, caption = 'The electric ({} arrow) and magnetic ({} arrow) {} plotted against the aligned molecular geometry with a scale of 1 Å = {:.1f} D = {:.1f} au'.format(report.images[image_name].electric_arrow_colour, report.images[image_name].magnetic_arrow_colour, capture(dipole_titles.dipole_name, dipole_moment), 1 / report.images[image_name].scaling, 1 / report.images[image_name].magnetic_scaling,), report = report" />
    %endif
</div>