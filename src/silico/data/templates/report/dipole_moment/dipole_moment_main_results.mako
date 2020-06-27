## -*- coding: utf-8 -*-

<%page args="dipole_moment"/>

<tr>
	<td class="results__name">Total:</td>
	<td class="results__value">${"{:0.2f}".format(dipole_moment.magnitude)} D</td>
</tr>
%if dipole_moment.X_axis_angle is not None:
<tr>
	<td class="results__name">X axis angle:</td>
	<td class="results__value">${"{:0.2f} {}".format(dipole_moment.X_axis_angle.angle, dipole_moment.X_axis_angle.pretty_units)}</td>
</tr>
%endif
%if dipole_moment.XY_plane_angle is not None:
<tr>
	<td class="results__name">XY plane angle:</td>
	<td class="results__value">${"{:0.2f} {}".format(dipole_moment.XY_plane_angle.angle, dipole_moment.XY_plane_angle.pretty_units)}</td>
</tr>
%endif