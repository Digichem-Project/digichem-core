## -*- coding: utf-8 -*-

<%page args="dipole_moment, report"/>

<%namespace name="dipole_titles" file="/dipole_moment/title.mako"/>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", dipole_moment.name + " summary")}:</div>
		Summary of the ${dipole_titles.dipole_name(dipole_moment)} properties.
	</div>
    <table class="resultsTable__table">
        <tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">Total</td>
		    <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.magnitude)} D</td>
		</tr>
		<tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">X axis angle</td>
		    %if dipole_moment.X_axis_angle is not None:
		    	<td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.X_axis_angle.angle, dipole_moment.X_axis_angle.pretty_units)}</td>
		    %else:
		    	<td class="resultsTable__value resultsTable__cell">N/A</td>
		    %endif
		</tr>
		<tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">XY plane angle</td>
		    %if dipole_moment.XY_plane_angle is not None:
		    	<td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.XY_plane_angle.angle, dipole_moment.XY_plane_angle.pretty_units)}</td>
		    %else:
		    	<td class="resultsTable__value resultsTable__cell">N/A</td>
		    %endif
		</tr>
    </table>
</div>