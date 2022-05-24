## -*- coding: utf-8 -*-

<%page args="dipole_moment, magnetic_dipole_moment, report"/>

<%namespace name="dipole_titles" file="/dipole_moment/title.mako"/>

<%!
	from silico.misc.text import text_float
%>

<%
    # The two dipole moments could be unrelated, but we generally expect them to be from the same transition.
    base_tdm = magnetic_dipole_moment
    if dipole_moment is not None:
        base_tdm = dipole_moment
        
    if base_tdm is None:
        raise ValueException("At least one of either dipole_moment or magnetic_dipole_moment must not be None")
%>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", base_tdm.name + " summary")}:</div>
		Summary of the ${dipole_titles.dipole_name(dipole_moment)} properties.
	</div>
    <table class="resultsTable__table">
##      Electric TDM.
        <tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">TEDM Total</td>
		    %if dipole_moment is not None:
                <td class="resultsTable__value resultsTable__cell">${text_float(dipole_moment.total)} D</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
		</tr>
		<tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">TEDM Total (Gaussian-CGS)</td>
            %if dipole_moment is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2e}".format(dipole_moment.gaussian_cgs)} esu⋅cm</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
		<tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">TEDM X axis angle</td>
		    %if dipole_moment is not None and dipole_moment.X_axis_angle is not None:
		    	<td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.X_axis_angle.angle, dipole_moment.X_axis_angle.pretty_units)}</td>
		    %else:
		    	<td class="resultsTable__value resultsTable__cell">N/A</td>
		    %endif
		</tr>
		<tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">TEDM XY plane angle</td>
		    %if dipole_moment is not None and dipole_moment.XY_plane_angle is not None:
		    	<td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.XY_plane_angle.angle, dipole_moment.XY_plane_angle.pretty_units)}</td>
		    %else:
		    	<td class="resultsTable__value resultsTable__cell">N/A</td>
		    %endif
		</tr>
##      Magnetic TDM.
		<tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">TMDM Total</td>
            %if magnetic_dipole_moment is not None:
                <td class="resultsTable__value resultsTable__cell">${text_float(magnetic_dipole_moment.total)} a.u.</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">TMDM Total (Gaussian-CGS)</td>
            %if magnetic_dipole_moment is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2e}".format(magnetic_dipole_moment.gaussian_cgs)} erg⋅G<sup>-1</sup></td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">TMDM X axis angle</td>
            %if magnetic_dipole_moment is not None and magnetic_dipole_moment.X_axis_angle is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(magnetic_dipole_moment.X_axis_angle.angle, magnetic_dipole_moment.X_axis_angle.pretty_units)}</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">TMDM XY plane angle</td>
            %if magnetic_dipole_moment is not None and magnetic_dipole_moment.XY_plane_angle is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(magnetic_dipole_moment.XY_plane_angle.angle, magnetic_dipole_moment.XY_plane_angle.pretty_units)}</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
    </table>
</div>