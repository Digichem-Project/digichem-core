## -*- coding: utf-8 -*-

<%page args="dipole_moment, report"/>

<%namespace name="dipole_titles" file="/dipole_moment/title.mako"/>

<%!
	from silico.misc.text import text_float
%>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", dipole_moment.name + " summary")}:</div>
		Summary of the ${dipole_titles.dipole_name(dipole_moment)} properties.
	</div>
    <table class="resultsTable__table">
##      Electric TDM.
        <tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">μ</td>
		    %if dipole_moment.electric is not None:
                <td class="resultsTable__value resultsTable__cell">${text_float(dipole_moment.electric.total)} D</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
		</tr>
		<tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">μ (Gaussian-CGS)</td>
            %if dipole_moment.electric is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2e}".format(dipole_moment.electric.gaussian_cgs)} esu⋅cm</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
		<tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">θ<sub>μ,X</sub></td>
		    %if dipole_moment.electric is not None and dipole_moment.electric.X_axis_angle is not None:
		    	<td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.electric.X_axis_angle.angle, dipole_moment.electric.X_axis_angle.pretty_units)}</td>
		    %else:
		    	<td class="resultsTable__value resultsTable__cell">N/A</td>
		    %endif
		</tr>
		<tr class="resultsTable__row">
		    <td class="resultsTable__title resultsTable__cell">θ<sub>μ,XY</sub></td>
		    %if dipole_moment.electric is not None and dipole_moment.electric.XY_plane_angle is not None:
		    	<td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.electric.XY_plane_angle.angle, dipole_moment.electric.XY_plane_angle.pretty_units)}</td>
		    %else:
		    	<td class="resultsTable__value resultsTable__cell">N/A</td>
		    %endif
		</tr>
##      Magnetic TDM.
		<tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">m</td>
            %if dipole_moment.magnetic is not None:
                <td class="resultsTable__value resultsTable__cell">${text_float(dipole_moment.magnetic.total)} a.u.</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">m (Gaussian-CGS)</td>
            %if dipole_moment.magnetic is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2e}".format(dipole_moment.magnetic.gaussian_cgs)} erg⋅G<sup>-1</sup></td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">θ<sub>m,X</sub></td>
            %if dipole_moment.magnetic is not None and dipole_moment.magnetic.X_axis_angle is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.magnetic.X_axis_angle.angle, dipole_moment.magnetic.X_axis_angle.pretty_units)}</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">θ<sub>m,XY</sub></td>
            %if dipole_moment.magnetic is not None and dipole_moment.magnetic.XY_plane_angle is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.magnetic.XY_plane_angle.angle, dipole_moment.magnetic.XY_plane_angle.pretty_units)}</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">θ<sub>μ,m</sub></td>
            %if dipole_moment.angle() is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2f} {}".format(dipole_moment.angle().angle, dipole_moment.angle().pretty_units)}</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">cos(θ<sub>μ,m</sub>)</td>
            %if dipole_moment.cos_angle() is not None:
                <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.cos_angle())}</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">g<sub>lum</sub></td>
            %if dipole_moment.angle is not None:
                <td class="resultsTable__value resultsTable__cell">${text_float(dipole_moment.g_value, 3)}</td>
            %else:
                <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
    </table>
</div>