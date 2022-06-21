## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
    from silico.misc.text import text_float
    from silico.result.angle import Angle
%>

<%
    # TODO: as in many other places, this way of determining angle units is not very smart.
    angle_units = Angle.units_to_pretty_units(Angle._default_angle_units)
%>

<div class="resultsTable">
    <div class="resultsTable__caption">
        <div class="caption">Table ${report.captions("table", "transition_dipole_moments")}:</div> Properties of the calculated transition dipole moments.
        [a]: The electric transition dipole moment (TEDM), in Debye (D).
        [b]: Angle between the TEDM and the x-axis of the molecule.
        [c]: Angle between the TEDM and xy-plane of the molecule.
        [d]: The magnetic transition dipole moment (TMDM), in atomic units (au).
        [e]: Angle between the TMDM and the x-axis of the molecule.
        [f]: Angle between the TMDM and xy-plane of the molecule.
        [g]: The TEDM, in Gaussian CGS (centimetre, gram, second) units.
        [h]: The TMDM, in Gaussian CGS (centimetre, gram, second) units.
        [i]: The angle between the electric and magnetic transition dipole moments, in Gaussian CGS units.
        [j]: The cosine of the angle between the electric and magnetic transition dipole moments, in Gaussian CGS units.
        [k]: The dissymmetry factor of the transition dipole moment.
        
    </div>
    <table class="resultsTable__table">
        ##############
        ## Headers. ##
        ##############
        <tr class="resultTable__row">
##          Electric
            ##<td class="resultsTable__title resultsTable__cell">Number</td>
            <td class="resultsTable__title resultsTable__cell">Excited<br>State</td>
            <td class="resultsTable__title resultsTable__cell">μ<sup>[a]</sup> Vector<br>/D</td>
            <td class="resultsTable__title resultsTable__cell">μ<sup>[a]</sup><br>/D</td>
            <td class="resultsTable__title resultsTable__cell">θ<sub>μ,x</sub><sup>[b]</sup><br>/${angle_units}</td>
            <td class="resultsTable__title resultsTable__cell">θ<sub>μ,xy</sub><sup>[c]</sup><br>/${angle_units}</td>
##          Magnetic
            <td class="resultsTable__title resultsTable__cell">m<sup>[d]</sup> Vector<br>/au</td>
            <td class="resultsTable__title resultsTable__cell">m<sup>[d]</sup><br>/au</td>
            <td class="resultsTable__title resultsTable__cell">θ<sub>m,x</sub><sup>[e]</sup><br>/${angle_units}</td>
            <td class="resultsTable__title resultsTable__cell">θ<sub>m,xy</sub><sup>[f]</sup><br>/${angle_units}</td>
##          Comparison
            <td class="resultsTable__title resultsTable__cell">μ<sup>[g]</sup><br>/esu⋅cm</td>
            <td class="resultsTable__title resultsTable__cell">m<sup>[h]</sup><br>/erg⋅G<sup>-1</sup></td>
            <td class="resultsTable__title resultsTable__cell">θ<sub>μ,m</sub><sup>[i]</sup><br>/${angle_units}</td>
            <td class="resultsTable__title resultsTable__cell">cos(θ<sub>μ,m</sub>)<sup>[j]</sup></td>
            <td class="resultsTable__title resultsTable__cell">g<sub>lum</sub><sup>[k]</sup></td>
        </tr>
        
        ###########
        ## Data. ##
        ###########
        %for excited_state, dipole_moment in [(excited_state, excited_state.transition_dipole_moment) for excited_state in report.result.excited_states]:
        <tr class="resultTable__row">
##          Electric
            %if dipole_moment.electric is not None:
            <td class="resultsTable__value resultsTable__cell">${excited_state.multiplicity_symbol}<sub>${excited_state.multiplicity_level}</sub></td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}, {:0.2f}, {:0.2f}".format(*dipole_moment.electric.vector_coords)}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(dipole_moment.electric.total)}</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.electric.X_axis_angle.angle)}</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.electric.XY_plane_angle.angle)}</td>
            %else:
                %for i in range(0, 4):
                <td class="resultsTable__value resultsTable__cell">N/A</td>
                %endfor
            %endif
##          Magnetic
            %if dipole_moment.magnetic is not None:
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}, {:0.2f}, {:0.2f}".format(*dipole_moment.magnetic.vector_coords)}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(dipole_moment.magnetic.total)}</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.magnetic.X_axis_angle.angle)}</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.magnetic.XY_plane_angle.angle)}</td>
            %else:
                %for i in range(0, 4):
                <td class="resultsTable__value resultsTable__cell">N/A</td>
                %endfor
            %endif
##          Comparison
            <td class="resultsTable__value resultsTable__cell">${"{:0.2e}".format(dipole_moment.electric.gaussian_cgs) if dipole_moment.electric is not None else "N/A"}</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2e}".format(dipole_moment.magnetic.gaussian_cgs) if dipole_moment.magnetic is not None else "N/A"}</td>
            %if dipole_moment.magnetic is not None and dipole_moment.electric is not None:
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.angle().angle)}</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(dipole_moment.cos_angle())}</td>
            <td class="resultsTable__value resultsTable__cell">${text_float(dipole_moment.g_value, 3)}</td>
            %else:
                %for i in range(0, 3):
                <td class="resultsTable__value resultsTable__cell">N/A</td>
                %endfor
            %endif
        </tr>
        %endfor
    </table>
</div>