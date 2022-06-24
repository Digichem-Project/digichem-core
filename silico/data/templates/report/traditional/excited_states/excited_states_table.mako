## -*- coding: utf-8 -*-

<%page args="excited_states"/>

<%!
    import colour
    import numpy
    from silico.result.excited_state import Excited_state
%>

<div class=tableBorder>
    <table class="table">
        <tr class="table__row table__row--header">
            %for column in ["Level", "Symbol", "Symmetry", "Energy<br>/eV", "Wavelength<br>/nm", "Colour,<br>CIE (x,y)", "Oscillator<br>Strength", "Transitions<br>(probability)"]:
            <th class="table__header">${column}</th>
            %endfor
        </tr>
        %for excited_state in excited_states:
        <tr class="table__row">
            <td class="table__cell">
                ${excited_state.level}
            </td>
            <td class="table__cell">
                ${excited_state.multiplicity_symbol}<sub>${excited_state.multiplicity_level}</sub>
            </td>
            <td class="table__cell">
                ${excited_state.symmetry}
            </td>
            <td class="table__cell">
                ${"{:0.4f}".format(excited_state.energy)}
            </td>
            <td class="table__cell">
                ${"{:0.2f}".format(excited_state.wavelength)}
            </td>
            <td class="table__cell">
                ${excited_state.color}<br>
                <%include file="color_box.mako" args="rgb = excited_state.rgb"/><br>
                (${"{:0.2f}".format(excited_state.CIE_xy[0])}, ${"{:0.2f}".format(excited_state.CIE_xy[1])})
            </td>
##             <td class="table__cell">
##                 (${"{:0.2f}".format(excited_state.CIE_xy[0])}, ${"{:0.2f}".format(excited_state.CIE_xy[1])})
##             </td>
            <td class="table__cell">
                ${"{:0.4f}".format(excited_state.oscillator_strength)}
            </td>
            <td class="table__cell">
                %for transition in excited_state.transitions:
                <div class="etTable__transition">
                    ${transition.starting_mo.label} → ${transition.ending_mo.label} (${"{:0.2f}".format(transition.probability)})
                </div>
                %endfor
            </td>
        </tr>
        %endfor
    </table>
</div>