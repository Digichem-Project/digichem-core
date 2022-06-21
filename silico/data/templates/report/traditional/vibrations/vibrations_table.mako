## -*- coding: utf-8 -*-
<%!
    import math
%>

<%page args="vibrations, min_frequency = None, max_frequency = None, max_num = None"/>

<%
    # Set some defaults.
    if min_frequency is None:
        min_frequency = -math.inf
    if max_frequency is None:
        max_frequency = math.inf
    
    # Filter our big list of vibrations by the given criteria:
    vibrations = [vibration for vibration in vibrations if vibration.frequency > min_frequency and vibration.frequency < max_frequency][:max_num]
%>
<div class="section">
    <h2 class="section__header">Table of Vibrational Frequencies</h2>
    <div class="section__body section__body--table">
        <div class="tableBorder">
            <table class="table">
                <tr class="table__row table__row--header">
                    %for column in ["Level", "Symmetry", "Frequency /cm<sup>-1</sup>", "Intensity /km mol<sup>-1</sup>"]:
                    <th class="table__header">${column}</th>
                    %endfor
                </tr>
                %for vibration in vibrations:
                <tr class="table__row">
                    <td class="table__cell">
                        ${vibration.level}
                    </td>
                    <td class="table__cell">
                        ${vibration.symmetry}
                    </td>
                    <td class="table__cell">
                        ${"{:0.4f}".format(vibration.frequency)}
                    </td>
                    <td class="table__cell">
                        ${"{:0.4f}".format(vibration.intensity)}
                    </td>            
                </tr>
                %endfor
            </table>
        </div>
    </div>
</div>