## -*- coding: utf-8 -*-

<%page args="excited_states" />

<%
    # Singlet/triplet splitting energy, only available if both singlets and triplets have been calculated.
    try:
        dest = excited_states.singlet_triplet_energy
    except Exception:
        dest = None
        
    # S1, the lowest energy singlet excited state. We include this in our summary because it is the most important excited state for most structures re. fluorescence.
    try:
        S1 = excited_states.get_state("S(1)")
    except Exception:
        S1 = None
        
    # T1, the lowest energy singlet excited state. We include this in our summary because it is the most important excited state for most structures re. phosphorescence.
    try:
        T1 = excited_states.get_state("T(1)")
    except Exception:
        T1 = None
    
%>

%if dest is not None:
<tr>
    <td class="results__name">ΔE<sub>ST</sub>:</td>
    <td class="results__value">${"{:0.2f}".format(dest)} eV</td>
</tr>
%endif
%if S1 is not None:
<%include file="state_result_rows.mako" args="state=S1"/>
%endif
%if T1 is not None:
<%include file="state_result_rows.mako" args="state=T1"/>
%endif
