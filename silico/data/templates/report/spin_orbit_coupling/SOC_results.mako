## -*- coding: utf-8 -*-

<%!
    import math
%>

<%page args="spin_orbit_coupling, excited_states"/>

<%
    soc_list = []
    for state1, state2 in [("S(0)", "T(1)"), ("S(1)", "T(1)")]:
        # Get parameters.
        try:
            soc_list.append(spin_orbit_coupling.between(state1, state2))
        except Exception:
            # Skip.
            raise
%>

<div class="resultsContainer">
    <div class="reportHeader reportHeader--minor reportHeader--results">Spin-Orbit Coupling</div>
    <table class="results">
        %for soc in soc_list:
        <tr>
            <td class="results__name">&lt;${soc.singlet_state.multiplicity_symbol}<sub>${soc.singlet_state.multiplicity_level}</sub>|H<sub>SO</sub>|${soc.triplet_state.multiplicity_symbol}<sub>${soc.triplet_state.multiplicity_level}</sub>&gt;:</td>
            <td class="results__value">${"{:0.2f}".format(soc.wavenumbers)} cm<sup>-1</sup></td>
        </tr>
        ## Also coupling constant.
        <tr>
            <td class="results__name">&lt;${soc.singlet_state.multiplicity_symbol}<sub>${soc.singlet_state.multiplicity_level}</sub>|λ|${soc.triplet_state.multiplicity_symbol}<sub>${soc.triplet_state.multiplicity_level}</sub>&gt;:</td>
            <td class="results__value">${"{:0.2f}".format(soc.mixing_coefficient) if soc.mixing_coefficient != math.inf else "∞"}</td>
        </tr>
        %endfor
    </table>
</div>
