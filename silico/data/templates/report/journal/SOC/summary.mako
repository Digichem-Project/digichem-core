## -*- coding: utf-8 -*-

<%!
    import math
%>

<%page args="spin_orbit_coupling"/>

<%
    soc_list = []
    for state1, state2 in [("S(0)", "T(1)"), ("S(1)", "T(1)")]:
        # Get parameters.
        try:
            soc_list.append(spin_orbit_coupling.between(state1, state2))
        except Exception:
            # Skip.
            pass
%>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "SOC summary")}:</div>
		Summary of the calculated spin-orbit coupling values.
		&lt;S|H<sub>SO</sub>|T&gt;: SOC between singlet state S and triplet state T.
		&lt;S|λ|T&gt;: First-order mixing coefficient between the same.
	</div>
    <table class="resultsTable__table">
        %for soc in soc_list:
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">&lt;${soc.singlet_state.multiplicity_symbol}<sub>${soc.singlet_state.multiplicity_level}</sub>|H<sub>SO</sub>|${soc.triplet_state.multiplicity_symbol}<sub>${soc.triplet_state.multiplicity_level}</sub>&gt;</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(soc.wavenumbers)} cm<sup>-1</sup></td>
        </tr>
        ## Also coupling constant.
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">&lt;${soc.singlet_state.multiplicity_symbol}<sub>${soc.singlet_state.multiplicity_level}</sub>|λ|${soc.triplet_state.multiplicity_symbol}<sub>${soc.triplet_state.multiplicity_level}</sub>&gt;</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(soc.mixing_coefficient) if soc.mixing_coefficient != math.inf else "∞"}</td>
        </tr>
        %endfor
    </table>
</div>
