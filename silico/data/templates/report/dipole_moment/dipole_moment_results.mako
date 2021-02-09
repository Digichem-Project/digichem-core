## -*- coding: utf-8 -*-

<%page args="dipole_moment"/>

<%
    # First work out our title, which changes slightly depending on whether this is the ground or excited state dipole.
    if dipole_moment.dipole_type == "permanent":
        # This is the ground state dipole.
        dipole_title = 'Permanent Dipole Moment'
    else:
        # This is an excited states dipole.
        dipole_title = 'Transition ({}<sub>{}</sub>) Dipole Moment'.format(dipole_moment.excited_state.multiplicity_symbol, dipole_moment.excited_state.multiplicity_level) 
%>

<div class="resultsContainer">
    <div class="reportHeader reportHeader--minor reportHeader--results">${dipole_title}</div>
    <table class="results">
        <%include file="dipole_moment_main_results.mako" args="dipole_moment = dipole_moment"/>
    </table>
</div>