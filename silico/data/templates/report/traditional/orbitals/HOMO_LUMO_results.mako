## -*- coding: utf-8 -*-

<%!
    from silico.exception import Result_unavailable_error
%>

<%page args="molecular_orbitals"/>

<%
    # First work out the type of orbitals we are dealing with (alpha, beta, or restricted).
    spin_type = molecular_orbitals.spin_type
    
    # The base of our title.
    title = "HOMO & LUMO"
    
    # Add a bit more if we've got some funky orbitals.
    if spin_type != "none":
        title = "{} ({})".format(title, spin_type)
        
    # Try and get our FMOs.
    try:
        HOMO_energy = molecular_orbitals.get_orbital(HOMO_difference = 0).energy
    except Result_unavailable_error:
        HOMO_energy = None
        
    try:
        LUMO_energy = molecular_orbitals.get_orbital(HOMO_difference = 1).energy
    except Result_unavailable_error:
        LUMO_energy = None
%>

<div class="resultsContainer">
    <div class="reportHeader reportHeader--minor reportHeader--results">${title}</div>
    <table class="results">
        <tr>
            <td class="results__name">E<sub>HOMO/LUMO</sub>:</td>
            %if HOMO_energy is not None and LUMO_energy is not None:
            <td class="results__value">${"{:0.2f}".format(LUMO_energy - HOMO_energy)} eV</td>
            %else:
            <td class="results__value results__value--bad">None</td>
            %endif
        </tr>
        <tr>
            <td class="results__name">E<sub>HOMO</sub>:</td>
            %if HOMO_energy is not None:
            <td class="results__value">${"{:0.2f}".format(HOMO_energy)} eV</td>
            %else:
            <td class="results__value results__value--bad">None</td>
            %endif
        </tr>
        <tr>
            <td class="results__name">E<sub>LUMO</sub>:</td>
            %if LUMO_energy is not None:
            <td class="results__value">${"{:0.2f}".format(LUMO_energy)} eV</td>
            %else:
            <td class="results__value results__value--bad">None</td>
            %endif
        </tr>
    </table>
</div>