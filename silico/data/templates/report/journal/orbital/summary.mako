## -*- coding: utf-8 -*-

<%!
    from silico.exception import Result_unavailable_error
%>

<%page args="molecular_orbitals, report"/>

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

<div class="resultsTable resultsTable--summary">
	<div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "MO summary")}:</div>
		Summary of ${title} properties.
	</div>
	<table class="resultsTable__table">
        <tr class="resultTable__row">
            <td class="resultsTable__title resultsTable__cell">E<sub>HOMO,LUMO</sub></td>
            %if HOMO_energy is not None and LUMO_energy is not None:
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(LUMO_energy - HOMO_energy)} eV</td>
            %else:
            <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultTable__row">
            <td class="resultsTable__title resultsTable__cell">E<sub>HOMO</sub></td>
            %if HOMO_energy is not None:
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(HOMO_energy)} eV</td>
            %else:
            <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
        <tr class="resultTable__row">
            <td class="resultsTable__title resultsTable__cell">E<sub>LUMO</sub></td>
            %if LUMO_energy is not None:
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(LUMO_energy)} eV</td>
            %else:
            <td class="resultsTable__value resultsTable__cell">N/A</td>
            %endif
        </tr>
    </table>
</div>