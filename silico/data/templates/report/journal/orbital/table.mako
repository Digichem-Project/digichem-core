## -*- coding: utf-8 -*-

<%!
    import math
    from itertools import zip_longest
    from silico.result.orbital import Molecular_orbital_list
%>

<%page args="report"/>
## We'll only display orbitals between min_HOMO_difference and max_LUMO_difference. These defaults will span HOMO-10 to LUMO+10 (for a total of 12 orbitals)

<%
    molecular_orbitals = report.result.molecular_orbitals
    beta_orbitals = report.result.beta_orbitals
    min_HOMO_difference = report.options['report']['orbital_table']['min']
    max_HOMO_difference = report.options['report']['orbital_table']['max']

    # Set our defaults.
    if min_HOMO_difference is None:
        min_HOMO_difference = -math.inf
    if max_HOMO_difference is None:
        max_HOMO_difference = math.inf
    
    # We need to work out which orbitals we're going to display.
    # This isn't quite as simple as searching by HOMO-LUMO difference if we have both alpha and beta orbitals, because the orbital lists might not align (ie, AHOMO-10 = BHOMO-12).
    # Use handy method to determine our limits.
    min_orbital_level = molecular_orbitals.find_common_level(beta_orbitals, HOMO_difference = min_HOMO_difference)
    max_orbital_level = molecular_orbitals.find_common_level(beta_orbitals, HOMO_difference = max_HOMO_difference)
    
    # The column headers of our table.
    if len(molecular_orbitals) > 0 and len(beta_orbitals) > 0:
        columns = ["Level", "Label", "Symmetry", "Energy /eV", "Label", "Symmetry", "Energy /eV"]
    else:
        columns = ["Level", "Label", "Symmetry", "Energy /eV"]
    
%>
<div class="resultsTable">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "orbitals")}:</div>
		Energies of the calculated molecular orbitals.
	</div>
    <table class="resultsTable__table">
        <tr class="resultsTable__row">
            %for column in columns:
            <td class="resultsTable__title resultsTable__cell">${column}</td>
            %endfor
        </tr>
        ## We can use our min and max orbital levels as index limits to traverse, we just need to be careful because there's no gaurantee these indices are valid for both lists.
        ## We traverse backwards, because orbitals are typically presented in increasing energy.
        %for orbital_index in range(max_orbital_level-1, min_orbital_level-2, -1):
        <%
        # Now try and get our two orbitals.
        try:
            molecular_orbital = molecular_orbitals[orbital_index]
        except IndexError:
            molecular_orbital = None
        try:
            beta_orbital = beta_orbitals[orbital_index]
        except IndexError:
            beta_orbital = None
        %>
        <tr class="resultsTable__row">
            <td class="resultsTable__value resultsTable__cell">
                ${orbital_index +1}
            </td>
            %if molecular_orbital is not None:
                %if molecular_orbital.HOMO_difference == 0 or molecular_orbital.HOMO_difference == 1:
                <td class="resultsTable__value resultsTable__cell resultsTable__cell--fmo">
                %else:
                <td class="resultsTable__value resultsTable__cell">
                %endif
                    ${molecular_orbital.label}
                </td>
                %if molecular_orbital.HOMO_difference == 0 or molecular_orbital.HOMO_difference == 1:
                <td class="resultsTable__value resultsTable__cell resultsTable__cell--fmo">
                %else:
                <td class="resultsTable__value resultsTable__cell">
                %endif
                    ${molecular_orbital.symmetry}
                </td>
                %if molecular_orbital.HOMO_difference == 0 or molecular_orbital.HOMO_difference == 1:
                <td class="resultsTable__value resultsTable__cell resultsTable__cell--fmo">
                %else:
                <td class="resultsTable__value resultsTable__cell">
                %endif
                    ${"{:0.4f}".format(molecular_orbital.energy)}
                </td>
            %endif
            %if beta_orbital is not None:
                %if beta_orbital.HOMO_difference == 0 or beta_orbital.HOMO_difference == 1:
                <td class="resultsTable__value resultsTable__cell resultsTable__cell--fmo">
                %else:
                <td class="resultsTable__value resultsTable__cell">
                %endif
                    ${beta_orbital.label}
                </td>
                %if beta_orbital.HOMO_difference == 0 or beta_orbital.HOMO_difference == 1:
                <td class="resultsTable__value resultsTable__cell resultsTable__cell--fmo">
                %else:
                <td class="resultsTable__value resultsTable__cell">
                %endif
                    ${beta_orbital.symmetry}
                </td>
                %if beta_orbital.HOMO_difference == 0 or beta_orbital.HOMO_difference == 1:
                <td class="resultsTable__value resultsTable__cell resultsTable__cell--fmo">
                %else:
                <td class="resultsTable__value resultsTable__cell">
                %endif
                    ${"{:0.4f}".format(beta_orbital.energy)}
                </td>
            %endif
        </tr>
        %endfor
    </table>
</div>