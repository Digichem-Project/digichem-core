## -*- coding: utf-8 -*-

<%page args="energies, report"/>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", energies.energy_type + " summary")}:</div>
		Summary of ${energies.energy_type} energy properties.
	</div>
    <table class="resultsTable__table">
        <tr class="resultTable__row">
            <td class="resultsTable__title resultsTable__cell">No. of steps</td>
            <td class="resultsTable__value resultsTable__cell">${len(energies)}</td>
        </tr>
        <tr>
            <td class="resultsTable__title resultsTable__cell">Final energy</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.4f}".format(energies.final)} eV</td>
        </tr>
        <tr>
            <td class="resultsTable__title resultsTable__cell">Final energy</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0,.0f}".format(energies.eV_to_kJmol(energies.final))} kJ⋅mol<sup>-1</sup></td>
        </tr>
    </table>
</div>