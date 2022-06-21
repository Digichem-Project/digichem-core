## -*- coding: utf-8 -*-

<%page args="energies"/>

<div class="resultsContainer">
    <div class="reportHeader reportHeader--minor reportHeader--results">${energies.energy_type} Energies</div>
    <table class="results results--et">
        <tr>
            <td class="results__name">No. of steps:</td>
            <td class="results__value">${len(energies)}</td>
        </tr>
        <tr>
            <td class="results__name">Final energy:</td>
            <td class="results__value">${"{:0.4f}".format(energies.final)} eV</td>
        </tr>
        <tr>
            <td class="results__name">Final energy:</td>
            <td class="results__value">${"{:0,.0f}".format(energies.eV_to_kJmol(energies.final))} kJmol<sup>-1</sup></td>
        </tr>
    </table>
</div>