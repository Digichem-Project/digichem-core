## -*- coding: utf-8 -*-

<%page args="alignment, report"/>

<div class="resultsTable resultsTable--summary">
    <div class="resultsTable__caption">
		<div class="caption">Table ${report.captions("table", "geometry summary")}:</div>
		Summary of geometry properties.
	</div>
    <table class="resultsTable__table">
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Formula</td>
            <td class="resultsTable__value resultsTable__cell">
            	<%include file="/geometry/formula.mako" args="atoms = alignment"/>
            </td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">SMILES</td>
            <td class="resultsTable__value resultsTable__cell">${alignment.smiles}</td>
        </tr>
        %if alignment.safe_get('mass') is not None:
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Exact mass</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.4f}".format(alignment.mass)} g⋅mol<sup>-1</sup></td>
        </tr>
        %endif
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Molar mass</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.4f}".format(alignment.molar_mass)} g⋅mol<sup>-1</sup></td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Alignment method</td>
            <td class="resultsTable__value resultsTable__cell">${alignment.CLASS_HANDLE[0]}</td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">X extension</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(alignment.X_length)} Å</td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Y extension</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(alignment.Y_length)} Å</td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Z extension</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(alignment.Z_length)} Å</td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Linearity ratio</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(alignment.get_linear_ratio())}</td>
        </tr>
        <tr class="resultsTable__row">
            <td class="resultsTable__title resultsTable__cell">Planarity ratio</td>
            <td class="resultsTable__value resultsTable__cell">${"{:0.2f}".format(alignment.get_planar_ratio())}</td>
        </tr>
    </table>
</div>