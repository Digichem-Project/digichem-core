## -*- coding: utf-8 -*-

## A template for displaying summary table rows of a particular excited state.

<%page args="state" />

<tr>
    <td class="results__name">${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub> energy:</td>
    <td class="results__value">${"{:0.2f}".format(state.energy)} eV</td>
</tr>
<tr>
    <td class="results__name">${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub> wavelength:</td>
    <td class="results__value">${"{:0.0f} nm".format(state.wavelength) if state.safe_get('wavelength') is not None else "N/A"}</td>
</tr>
<tr>
    <td class="results__name">${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub> colour:</td>
    <td class="results__value results__value--color">${state.color if state.safe_get('color') is not None else "N/A"}
    %if state.safe_get('color') is not None:
    <%include file="color_box.mako" args="rgb = state.rgb"/></td>
    %endif
</tr>
<tr>
    <td class="results__name">${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub> CIE (x,y):</td>
    <td class="results__value">${"({:0.2f}, {:0.2f})".format(state.CIE_xy[0], state.CIE_xy[1]) if state.safe_get('CIE_xy') is not None else "N/A"}</td>
</tr>
<tr>
    <td class="results__name">${state.multiplicity_symbol}<sub>${state.multiplicity_level}</sub> oscillator strength:</td>
    <td class="results__value">${"{:0.2f}".format(state.oscillator_strength) if state.safe_get('oscillator_strength') is not None else "N/A"}</td>
</tr>