## -*- coding: utf-8 -*-

<%page args="alignment"/>

<%
	# Build our formula string here becase we don't want any spaces between the letters, which can be tricky to do in the template.
	formula_string = ""
	
	# Go through each element and add to the string.
	for element in alignment.element_dict:
		# Add the element symbol (which is our key).
		formula_string += element
		# Add the number, unless it is 1.
		if alignment.element_dict[element] > 1:
			formula_string += "<sub>{}</sub>".format(alignment.element_dict[element])
			
	# Finally, add on our charge if we have one, but don't include the number if we are +/- 1.
	if alignment.charge == 1:
		formula_string += "<sup>+</sup>"
	elif alignment.charge == -1:
		formula_string += "<sup>-</sup>"
	elif alignment.charge != 0:
		formula_string += "<sup>{:+}</sup>".format(alignment.charge)
%>

<div class="resultsContainer">
	<div class="reportHeader reportHeader--minor reportHeader--results">Geometry</div>
	<table class="results results--et">
		<tr>
			<td class="results__name">Formula:</td>
			<td class="results__value">${formula_string}</td>
		</tr>
		%if alignment.safe_get('mass') is not None:
		<tr>
			<td class="results__name">Exact mass:</td>
			<td class="results__value">${"{:0.4f}".format(alignment.mass)} gmol<sup>-1</sup></td>
		</tr>
		%endif
		<tr>
			<td class="results__name">Molar mass:</td>
			<td class="results__value">${"{:0.4f}".format(alignment.molar_mass)} gmol<sup>-1</sup></td>
		</tr>
		<tr>
			<td class="results__name">Alignment method:</td>
			<td class="results__value">${alignment.CLASS_HANDLE[0]}</td>
		</tr>
		<tr>
			<td class="results__name">X extension:</td>
			<td class="results__value">${"{:0.2f}".format(alignment.X_length)} Å</td>
		</tr>
		<tr>
			<td class="results__name">Y extension:</td>
			<td class="results__value">${"{:0.2f}".format(alignment.Y_length)} Å</td>
		</tr>
		<tr>
			<td class="results__name">Z extension:</td>
			<td class="results__value">${"{:0.2f}".format(alignment.Z_length)} Å</td>
		</tr>
		<tr>
			<td class="results__name">Linearity ratio:</td>
			<td class="results__value">${"{:0.2f}".format(alignment.get_linear_ratio())}</td>
		</tr>
		<tr>
			<td class="results__name">Planarity ratio:</td>
			<td class="results__value">${"{:0.2f}".format(alignment.get_planar_ratio())}</td>
		</tr>
	</table>
</div>