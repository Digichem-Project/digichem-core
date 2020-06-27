<%page args="atoms"/>

<div class="section section--fullPage">
	<h2 class="section__header">Table of Atoms</h2>
	<div class="section__body section__body--table">
		<div class="tableBorder">
			<table class="table">
				<tr class="table__row table__row--header">
					%for column in ["Element", "Isotope Mass", "X Coord", "Y Coord", "Z Coord"]:
					<th class="table__header">${column}</th>
					%endfor
				</tr>
				%for atom in atoms:
				<tr class="table__row">
					<td class="table__cell">
						${atom.element.symbol}
					</td>
					<td class="table__cell">
						${atom.mass if atom.safe_get('mass') is not None else "N/A"}
					</td>
					<td class="table__cell">
						${"{:0.7f}".format(atom.coords[0])}
					</td>
					<td class="table__cell">
						${"{:0.7f}".format(atom.coords[1])}
					</td>
					<td class="table__cell">
						${"{:0.7f}".format(atom.coords[2])}
					</td>		
				</tr>
				%endfor
			</table>
		</div>
	</div>
</div>