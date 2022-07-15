## -*- coding: utf-8 -*-

<%page args="energy, report"/>

<%!
	import inflect
	from silico.misc.text import text_integer
%>

<%
	inflector = inflect.engine()
%>

<div class="content">
	<h5>Total <div class="nocap">${energy.energy_type}</div> Energy</h5>
	The total energy of the system was calculated at the <div class="result"><div class="result__title">${energy.human_energy_type} (${energy.energy_type})</div> level${"," if energy.energy_type == "SCF" else ""}
	%if energy.energy_type == "SCF":
		%if "DFT" in report.result.metadata.methods:
		corresponding to the energy calculated by the density-functional theory (DFT) method,
		%else:
		corresponding to the energy calculated by the Hartree-Fock (HF) method,
		%endif
	%endif
	%if len(energy) > 1:
	over a total of ${text_integer(len(energy))} steps, the results of which are displayed in figure ${report.captions("figure", '{}_convergence_graph'.format(energy.energy_type))}. The energy calculated by the final step was <div class="result__value">${"{:.2f}".format(energy.final)} eV</div>, corresponding to <div class="result__value">${"{:0,.0f}".format(energy.eV_to_kJmol(energy.final))} KJmol</div></div><sup>-1</sup>.
	%else:
	with a value of <div class="result__value">${"{:.2f}".format(energy.final)} eV</div>, corresponding to <div class="result__value">${"{:0,.0f}".format(energy.eV_to_kJmol(energy.final))} KJmol</div></div><sup>-1</sup>.
	%endif
	%if energy.energy_type in report.images:
	A plot of the total ${energy.energy_type} electron density is shown in figure ${report.captions("figure", energy.energy_type)}.
	%endif 
	%if len(energy) > 1 and '{}_convergence_graph'.format(energy.energy_type) in report.images:
	<div class="resultImage resultImage--graph">
		<img class="resultImage__image" src="${report.relative_image('{}_convergence_graph'.format(energy.energy_type))}">
		<div class="resultImage__caption"><div class="caption">
			Figure ${report.captions("figure", '{}_convergence_graph'.format(energy.energy_type))}:</div>
			Graph of calculated energies at the ${energy.human_energy_type} (${energy.energy_type}) level.
		</div>
	</div>
	%endif
	%if energy.energy_type in report.images:
	<%include file="/geometry/image.mako" args="image_name = energy.energy_type, caption = 'Plot of the total {} electron density, plotted with an isovalue of {}'.format(energy.energy_type, report.images[energy.energy_type].isovalue), report = report" />
	%endif
</div>