## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	import inflect
	from silico.misc.text import text_integer, andjoin, text_join, listjoin, text_float
%>

<%
	inflector = inflect.engine()
	mult_pairs, mults = report.result.excited_states.group_pairs()
	
	# A list of excited states that have NTOs for.
	NTOs = []
	# A list of excited states that have difference density plots.
	diff_dens = []
	for excited_state in report.result.excited_states:
		if excited_state.state_symbol + "_NTO" in report.images:
			NTOs.append(excited_state)
		if excited_state.state_symbol + "_difference_density" in report.images:
			diff_dens.append(excited_state)
		
	
	
%>

<%def name="excited_state_text(state)">
<div class="result__value">${"{:.2f}".format(state.energy)} eV (${"{:.0f}".format(state.wavelength)} nm, ${state.color.lower()} <%include file="/excited_states/color.mako" args="colorRGB = state.rgb"/>, CIE: ${"({:.2f}, {:.2f})".format(state.CIE_xy[0], state.CIE_xy[1])})</div>\
</%def> 

<div class="content">
	<h5>Excited States</h5>
	%if len(mults) == 1:
		In total, the energies of ${text_integer(len(report.result.excited_states))} ${report.result.excited_states[0].multiplicity_string} electronic excited states were calculated, which are shown in figure ${report.captions("figure", "excited_state_energies")}.
	%else:
		In total, the energies of ${text_integer(len(report.result.excited_states))} electronic excited states were calculated (figure ${report.captions("figure", "excited_state_energies")}), consisting of
		%for index, (mult, grouped_states) in enumerate(mults.items()):
			%if index == 0:
				${text_integer(len(grouped_states))} states with a multiplicity of ${grouped_states[0].multiplicity_string}${"," if index < len(mults)-2 else ""}
			%else:
				${text_integer(len(grouped_states))} of multiplicity ${grouped_states[0].multiplicity_string}${"," if index < len(mults)-2 else "." if index == len(mults)-1 else ""}
			%endif
			%if index == len(mults)-2:
			and
			%endif
		%endfor
		
	%endif
	## Specific energies.
	%for index, (mult, grouped_states) in enumerate(mults.items()):
		%if index == 0:
			The energy of the lowest <div class="result"><div class="result__title">${grouped_states[0].multiplicity_string} excited state (${grouped_states[0].multiplicity_symbol}<sub>${grouped_states[0].multiplicity_level}</sub>)</div> was <div class="result__value">${"{:.2f}".format(grouped_states[0].energy)} eV</div></div>, corresponding to absorption by a photon with a wavelength of ${"{:.0f}".format(grouped_states[0].wavelength)} nm, ${inflector.a(grouped_states[0].color.lower())} 'color' <%include file="/excited_states/color.mako" args="colorRGB = grouped_states[0].rgb"/> and CIE coordinates of ${"({:.2f}, {:.2f})".format(grouped_states[0].CIE_xy[0], grouped_states[0].CIE_xy[1])}${text_join(index, len(mults), "while")}
		%elif index == 1:
			the energy of the <div class="result"><div class="result__title">${grouped_states[0].multiplicity_symbol}<sub>${grouped_states[0].multiplicity_level}</sub></div> was ${excited_state_text(grouped_states[0])}</div>${text_join(index, len(mults))}
		%else:
			the <div class="result"><div class="result__title">${grouped_states[0].multiplicity_symbol}<sub>${grouped_states[0].multiplicity_level}</sub> was ${excited_state_text(grouped_states[0])}</div>${text_join(index, len(mults))}
		%endif
	%endfor
	## Difference in energies between different mults.
	%for index, mult_pair in enumerate(mult_pairs):
		<%
		state1 = mults[mult_pair[0]][0]
		state2 = mults[mult_pair[1]][0]
		%>
		%if index == 0:
			The difference in energy between the ${state1.multiplicity_symbol}<sub>${state1.multiplicity_level}</sub> and ${state2.multiplicity_symbol}<sub>${state2.multiplicity_level}</sub> excited states <div class="result"><div class="result__title">(ΔE<sub>${state1.multiplicity_symbol}${state2.multiplicity_symbol}</sub>)</div> was therefore <div class="result__value">${"{:.2f}".format(state1.energy - state2.energy)} eV</div></div>${text_join(index, len(mult_pairs))}
		%else:
			ΔE<sub>${state1.multiplicity_symbol}${state2.multiplicity_symbol}</sub>  was ${text_float(state1.energy - state2.energy)} eV${text_join(index, len(mult_pairs))} 
		%endif
	%endfor
	A complete table of the calculated excited state properties is available in table ${report.captions("table", "excited_states")}.
	## This can be none if all oscillator strengths are zero (eg, triplets only).
	%if 'simulated_absorption_graph' in report.images and report.images['simulated_absorption_graph'].safe_get_file() is not None:
	In addition, an electronic transition spectrum was simulated using a gaussian function with full-width at half maximum (FWHM) of ${"{:.2f}".format(report.images['simulated_absorption_graph'].fwhm)} eV, from which the <div class="result"><div class="result__title">${text_integer(len(report.images['simulated_absorption_graph'].selected_peaks(0, 5)))} most intense ${inflector.plural("peak", len(report.images['simulated_absorption_graph'].selected_peaks(0, 5)))}</div> ${inflector.plural("was", len(report.images['simulated_absorption_graph'].selected_peaks(0, 5)))} found at <div class="result__value">${andjoin(report.images['simulated_absorption_graph'].selected_peaks(0, 5))} nm</div></div>.
	The full simulated absorption spectrum is shown in figure ${report.captions("figure", "simulated_absorption_graph")}.
	%endif
	## NTOs.
	%if len(NTOs) > 0:
		%if len(diff_dens) == 0:
			Finally, 
		%else:
			Also,
		%endif
		<div class="result"><div class="result__title">natural transition orbitals (NTOs)</div> were calculated for
		%if len(NTOs) == len(report.result.excited_states):
			each excited state
		%else:
			the ${andjoin(["{}<sub>{}</sub>".format(state.multiplicity_symbol, state.multiplicity_level) for state in NTOs])} excited states
		%endif
		<%
			nto_captions = []
			for nto_state in NTOs:
				nto_captions.append(report.captions("figure", nto_state.state_symbol + "_NTO"))
		%>
		and are shown in <div class="result__value">${inflector.plural("figure", len(nto_captions))} ${listjoin(nto_captions)}</div></div>.
	%endif
	## Differential density.
	%if len(diff_dens) > 0:
		Finally, the <div class="result"><div class="result__title">difference in density</div> between
		%if len(diff_dens) == len(report.result.excited_states):
			each excited state
		%else:
			the ${andjoin(["{}<sub>{}</sub>".format(state.multiplicity_symbol, state.multiplicity_level) for state in diff_dens])} excited states
		%endif
		and the ground was calculated
		<%
			captions = []
			for state in diff_dens:
				captions.append(report.captions("figure", state.state_symbol + "_difference_density"))
		%>
		and is shown in <div class="result__value">${inflector.plural("figure", len(captions))} ${listjoin(captions)}</div></div>.
	%endif
	<div class="resultImage resultImage--graph">
		<img class="resultImage__image" src="${report.relative_image('excited_state_energies')}">
		<div class="resultImage__caption"><div class="caption">Figure ${report.captions("figure", 'excited_state_energies')}:</div> Graph of the calculated excited states. f: oscillator strength of the relevant ground to excited state transition.</div>
	</div>
	%if 'simulated_absorption_graph' in report.images and report.images['simulated_absorption_graph'].safe_get_file() is not None:
	<div class="resultImage resultImage--graph">
		<img class="resultImage__image" src="${report.relative_image('simulated_absorption_graph')}">
		<div class="resultImage__caption">
			<div class="caption">Figure ${report.captions("figure", 'simulated_absorption_graph')}:</div>
			Graph of the simulated absorption spectrum.
			Excited states are shown as vertical black lines, while peaks simulated with a gaussian function with FWHM: ${"{:.2f}".format(report.images['simulated_absorption_graph'].fwhm)} eV are shown as a blue line.
			%if report.images['simulated_absorption_graph'].hide_y:
				The oscillator strength of each excited state has arbitrarily been set to 1 because all oscillator strengths were 0.
			%endif
			Peaks can be found at: ${andjoin(report.images['simulated_absorption_graph'].selected_peaks())} nm.
		</div>
	</div>
	%endif
</div>
## Spin-orbit coupling.
%if len(report.result.spin_orbit_coupling) > 0:
	<div class="content">
		<h5>Spin-Orbit Coupling</h5>
		## TODO: This section is fine for now because PySOC is the only method implemented for calculating SOC, but in the future when further methods have been added this will need expanding.
		The <div class="result"><div class="result__title">spin-orbit coupling</div>  between each singlet state
		%if report.result.ground_state.multiplicity == 1:
			(including the ground state)
		%endif
		and each triplet excited state was then calculated using a custom implementation of the PySOC program,
		the results of which are displayed in <div class="result__value">table ${report.captions("table", "SOC")}</div></div>.
		<%
			S0_T1 = report.result.spin_orbit_coupling.between("S(0)", "T(1)", default = None)
			S1_T1 = report.result.spin_orbit_coupling.between("S(1)", "T(1)", default = None)
		%>
		%if S0_T1 is not None or S1_T1 is not None:
			From this analysis, the H<sub>SO</sub> between
			%if S0_T1 is not None:
				the <div class="result"><div class="result__title">S<sub>0</sub> and T<sub>1</sub></div> states was found to be <div class="result__value">${text_float(S0_T1.wavenumbers)} cm <sup>-1</sup></div></div>, while the H<sub>SO</sub> between
			%endif
			%if S1_T1 is not None:
				the <div class="result"><div class="result__title">S<sub>1</sub> and T<sub>1</sub></div> excited states was <div class="result__value">${text_float(S1_T1.wavenumbers)} cm <sup>-1</sup></div></div>.
			%endif
			%if S0_T1 is not None and S1_T1 is not None:
				These values correspond
			%else:
				This value corresponds
			%endif
			 to a first-order mixing coefficient (λ = H<sub>SO</sub>/ΔE<sub>ST</sub>) of 
			%if S0_T1 is not None and S1_T1 is not None:
				${text_float(S0_T1.mixing_coefficient)} and ${text_float(S1_T1.mixing_coefficient)} respectively.
			%elif S0_T1 is not None:
				${text_float(S0_T1.mixing_coefficient)}.
			%else:
				${text_float(S1_T1.mixing_coefficient)}.
			%endif
		%endif
	</div>
%endif
<div class="content">
	%for state in NTOs:
	<%include file="/geometry/image.mako" args="image_name = state.state_symbol + '_NTO', caption = 'Density plot of the NTO hole (' + report.images[state.state_symbol + '_NTO'].primary_colour + ') & electron (' + report.images[state.state_symbol + '_NTO'].secondary_colour + ') of the {}<sub>{}</sub> state, plotted with isovalue: {}'.format(state.multiplicity_symbol, state.multiplicity_level, report.images[state.state_symbol + '_NTO'].isovalue), report = report" />
	%endfor
	%for state in diff_dens:
	<%include file="/geometry/image.mako" args="image_name = state.state_symbol + '_difference_density', caption = 'Differential density plot of the hole (' + report.images[state.state_symbol + '_difference_density'].primary_colour + ') & electron (' + report.images[state.state_symbol + '_difference_density'].secondary_colour + ') of the {}<sub>{}</sub> state, plotted with isovalue: {}'.format(state.multiplicity_symbol, state.multiplicity_level, report.images[state.state_symbol + '_difference_density'].isovalue), report = report" />
	%endfor
</div>