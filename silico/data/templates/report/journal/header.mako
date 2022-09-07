## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	from silico.exception import Result_unavailable_error
	from silico.misc.text import andjoin, ordinal_suffix, text_integer
	import inflect
%>

<%
    inflector = inflect.engine()

    # Singlet/triplet splitting energy, only available if both singlets and triplets have been calculated.
    try:
        dest = report.result.excited_states.singlet_triplet_energy
    except Exception:
        dest = None
        
    # S1, the lowest energy singlet excited state. We include this in our summary because it is the most important excited state for most structures re. fluorescence.
    try:
        S1 = report.result.excited_states.get_state("S(1)")
    except Exception:
        S1 = None
        
    # T1, the lowest energy singlet excited state. We include this in our summary because it is the most important excited state for most structures re. phosphorescence.
    try:
        T1 = report.result.excited_states.get_state("T(1)")
    except Exception:
        T1 = None
    
%>

<div class="header">
	<div class="header__titleArea">
		<img class="header__logo" src="static/image/logo.png">
		<div class="header__info">
			<div class="header__title h1">
				A report on the calculation of the ${report.result.metadata.human_calculations_string}
				of ${report.result.metadata.molecule_name}
				at the <div class="nocap">${report.result.level_of_theory}</div> level
			</div>
		 	%if report.result.metadata.user is not None:
		 	<div class="header__authors">
		 	${report.result.metadata.user}
		 	</div>
		 	%endif
		 	%if report.result.metadata.date is not None:
		 	<div class="header__date">
			 	%if report.result.metadata.user is not None and report.result.metadata.date is not None:
			 	-
			 	%endif
		 	${report.result.metadata.date.strftime("%d<sup>{}</sup> %B %Y".format(ordinal_suffix(report.result.metadata.date.strftime("%d"))))}
		 	</div>
		 	% endif
	 	</div>
	</div>
	<div class="abstract">
		<div class="abstract__imageArea">
			% if report.options.report['front_page_image'] == 'skeletal' and 'skeletal' in report.images:
		    <img class="abstract__image" src="${report.relative_image('skeletal')}">
		    % elif report.options.report['front_page_image'] == 'rendered' and 'structure' in report.images:
		    <img class="abstract__image" src="${report.relative_image('structure', 'x0y0z0')}">
			% endif
		</div><!--
 		--><div class="abstract__body content">
	 		<h5>Abstract</h5>
	 		The calculation of ${report.result.metadata.human_calculations_string} for the system '${report.result.metadata.molecule_name}' is presented, accompanied by automated analysis and image generation provided by the Silico software package.
	 		% if len(report.result.results) == 1:
	 		## Only one calc.
	 		The calculation was performed using the ${report.result.metadata.package} software package at the ${report.result.level_of_theory} level of theory.
	 		% else:
	 		## Multiple calcs.
	 		The calculations were performed using the ${report.result.metadata.package} software package(s) at the ${report.result.level_of_theory} level of theory.
	 		% endif
			## Now highlight main results.
			##
			## Final energies
			% for energy in report.result.energies:
				% if len(energy) > 0:
				The total ${energy.human_energy_type} (${energy.energy_type}) energy of the system was found to be ${"{:.2f}".format(energy.final)} eV after ${len(energy)} ${inflector.plural("step", len(energy))}.
				% endif
			% endfor
			##
			## Orbitals.
			% try:
				% if len(report.result.beta_orbitals) > 0:
				## Alpha and beta.
				The alpha and beta highest-occupied molecular orbitals (HOMO) were calculated to be ${"{:.2f}".format(report.result.orbitals.HOMO_energy)} and ${"{:.2f}".format(report.result.beta_orbitals.HOMO_energy)} eV respectively, while the alpha and beta lowest-unoccupied molecular orbitals (LUMO) were ${"{:.2f}".format(report.result.orbitals.LUMO_energy)} and ${"{:.2f}".format(report.result.beta_orbitals.LUMO_energy)} eV. These values correspond to a calculated HOMO-LUMO band gap of ${"{:.2f}".format(report.result.orbitals.HOMO_LUMO_energy)} and ${"{:.2f}".format(report.result.beta_orbitals.HOMO_LUMO_energy)} eV for the alpha and beta case respectively. 
				% elif len(report.result.orbitals) > 0:
				## Restricted.
				The highest-occupied molecular orbital (HOMO) and lowest-unoccupied molecular orbital (LUMO) were calculated to be ${"{:.2f}".format(report.result.orbitals.HOMO_energy)} and ${"{:.2f}".format(report.result.orbitals.LUMO_energy)} eV respectively, corresponding to a HOMO-LUMO band gap of ${"{:.2f}".format(report.result.orbitals.HOMO_LUMO_energy)} eV.
				% endif
			##
			% except Result_unavailable_error:
			## Do nothing.
			% endtry
			##
			## PDM.
			% if report.result.dipole_moment is not None:
			The permanent dipole moment (PDM) was calculated to be ${"{:.2f}".format(report.result.dipole_moment.total)} D.
			% endif
			##
			## Vibrations.
			% if len(report.result.vibrations) > 0:
			## TODO: This will fail if there are no vibrations with intensity above 0 (unlikely but probable); need to safe_get_file() first.
			The most intense vibrational frequencies ${inflector.plural("was", len(report.images['simulated_IR_graph'].selected_peaks(0, 5)))} calculated to be at ${andjoin(report.images['simulated_IR_graph'].selected_peaks(0, 5))} cm<sup>-1</sup>,
			and there ${inflector.plural("was", len(report.result.vibrations.negative))} ${text_integer(len(report.result.vibrations.negative))} negative frequencies.
			% endif
			##
			## Excited states.
			% if len(report.result.excited_states) > 0:
				In total, ${len(report.result.excited_states)} excited states ${inflector.plural("was", len(report.result.excited_states))} calculated with ${report.result.excited_states.multiplicity_strings} multiplicity.
				% if 'simulated_absorption_graph' in report.images and report.images['simulated_absorption_graph'].safe_get_file() is not None:
					<%
					abs_peaks = report.images['simulated_absorption_graph'].selected_peaks(0, 5)
					%>
					The most intense absorption ${inflector.plural("peak", len(abs_peaks))} ${inflector.plural("was", len(abs_peaks))} calculated to be at ${andjoin(abs_peaks)} nm.
				% endif
				## S(1) and T(1).
				% if S1 is not None and T1 is not None:
					The lowest energy singlet and triplet excited states (S<sub>1</sub> and T<sub>1</sub>) were calculated to be ${"{:.2f}".format(S1.energy)} and  ${"{:.2f}".format(T1.energy)} eV (${"{:.0f}".format(S1.wavelength)} and ${"{:.0f}".format(T1.wavelength)} nm) respectively, corresponding to a singlet/triplet splitting energy (ΔE<sub>ST</sub>) of ${"{:.2f}".format(dest)} eV.
				% elif S1 is not None:
					The lowest energy singlet excited state (S<sub>1</sub>) was calculated to be ${"{:.2f}".format(S1.energy)} eV (${"{:.0f}".format(S1.wavelength)} nm).
				% elif T1 is not None:
					The lowest energy triplet excited state (S<sub>1</sub>) was calculated to be ${"{:.2f}".format(T1.energy)} eV (${"{:.0f}".format(T1.wavelength)} nm).
				% endif
			% endif
			##
			## Emission.
		</div>
	</div>
</div>