## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	import inflect
	from silico.misc.text import andjoin, text_integer
%>

<%
	inflector = inflect.engine()
%>

<%def name="emission_section(emissions, emissions_type)">
	<%
		state_names = ["{}<sub>{}</sub>".format(emission.multiplicity_symbol, emission.multiplicity_level) for emission in emissions.values()]
	%>
	The ${emissions_type} emission energy,
	corresponding to the difference in energy between an excited state at an excited state geometry, and the ground state at the ground state geometry,
	from the <div class="result"><div class="result__title">${andjoin(state_names)} ${inflector.plural("state", len(state_names))}</div> to the ground state was calculated and found to be
	<div class="result__value">${andjoin(["{:.2f}".format(emission.energy) for emission in emissions.values()])} eV</div></div>${" respectively." if len(state_names) > 1 else "."}
	%for emission in emissions.values():
		%if len(state_names) > 1:
			For the ${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> emission, this
		%else:
			This
		%endif
		energy is equivalent to emission of a photon with a wavelength of ${"{:.0f}".format(emission.wavelength)} nm, corresponding to a colour of ${emission.color} <%include file="/excited_states/color.mako" args="colorRGB = emission.rgb"/> and CIE coordinates (x,y) of ${"({}, {})".format(*emission.CIE_xy)}.
		The excited state had a total energy of ${"{:.2f}".format(emission.excited_energy)} eV and a multiplicity of ${text_integer(emission.excited_multiplicity)}, while the ground state had a total energy of ${"{:.2f}".format(emission.ground_energy)} eV and a multiplicity of ${text_integer(emission.ground_multiplicity)}.
		This emission is therefore a ${emission.emission_type} type process, because
		%if emission.emission_type == "fluorescence":
			both the ground and excited state have the same multiplicity.
		%else:
			the ground and excited state have different multiplicity.
		%endif
	%endfor
</%def> 

%if len(report.result.adiabatic_emission) != 0 or len(report.result.vertical_emission):
<div class="content">
	<h5>Emission energy</h5>
	%if len(report.result.adiabatic_emission) > 0:
		## Adiabatic.
		${emission_section(report.result.adiabatic_emission, "adiabatic")}
	%endif
	%if len(report.result.vertical_emission) > 0:
		## Adiabatic.
		${emission_section(report.result.vertical_emission, "vertical")}
	%endif
</div>
%endif