## -*- coding: utf-8 -*-

<%!
	from silico.exception import Result_unavailable_error
	import math
%>

<%page args="result"/>

<%
	# Split our list of orbitals to render images for into two, one for HOMO-n, another for LUMO+n.
	pre_HOMO_orbitals = [orbital for orbital in result.orbitals_to_render if orbital.HOMO_difference < 0]
	post_LUMO_orbitals = [orbital for orbital in result.orbitals_to_render if orbital.HOMO_difference > 1]
%>
<!DOCTYPE html>

<html>
	<head>
		<meta charset="utf-8"/>
		<% 
		stylesheets = [
			'font.css',
			'report.css',
			'table.css',
			'image.css',
			'results.css',
			'front_page.css',
			'metadata.css',
			'geometry.css',
			'mo.css',
			'vibrations.css',
			'excited_states.css',
			'excited_states_table.css',
			'reference.css',
			'summary.css',
			'color_box.css',
			'energies.css',
			'absorptions.css',
			'about.css'
		]
		%>
		%for stylesheet in stylesheets:
		<link rel="stylesheet" type="text/css" href="static/css/${stylesheet}">
		%endfor
	</head>
	<body>
		<%include file="/front_page/front_page.mako" args="result = result"/>
		<%include file="/summary/summary_section.mako" args="result = result"/>
		## We don't need these sections unless we're doing an opt.
		%if len(result.CC_energies) > 1:
			<%include file="/energy/energy_section.mako" args="energies = result.CC_energies"/>
		%endif
		%if len(result.MP_energies) > 1:
			<%include file="/energy/energy_section.mako" args="energies = result.MP_energies"/>
		%endif
		%if len(result.SCF_energies) > 1:
			<%include file="/energy/energy_section.mako" args="energies = result.SCF_energies"/>
		%endif
		%if len(result.atoms) > 0:
			<%include file="/geometry/geometry_section.mako" args="alignment = result.alignment"/>
		%endif
		%if result.dipole_moment is not None:
			<%include file="/dipole_moment/dipole_moment_section.mako" args="dipole_moment = result.dipole_moment"/>
		%endif
		%if result.transition_dipole_moment is not None:
			<%include file="/dipole_moment/dipole_moment_section.mako" args="dipole_moment = result.transition_dipole_moment"/>
		%endif
		%if result.metadata.system_multiplicity != 1:
			<%include file="/spin/spin_density_section.mako" args="result = result"/>
		%endif
		%if len(pre_HOMO_orbitals) > 0:
			<%include file="/orbitals/orbitals_section.mako" args="molecular_orbitals = pre_HOMO_orbitals"/>
		%endif
		%if len(result.molecular_orbitals) > 0:
			<%include file="/orbitals/HOMO_LUMO_section.mako" args="molecular_orbitals = result.molecular_orbitals"/>
		%endif
		%if len(result.beta_orbitals) > 0:
			<%include file="/orbitals/HOMO_LUMO_section.mako" args="molecular_orbitals = result.beta_orbitals"/>
		%endif
		%if len(post_LUMO_orbitals) > 0:
			<%include file="/orbitals/orbitals_section.mako" args="molecular_orbitals = post_LUMO_orbitals"/>
		%endif
		%if result.vertical_emission is not None:
			<%include file="/emission/emission_section.mako" args="relaxed_excited_state = result.vertical_emission"/>
		%endif
		%if result.adiabatic_emission is not None:
			<%include file="/emission/emission_section.mako" args="relaxed_excited_state = result.adiabatic_emission"/>
		%endif
		%if len(result.excited_states) > 0:
			<%include file="/excited_states/excited_states_section.mako" args="excited_states = result.excited_states"/>
		%endif
		%if len(result.spin_orbit_coupling) > 0:
			<%include file="/spin_orbit_coupling/SOC_table.mako" args="spin_orbit_coupling = result.spin_orbit_coupling"/>
		%endif
		%if len(result.vibrations) > 0:
			<%include file="/vibrations/vibrations_section.mako" args="vibrations = result.vibrations" />
			<%include file="/vibrations/vibrations_table.mako" args="vibrations = result.vibrations, min_frequency = result.options['report']['frequency_table']['min_frequency'], max_frequency = result.options['report']['frequency_table']['max_frequency'], max_num = result.options['report']['frequency_table']['max_num']" />
		%endif
		%if len(result.molecular_orbitals) > 0 or len(result.beta_orbitals) > 0:
			##<%include file="/orbitals/select_mo_table.mako" args="molecular_orbitals = result.molecular_orbitals, beta_orbitals = result.beta_orbitals, min_HOMO_difference = -10, max_LUMO_difference = 10 "/>
			<%include file="/orbitals/select_mo_table.mako" args="molecular_orbitals = result.molecular_orbitals, beta_orbitals = result.beta_orbitals, min_HOMO_difference = result.options['report']['orbital_table']['min'], max_HOMO_difference = result.options['report']['orbital_table']['max']"/>
		%endif
		%if len(result.alignment) > 0:
			<%include file="/geometry/atom_list_section.mako" args="atoms = result.alignment"/>
		%endif
		<%include file="/about/about_section.mako"/>
		<%include file="/references/references_section.mako"/>
	</body>
</html>