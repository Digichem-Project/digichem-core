## -*- coding: utf-8 -*-

<%!
    from silico.exception import Result_unavailable_error
    import math
%>

<%page args="report"/>

<%
    # Split our list of orbitals to render images for into two, one for HOMO-n, another for LUMO+n.
    pre_HOMO_orbitals = [orbital for orbital in report.orbitals_to_render if orbital.HOMO_difference < 0]
    post_LUMO_orbitals = [orbital for orbital in report.orbitals_to_render if orbital.HOMO_difference > 1]
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
            'about.css',
            'general_warning.css'
        ]
        %>
        %for stylesheet in stylesheets:
        <link rel="stylesheet" type="text/css" href="static/css/${stylesheet}">
        %endfor
    </head>
    <body>
        <%include file="/front_page/front_page.mako" args="report = report"/>
        <%include file="/warning/warning_section.mako" args="result = report.result"/>
        <%include file="/summary/summary_section.mako" args="result = report.result"/>
        ## We don't need these sections unless we're doing an opt.
        %if len(report.result.energies.scf) > 1:
            <%include file="/energy/energy_section.mako" args="energies = report.result.energies.scf, report = report"/>
        %endif
        %if len(report.result.energies.mp) > 1:
            <%include file="/energy/energy_section.mako" args="energies = report.result.energies.mp, report = report"/>
        %endif
        %if len(report.result.energies.cc) > 1:
            <%include file="/energy/energy_section.mako" args="energies = report.result.energies.cc, report = report"/>
        %endif
        %if len(report.result.atoms) > 0:
            <%include file="/geometry/geometry_section.mako" args="report = report"/>
        %endif
        %if 'SCF' in report.images:
        	<%include file="/total_density/section.mako" args="report = report, density_image_name = 'SCF'"/>
        %endif
        %if report.result.pdm is not None:
            <%include file="/dipole_moment/dipole_moment_section.mako" args="dipole_moment = report.result.pdm, report = report, image_name = 'dipole_moment'"/>
        %endif
        %if report.result.transition_dipole_moment is not None:
            <%include file="/dipole_moment/dipole_moment_section.mako" args="dipole_moment = report.result.transition_dipole_moment, report = report, image_name = '{}_dipole'.format(report.result.transition_dipole_moment.excited_state.state_symbol)"/>
        %endif
        %if "spin_density" in report.images:
            <%include file="/spin_density/section.mako" args="result = report.result, report = report"/>
        %endif
        %if len(pre_HOMO_orbitals) > 0:
            <%include file="/orbitals/section.mako" args="orbitals = pre_HOMO_orbitals, report = report"/>
        %endif
        %if len(report.result.orbitals) > 0:
            <%include file="/orbitals/HOMO_LUMO.mako" args="orbitals = report.result.orbitals, report = report"/>
        %endif
        %if len(report.result.beta_orbitals) > 0:
            <%include file="/orbitals/HOMO_LUMO.mako" args="orbitals = report.result.beta_orbitals, report = report"/>
        %endif
        %if len(post_LUMO_orbitals) > 0:
            <%include file="/orbitals/section.mako" args="orbitals = post_LUMO_orbitals, report = report"/>
        %endif
        %for multiplicity, vertical in report.result.emission.vertical.items():
        	<%include file="/emission/emission_section.mako" args="emission = vertical, report = report"/>
        %endfor
        %for multiplicity, adiabatic in report.result.emission.adiabatic.items():
        	<%include file="/emission/emission_section.mako" args="emission = adiabatic, report = report"/>
        %endfor
        %if len(report.result.excited_states) > 0:
            <%include file="/excited_states/excited_states_section.mako" args="excited_states = report.result.excited_states, report = report"/>
        %endif
        <%include file="/difference_density/section.mako" args="excited_states = report.result.excited_states, report = report"/>
        <%include file="/NTO/section.mako" args="excited_states = report.result.excited_states, report = report"/>
        %if len(report.result.soc) > 0:
            <%include file="/soc/SOC_table.mako" args="soc = report.result.soc"/>
        %endif
        %if len(report.result.vibrations) > 0:
            <%include file="/vibrations/vibrations_section.mako" args="vibrations = report.result.vibrations, report = report" />
            <%include file="/vibrations/vibrations_table.mako" args="vibrations = report.result.vibrations, min_frequency = report.options['report']['frequency_table']['min_frequency'], max_frequency = report.options['report']['frequency_table']['max_frequency'], max_num = report.options['report']['frequency_table']['max_num']" />
        %endif
        %if len(report.result.orbitals) > 0 or len(report.result.beta_orbitals) > 0:
            <%include file="/orbitals/table.mako" args="orbitals = report.result.orbitals, beta_orbitals = report.result.beta_orbitals, min_HOMO_difference = report.options['report']['orbital_table']['min'], max_HOMO_difference = report.options['report']['orbital_table']['max']"/>
        %endif
        %if len(report.result.alignment) > 0:
            <%include file="/geometry/atom_list_section.mako" args="atoms = report.result.alignment"/>
        %endif
        <%include file="/about/about_section.mako"/>
        <%include file="/references/references_section.mako"/>
    </body>
</html>