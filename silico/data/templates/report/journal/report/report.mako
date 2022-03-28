## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	from silico.misc.text import andjoin
%>

<!DOCTYPE html>

<html>
    <head>
        <meta charset="utf-8"/>
        <%
        stylesheets = [
        	"general.css",
        	"results.css",
        	"image.css",
        	"table.css",
        	"header.css",
        	"color_box.css",
        	"citation.css",
        	"reference.css"
        ]
        %>
        %for stylesheet in stylesheets:
        <link rel="stylesheet" type="text/css" href="static/css/${stylesheet}">
        %endfor
    </head>
    <body>
    	<%include file="/header.mako" args="report = report"/>
    	<%include file="/metadata/table.mako" args="report = report"/>
    	<div class="reportBody">
    		## Overall metadata.
    		<%include file="/metadata/metadata.mako" args="metadata = report.result.metadata, report = report, primary = True"/>
    		## Individual metadatas.
    		%if len(report.result.metadatas) > 1:
	    		%for metadata in report.result.metadatas:
	    		<%include file="/metadata/metadata.mako" args="metadata = metadata, report = report, primary = False"/>
	    		%endfor
    		%endif
    		%for energy in [report.result.SCF_energies, report.result.MP_energies, report.result.CC_energies]:
	    		%if len(energy) > 0:
	    		<%include file="/energy.mako" args="energy = energy, report = report"/>
	    		%endif
    		%endfor
    		%if "spin_density" in report.images:
    		<%include file="/spin.mako" args="report = report"/>
    		%endif
	    	%if len(report.result.alignment) > 0:
	    	<%include file="/geometry/geometry.mako" args="report = report"/>
	    	%endif
	    	%if report.result.dipole_moment is not None:
	    	<%include file="/dipole_moment.mako" args="dipole_moment = report.result.dipole_moment, report = report, image_name = 'dipole_moment'"/>
	    	%endif
	    	%if report.result.transition_dipole_moment is not None:
	    	<%include file="/dipole_moment.mako" args="dipole_moment = report.result.transition_dipole_moment, report = report, image_name = '{}_dipole'.format(report.result.transition_dipole_moment.excited_state.state_symbol)"/>
	    	%endif
	    	%if len(report.result.molecular_orbitals) > 0:
	    	<%include file="/orbital/orbitals.mako" args="report = report"/>
	    	%endif
	    	%if len(report.result.vibrations) > 0:
	    	<%include file="/vibrations.mako" args="report = report"/>
	    	%endif
	    	%if len(report.result.excited_states) > 0:
	    	<%include file="/excited_states/excited_states.mako" args="report = report"/>
	    	%endif
	        <%include file="/emission.mako" args="report = report"/>
    	</div>
    	%if len(report.result.excited_states) > 0:
    	<%include file="/excited_states/table.mako" args="report = report"/>
    	%endif
    	%if len(report.result.spin_orbit_coupling) > 0:
    	<%include file="/excited_states/SOC.mako" args="report = report"/>
    	%endif
    </body>
</html>