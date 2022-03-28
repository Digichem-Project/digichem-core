## -*- coding: utf-8 -*-

<%page args="report"/>

<%!
	from silico.misc.text import andjoin
	import itertools
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
    	<h4>Summary Of Results</h4>
    	<div class="reportBody">
    		## Summaries.
    		%for energies in [report.result.SCF_energies, report.result.MP_energies, report.result.CC_energies]:
    			%if len(energies) > 0:
    				<h5>${energies.energy_type} Energy</h5>
    				<%include file="/energy/summary.mako" args="energies = energies, report = report"/>
    			%endif
    		%endfor
    		%if len(report.result.alignment) > 0:
    			<h5>Geometry</h5>
    			<%include file="/geometry/summary.mako" args="alignment = report.result.alignment, report = report"/>
    		%endif
    		%if len(report.result.molecular_orbitals) > 0:
    			<h5>${"Molecular Orbitals" if report.result.molecular_orbitals.spin_type == 'none' else "Alpha Orbitals"}</h5>
    			<%include file="/orbital/summary.mako" args="molecular_orbitals = report.result.molecular_orbitals, report = report"/>
    		%endif
    		%if len(report.result.beta_orbitals) > 0:
    			<h5>Beta Orbitals</h5>
    			<%include file="/orbital/summary.mako" args="molecular_orbitals = report.result.beta_orbitals, report = report"/>
    		%endif
    		%if report.result.dipole_moment is not None:
    			<h5>Permanent Dipole Moment</h5>
    			<%include file="/dipole_moment/summary.mako" args="dipole_moment = report.result.dipole_moment, report = report"/>
    		%endif
    		%if report.result.transition_dipole_moment is not None:
    			<h5>${report.result.transition_dipole_moment.excited_state.multiplicity_symbol}<sub>${report.result.transition_dipole_moment.excited_state.multiplicity_level}</sub> Transition Dipole Moment</h5>
    			<%include file="/dipole_moment/summary.mako" args="dipole_moment = report.result.transition_dipole_moment, report = report"/>
    		%endif
    		%if len(report.result.vibrations) > 0:
    			<h5>Vibrations</h5>
    			<%include file="/vibrations/summary.mako" args="vibrations = report.result.vibrations, report = report"/>
    		%endif
    		%if len(report.result.excited_states) > 0:
    			<h5>Excited States</h5>
    			<%include file="/excited_states/summary.mako" args="excited_states = report.result.excited_states, report = report"/>
    		%endif
    		%for emission in itertools.chain(report.result.adiabatic_emission.values(), report.result.vertical_emission.values()):
    			<h5>${emission.transition_type} ${emission.multiplicity_symbol}<sub>${emission.multiplicity_level}</sub> Emission</h5>
    			<%include file="/emission/summary.mako" args="emission = emission"/>
    		%endfor
    		##
    		##
    		<h4>Methodology</h4>
    		## Overall metadata.
    		<%include file="/metadata/metadata.mako" args="metadata = report.result.metadata, report = report, primary = True"/>
    		## Individual metadatas.
    		%if len(report.result.metadatas) > 1:
	    		%for metadata in report.result.metadatas:
	    		<%include file="/metadata/metadata.mako" args="metadata = metadata, report = report, primary = False"/>
	    		%endfor
    		%endif
    		##
    		## Silico methods.
    		<%include file="/method.mako" args="report = report"/>
    		##
    		<h4>Discussion</h4>
    		%for energy in [report.result.SCF_energies, report.result.MP_energies, report.result.CC_energies]:
	    		%if len(energy) > 0:
	    		<%include file="/energy/energy.mako" args="energy = energy, report = report"/>
	    		%endif
    		%endfor
    		%if "spin_density" in report.images:
    		<%include file="/spin.mako" args="report = report"/>
    		%endif
	    	%if len(report.result.alignment) > 0:
	    	<%include file="/geometry/geometry.mako" args="report = report"/>
	    	%endif
	    	%if report.result.dipole_moment is not None:
	    	<%include file="/dipole_moment/dipole_moment.mako" args="dipole_moment = report.result.dipole_moment, report = report, image_name = 'dipole_moment'"/>
	    	%endif
	    	%if report.result.transition_dipole_moment is not None:
	    	<%include file="/dipole_moment/dipole_moment.mako" args="dipole_moment = report.result.transition_dipole_moment, report = report, image_name = '{}_dipole'.format(report.result.transition_dipole_moment.excited_state.state_symbol)"/>
	    	%endif
	    	%if len(report.result.molecular_orbitals) > 0:
	    	<%include file="/orbital/orbitals.mako" args="report = report"/>
	    	%endif
	    	%if len(report.result.vibrations) > 0:
	    	<%include file="/vibrations/vibrations.mako" args="report = report"/>
	    	%endif
	    	%if len(report.result.excited_states) > 0:
	    	<%include file="/excited_states/excited_states.mako" args="report = report"/>
	    	%endif
	        <%include file="/emission/emission.mako" args="report = report"/>
	    </div>
	    <div class="section section--separate reportBody">
	    	<h4>Tables of Results</h4>
	    	%if len(report.result.alignment) > 0:
	    		<h5 class="h5--tableHeader">Atom Coordinates</h5>
	    		<%include file="/geometry/table.mako" args="report = report"/>
	    	%endif
	    	%if len(report.result.molecular_orbitals) > 0:
	    		<h5 class="h5--tableHeader">Molecular Orbitals</h5>
	    		<%include file="/orbital/table.mako" args="report = report"/>
	    	%endif
	    	%if len(report.result.vibrations) > 0:
	    		<h5 class="h5--tableHeader">Vibrational Frequencies</h5>
	    		<%include file="/vibrations/table.mako" args="report = report"/>
	    	%endif
    	</div>
    	%if len(report.result.excited_states) > 0 or len(report.result.spin_orbit_coupling) > 0:
	    	<div class="section section--separate">
		    	%if len(report.result.excited_states) > 0:
		    		<h5 class="h5--tableHeader">Excited States</h5>
		    		<%include file="/excited_states/table.mako" args="report = report"/>
		    	%endif
		    	%if len(report.result.spin_orbit_coupling) > 0:
		    		<h5 class="h5--tableHeader">Spin-Orbit Coupling</h5>
		    		<%include file="/SOC/table.mako" args="report = report"/>
		    	%endif
	    	</div>
    	%endif
    </body>
</html>