## -*- coding: utf-8 -*-

<%!
    from silico.exception import Result_unavailable_error
%>

<%page args="result" />

<%
    # Get our transition_dipole_moment if we can.
    try:
        S1 = result.excited_states.get_state("S(1)")
        transition_dipole_moment = S1.transition_dipole_moment
    except Result_unavailable_error:
        # No S1 available.
        transition_dipole_moment = None
%>

<div class="section section--summary section--fullPage">
    <h2 class="section__header">Summary of Results</h2>
    <div class="section__body section__body--summary">
    %if result.metadata is not None:
        <%include file="/metadata/metadata_results.mako" args="metadata = result.metadata"/>
    %endif
    %if len(result.metadatas) > 1:
    	%for index,sub_metadata in enumerate(result.metadatas):
    		<%include file="/metadata/metadata_results.mako" args="metadata = sub_metadata, title = 'Calculation {}'.format(index +1)"/>
    	%endfor
    %endif
    %if len(result.energies.scf) > 0:
        <%include file="/energy/energy_results.mako" args="energies = result.energies.scf"/>
    %endif
    %if len(result.energies.mp) > 0:
        <%include file="/energy/energy_results.mako" args="energies = result.energies.mp"/>
    %endif
    %if len(result.energies.cc) > 0:
        <%include file="/energy/energy_results.mako" args="energies = result.energies.cc"/>
    %endif
    %if len(result.alignment) > 0:
        <%include file="/geometry/geometry_results.mako" args="alignment = result.alignment"/>
    %endif
    %if len(result.orbitals) > 0:
        <%include file="/orbitals/HOMO_LUMO_results.mako" args="orbitals = result.orbitals"/>
    %endif
    %if len(result.beta_orbitals) > 0:
        <%include file="/orbitals/HOMO_LUMO_results.mako" args="orbitals = result.beta_orbitals"/>
    %endif
    %if result.pdm is not None:
        <%include file="/dipole_moment/dipole_moment_results.mako" args="dipole_moment = result.pdm"/>
    %endif
    %if transition_dipole_moment is not None:
        <%include file="/dipole_moment/dipole_moment_results.mako" args="dipole_moment = transition_dipole_moment"/>
    %endif
    %if len(result.vibrations) > 0:
        <%include file="/vibrations/vibrations_results.mako" args="vibrations = result.vibrations"/>
    %endif
    %for multiplicity, vertical in result.emission.vertical.items():
        <%include file="/emission/emission_results.mako" args="emission = vertical"/><span></span>
    %endfor
    %for multiplicity, adiabatic in result.emission.adiabatic.items():
        <%include file="/emission/emission_results.mako" args="emission = adiabatic"/><span></span>
    %endfor
    %if len(result.excited_states) > 0 and len(result.soc) > 0:
        <%include file="/soc/SOC_results.mako" args="soc = result.soc, excited_states = result.excited_states"/>
    %endif
    %if len(result.excited_states) > 0:
        <%include file="/excited_states/excited_states_results.mako" args="excited_states = result.excited_states"/><span></span>
    %endif
    </div>
</div>