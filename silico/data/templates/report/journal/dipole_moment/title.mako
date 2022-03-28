## -*- coding: utf-8 -*-
##
## Two functions for formatting the name of a PDM or TDM section.
## They are the same, except of capitalisation.
##
## Without capitals.
##
<%def name="dipole_name(dipole_moment)">
	%if dipole_moment.dipole_type == "permanent":
permanent dipole moment
    %else:
transition (${dipole_moment.excited_state.multiplicity_symbol}<sub>${dipole_moment.excited_state.multiplicity_level}</sub>) dipole moment
    %endif
</%def>
##
##
## With capitals.
##
<%def name="dipole_title(dipole_moment)">
	%if dipole_moment.dipole_type == "permanent":
Permanent dipole moment
    %else:
Transition (${dipole_moment.excited_state.multiplicity_symbol}<sub>${dipole_moment.excited_state.multiplicity_level}</sub>) dipole moment
    %endif
</%def>