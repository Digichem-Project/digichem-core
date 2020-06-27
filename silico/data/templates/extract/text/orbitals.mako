## -*- coding: utf-8 -*-
##
<%
	from silico.exception import Result_unavailable_error
%>\
##
<%page args="orbitals, result_name"/>\
##
##
<%
	if len(orbitals) == 0:
		raise Result_unavailable_error("orbitals", "there are no orbitals of the requested type") 

	title = "Orbitals"
	if orbitals.spin_type != "none":
		title += " ({})".format(orbitals.spin_type)
%>\
##
<%include file="title.mako" args="title=title, result_name=result_name"/>
##
HOMO eV: ${"{:0.2f}".format(orbitals.HOMO_energy)}
LUMO eV: ${"{:0.2f}".format(orbitals.LUMO_energy)}
HOMO/LUMO energy /eV: ${"{:0.2f}".format(orbitals.HOMO_LUMO_energy)}
No. virtual: ${len([orbital for orbital in orbitals if orbital.HOMO_difference > 0])}
No. occupied: ${len([orbital for orbital in orbitals if orbital.HOMO_difference <= 0])}
