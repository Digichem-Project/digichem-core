## -*- coding: utf-8 -*-
##
<%!
    from silico.exception import Result_unavailable_error
%>\
##
<%page args="energy, result_name = ''"/>\
##
<%
    if len(energy) == 0:
        raise Result_unavailable_error(energy.energy_type, "there is no energy of the requested type")
%>\
##
<%include file="title.mako" args="title=energy.energy_type + ' Energy', result_name=result_name"/>
##
##
No. of steps: ${len(energy)}
Energy /eV: ${"{:0.4f}".format(energy.final)}
Energy /kJmol-1: ${"{:0.4f}".format(energy.eV_to_kJmol(energy.final))}
