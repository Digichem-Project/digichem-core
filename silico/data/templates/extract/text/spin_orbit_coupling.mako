## -*- coding: utf-8 -*-
##
<%
	from silico.exception import Result_unavailable_error
%>\
##
<%page args="SOC_list, result_name = ''"/>\
##
##
<%
	if len(SOC_list) == 0:
		raise Result_unavailable_error("SOC", "there is no spin-orbit coupling")
		
	soc_list = []
	for state1, state2 in [("S(0)", "T(1)"), ("S(1)", "T(1)")]:
		# Get parameters.
		try:
			soc_list.append(SOC_list.between(state1, state2))
		except Exception:
			# Skip.
			raise
			continue

%>\
##
<%include file="title.mako" args="title='Spin-Orbit Coupling', result_name=result_name"/>
##
%for soc in soc_list:
<${soc.singlet_state.state_symbol}|Hso|${soc.triplet_state.state_symbol}> /cm-1: ${"{:0.2f}".format(soc.wavenumbers)}
<${soc.singlet_state.state_symbol}|λ|${soc.triplet_state.state_symbol}>: ${"{:0.2f}".format(soc.mixing_coefficient)}
%endfor
