## -*- coding: utf-8 -*-
##
<%!
	from silico.exception import Result_unavailable_error
%>\
##
<%page args="alignment, result_name = ''"/>\
##
<%
	if len(alignment) == 0:
		raise Result_unavailable_error("geometry", "there are no atoms")
%>\
##
<%include file="title.mako" args="title='Geometry', result_name=result_name"/>
##
Formula: ${alignment.formula_string}
Exact mass /gmol-1: ${"{:0.4f}".format(alignment.safe_get('mass')) if alignment.safe_get('mass') is not None else None}
Molar mass /gmol-1: ${"{:0.4f}".format(alignment.molar_mass)}
No. atoms: ${len(alignment)}
Alignment method: ${alignment.CLASS_HANDLE[0]}
X extension /Å: ${"{:0.2f}".format(alignment.X_length)}
Y extension /Å: ${"{:0.2f}".format(alignment.Y_length)}
Z extension /Å: ${"{:0.2f}".format(alignment.Z_length)}
Linearity ratio: ${"{:0.2f}".format(alignment.get_linear_ratio())}
Planarity ratio: ${"{:0.2f}".format(alignment.get_planar_ratio())}
