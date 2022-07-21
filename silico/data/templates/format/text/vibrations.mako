## -*- coding: utf-8 -*-
##
<%
    from silico.exception import Result_unavailable_error
%>\
##
<%page args="vibrations, result_name"/>\
##
<%
    if len(vibrations) == 0:
        raise Result_unavailable_error("vibrations", "there are no vibrational frequencies")
%>\
##
<%include file="title.mako" args="title='Vibrational Frequencies', result_name=result_name"/>
##
##
No. vibrations: ${len(vibrations)}
No. negative frequency: ${len(vibrations.negative)}
##
## We'll now show all negative vibrations.
%for index, vibration in enumerate(vibrations.negative):
Negative frequency ${index+1} /cm-1: ${"{:0.2f}".format(vibration.frequency)}
%endfor
