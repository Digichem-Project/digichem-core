## -*- coding: utf-8 -*-
##
## Template for displaying part of a calculation results in text format.
##
<%page args="title, result_name = ''"/>\
##
##
<%
	# Amount of whitespace to add either side of our title.
	padding = 5	
	
	longest_word = max(len(title), len(result_name))
	
	# Width of the decorative border (+2 is for the two pipe characters at either end).
	border_width = longest_word + 2*padding +2
%>\
## Header top.
${"-" * border_width}
## The title itself.
|${" " * padding}${title}${" " * (longest_word - len(title))}${" " * padding}|
##
## Now the result name (if given).
%if result_name != "":
|${" " * padding}${result_name}${" " * (longest_word - len(result_name))}${" " * padding}|
%endif
##
##
## Header bottom.
${"-" * border_width}