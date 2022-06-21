## -*- coding: utf-8 -*-

## This template is for displaying a citation in text. For a list of references as would be found in a bibliography, see references templates.

<%page args="name, numbers" />

<%
	if not isinstance(numbers, list):
		numbers = [numbers]
%>

<div class="citation"><div class="citation__name">${name}</div><div class="citation__number">[${",".join((str(number) for number in numbers))}]</div></div>