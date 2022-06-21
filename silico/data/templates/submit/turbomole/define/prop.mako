## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
##
<%
	# There is a turbomole bug that means we cannot list more than a single line of orbitals to print, which maxes out around 78 characters.
	# To overcome this, we ask for all orbitals between the min and max.
	# This may print more orbitals than we need, and so is somewhat wasteful.
	if len(calculation.analysis['plt']['orbitals']) > 0:
		orbitals = "{}-{}".format(min(calculation.analysis['plt']['orbitals']), max(calculation.analysis['plt']['orbitals']))
	else:
		orbitals = None
%>\
##
## Define input for prop (orbital analysis etc).
##
## Only continue if we've got some orbitals to plot.
##%if len(calculation.analysis['plt']['orbitals']) > 0:
%if calculation.analysis['plt']['calculate']:
## First, enter prop
prop
##
## Next, plt
plt
##
## Answer yes to modify.
y
##
##
## Orbital list.
%if orbitals is not None:
m ${orbitals} ${"dens" if calculation.analysis['plt']['density'] else ""}
%endif
##
## Also set format.
f ${calculation.analysis['plt']['format']}
##
##
## Exit plt.
*
##
## Exit prop.
*
##
##
%endif
##
## Done.