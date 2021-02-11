## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
##
<%
	# There is a turbomole bug that means we cannot list more than a single line of orbitals to print, which maxes out around 78 characters.
	# To overcome this, we ask for all orbitals between the min and max.
	# This may print more orbitals than we need, and so is somewhat wasteful.
	if len(calculation.plt['orbitals']) > 0:
		orbitals = "{}- {}".format(min(calculation.plt['orbitals']), max(calculation.plt['orbitals']))
	else:
		orbitals = None
%>\
##
## Define input for prop (orbital analysis etc).
##
## Only continue if we've got some orbitals to plot.
##%if len(calculation.plt['orbitals']) > 0:
%if orbitals is not None:
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
##m ${" ".join(calculation.plt['orbitals'])}
m ${orbitals}
##
## Also set format.
f ${calculation.plt['format']}
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