## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## Define input for prop (orbital analysis etc).
##
## Only continue if we've got some orbitals to plot.
%if len(calculation.plt['orbitals']) > 0:
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
m ${" ".join(calculation.plt['orbitals'])}
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