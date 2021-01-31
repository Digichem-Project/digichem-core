## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## Define input for prop (orbital analysis etc).
##
## First, enter prop
prop
##
## Next, plt
plt
##
## Answer yes to modify.
y
##
## List our orbitals if we have some.
%if len(calculation.plt.orbitals) > 0:
##
## Orbital list.
m ${" ".join(calculation.plt.orbitals)}
##
## Also set format.
f ${calculation.plt.format}
##
## Done.
%endif
##
## Exit plt.
*
##
## Exit prop.
*
## Done.