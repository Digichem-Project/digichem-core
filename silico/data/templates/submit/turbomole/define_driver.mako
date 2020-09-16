## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## This file is used as input to Turbomole's 'define' program, which sets up a calculation.
## This is very hacky; define does not have a well defined API (?).
##
## The first question allows us to specify a default control file, we skip this for now but it might prove useful in the future.

## Next is the job title.
${calculation.safe_name(calculation.descriptive_name)}
##
## Load atom coordinates.
a ../Input/coord
##
## Load redundant coords.
%if calculation.redundant_internal_coordinates:
ired
%endif
##
## Next section.
*
## 
## Load basis set.
b all ${calculation.basis_set}
##
## Next section.
*
##
## Setup occupancy.
eht
##
## Accept defaults.

##
## Molecule charge.
${calculation.charge}
##
## Accept.

##
## Now our method(s).
%for method in methods:
## Enter method.
${method}
##
## Set memory
memory ${calculation.memory}
##
## Next keywords.
%for keyword in keywords
##
## Enter the keyword section.
${keyword}
##
## Now print each option in the keyword.
${"\n".join(calculation.keywords[keyword])}
##
## Exit keyword section.
*
%endfor
##
## Exit method section.
*
%endfor
##
## All done.
*
##