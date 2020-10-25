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
## What we do next depends on whether we are setting up a RHF or UHF calc.
%if calculation.multiplicity == 0:
##
## RHF, accept defaults.

%else:
##
## Reject defaults.
n
##
## Set the multiplicity (for some reason turbomole actually asks for the num of unpaired e-, which is m-1, but whatever)
u ${calculation.unpaired_electrons}
##
## Next.
*
##
## Confirm.

##
## Done EHT.
%endif
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
## ! SCF                     !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
<%include file="scf.mako" args="calculation = calculation"/>\
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
## ! DFT                     !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
<%include file="dft.mako" args="calculation = calculation"/>\
##
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
## ! CC                      !
## !!!!!!!!!!!!!!!!!!!!!!!!!!!
<%include file="cc.mako" args="calculation = calculation"/>\
##
## All done.
*
##