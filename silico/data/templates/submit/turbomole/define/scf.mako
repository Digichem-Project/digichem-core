## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## Define input for scf.
##
## Enter SCF.
scf
## 
## Set iterations.
%if calculation.scf['iterations'] is not None:
iter
${calculation.scf['iterations']}
%endif
##
## SCF convergence.
%if calculation.scf['convergence'] is not None:
conv
${calculation.scf['convergence']}
%endif
##
## SCF Damping.
damp
##
## Set all 3 options (they all work the same way).
%for option in ("weight", "step", "min"):
## Enter the value of the option.
${calculation.scf['damping'][option] if calculation.scf['damping'][option] is not None else ""}
## Done.
%endfor
##
## Exit scf.

## Done.