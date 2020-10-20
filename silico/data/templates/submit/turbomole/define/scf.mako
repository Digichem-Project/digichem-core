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
%if calculation.scf['iter'] is not None:
iter
${calculation.scf['iter']}
%endif
##
## SCF convergence.
%if calculation.scf['scfconv'] is not None:
conv
${calculation.scf['scfconv']}
%endif
##
## SCF Damping.
damp
##
## Set all 3 options (they all work the same way).
%for option in ("start", "step", "min"):
## Enter the value of the option.
${calculation.scf_damp[option] if calculation.scf_damp[option] is not None else ""}
## Done.
%endfor
##
## Exit scf.

## Done.