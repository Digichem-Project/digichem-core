## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## Define input for cc.
##
## We always enter CC because this is where we set memory usage and SCF density convergence, even if we're not using coupled-cluster...
## Enter cc.
cc
##
## Set memory (in MiB).
memory ${calculation.mem_per_cpu.MiB}
##
## Set SCF den conv.
%if calculation.scf['denconv'] is not None:
denconv 0.1E-${"{:02}".format(calculation.scf['denconv'])}
%endif
##
## Only set CC options if we're using CC.
%if calculation.ricc2['model'] is not None:
##
## Set aux. basis sets.
%for basis in ['cbas', 'cabs', 'jkbas']:
%if calculation.cc[basis] == "auto":
## Auto, enter nothing to use defaults.
${basis}
*
## Done.
%elif calculation.cc[basis] is not None:
## Explicit basis set, enter it.
${basis}
b all ${calculation.cc[basis]}
*
## Done.
%endif
## Done this basis set.
%endfor
##
## Done basis sets.
##
## ricc2 options next.
ricc2
##
## The method.
${calculation.ricc2['model']}
##
## SCS.
%if calculation.ricc2['scs'] is not None:
${calculation.ricc2['scs']}
%endif
##
## Optimisation options.
%if calculation.cc_geoopt['wavefunction'] is not None:
geoopt ${calculation.cc_geoopt['wavefunction']} ${calculation.cc_geoopt['state'] if calculation.cc_geoopt['state'] is not None else ""}
%endif
##
## End of ricc2.
*
##
## Excited state options.
%if calculation.ricc2_exci['nexc'] != 0:
##
## Enter exci menu.
exci
## 
## States to calculate.
irrep=${calculation.ricc2_exci['symmetry']} multiplicity=${calculation.ricc2_exci['multiplicity']} nexc=${calculation.ricc2_exci['nexc']}
##
## Whether to calculate oscillator strengths.
%if calculation.ricc2_exci['oscillator_strengths'] is not None:
spectrum states=all operators=${calculation.ricc2_exci['oscillator_strengths']}
%endif
##
## Whether to calculate excited state gradients.
%if calculation.ricc2_exci['gradients']:
xgrad states=(${calculation.ricc2_exci['symmetry']}{${calculation.ricc2_exci['multiplicity']}} 1-${calculation.ricc2_exci['nexc']})
%endif
##
## Done exci.
*
%endif
##
## Done real cc options.
%endif
##
## Exit cc menu.
*
## Done.