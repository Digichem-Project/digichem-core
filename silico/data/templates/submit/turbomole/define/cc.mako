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
%if calculation.scf['density_convergence'] is not None:
denconv 0.1E-${"{:02}".format(calculation.scf['density_convergence'])}
%endif
##
## Set aux. basis set, if we're using RI-C.
%if calculation.method['ri']['correlated']['calc']:
%if str(calculation.method['ri']['correlated']['basis_set']) == "auto":
## Auto, enter nothing to use defaults.
cbas
*
## Done.
%else:
## Explicit basis set, enter it.
cbas
b all ${calculation.method['ri']['correlated']['basis_set'].to_turbomole()}
*
## Done.
%endif
##
## Done this basis set.
%endif
##
## Done basis sets.
##
## Even though this menu is called 'ricc2', we actually use it for most higher order methods,
## including thise run with ccsdf12 (MP3, MP4, CCSD etc.)
%if calculation.post_HF_method:
##
## ricc2 options next.
ricc2
##
## The method. This needs to be lower case.
${calculation.post_HF_method.lower()}
##
## SCS.
%if calculation.method['scs']['calc']:
${calculation.scs_line}
%endif
##
## Optimisation options.
## Geoopt actually selects whether to calculate gradients or not.
## We need these for both opt and freq jobs.
## TODO: Need to go over how geopt really works, the input for excited states and ground are quite different.
## TOOD: HERE: Need to figure out geoopt syntax for excited state triplets...
%if calculation.properties['opt']['calc'] and calculation.properties['es']['calc']:
geoopt ${calculation.post_HF_method} (${calculation.properties['es']['symmetry']}{${3 if calculation.properties['es']['multiplicity'] == "Triplet" else 1}} ${calculation.properties['es']['state_of_interest']})
%elif calculation.properties['opt']['calc'] or calculation.properties['freq']['calc']:
geoopt ${calculation.post_HF_method} (${calculation.properties['opt']['ricc2']['optimise_symmetry']} ${calculation.optimise_multiplicity})
%endif
##
## End of ricc2.
*
##
##
%endif
##
## Excited state options.
%if calculation.properties['es']['calc'] and calculation.properties['es']['method'] in ["CIS(D)", "CIS(Dinf)", "ADC(2)", "CC2", "CCS", "CCSD"]:
##
## Enter exci menu.
exci
## 
## States to calculate.
${calculation.exci_irrep_line}
##
## Whether to calculate oscillator strengths.
%if calculation.properties['es']['ricc2']['spectrum_operators']:
spectrum states=all operators=${calculation.properties['es']['ricc2']['spectrum_operators']}
%endif
##
## Whether to calculate excited state sp_gradients.
%if calculation.properties['es']['ricc2']['sp_gradients']:
${calculation.exci_gradient_line}
%endif
##
## Done exci.
*
%endif
##
## Exit cc menu.
*
## Done.