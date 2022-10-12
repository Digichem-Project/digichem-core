## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## Define input for dft.
##
%if calculation.method['dft']['calc']:
##
## Enter dft.
dft
##
## Turn on.
on
##
## Set functional.
func ${calculation.method['dft']['functional'].to_turbomole()}
##
## Change grid if we've been asked to.
%if calculation.method['dft']['grid'] is not None:
grid ${calculation.method['dft']['grid']}
%endif
##
## Exit dft.
*
##
## Done DFT
%endif
##
## Next, excited state properties for DFT.
%if calculation.properties['es']['calc'] and calculation.properties['es']['method'] in ["TD-HF", "CIS", "TD-DFT", "TDA"]:
ex
##
## Select our type of excited state based on TDA or no TDA and multiplicity.
## Turbomole doesn't recognise multiplicity if we're unrestricted.
%if calculation.is_unrestricted:
##
##
${"urpa" if calculation.properties['es']['method'] in ["TD-HF", "TD-DFT"] else "ucis"}
##
##
%else:
##
##
${"rpa" if calculation.properties['es']['method'] in ["TD-HF", "TD-DFT"] else "cis"}${"s" if calculation.properties['es']['multiplicity'] == "Singlet" else "t"}
##
##
%endif
##
## Exit ex.
*
##
## Select symmetry and number of states.
${calculation.properties['es']['symmetry']} ${calculation.properties['es']['num_states']}
*
##
## Select memory.
rpacor ${calculation.performance['memory'].MB}
*
##
## Skip question about SCF convergence (we set this elsewhere).
n
##
## Done excited states.
%endif
##
## Now dispersion correction, which weirdly is in a different menu.
%if calculation.method['dft']['dispersion'] is not None:
## Enter dsp menu.
dsp
##
## Set type of dsp (some have weird names).
%if calculation.method['dft']['dispersion'].upper() == "GD3":
on
%elif calculation.method['dft']['dispersion'].upper() == "GD2":
old
%elif calculation.method['dft']['dispersion'].upper() == "GD3BJ":
bj
%elif calculation.method['dft']['dispersion'].upper() == "GD4":
d4
%else:
${calculation.method['dft']['dispersion']}
%endif
##
## Turn on three-body params if asked.
%if calculation.method['dft']['dispersion_abc']:
abc
%endif
##
## Quit dsp.
*
%endif
##
## Done.