## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## Define input for dft.
##
%if calculation.dft['functional'] is not None:
##
## Enter dft.
dft
##
## Turn on.
on
##
## Set functional.
##func ${calculation.dft['functional']}
func ${calculation.func}
##
## Change grid if we've been asked to.
%if calculation.dft['grid'] is not None:
grid ${calculation.dft['grid']}
%endif
##
## Exit dft.
*
##
## Done DFT
%endif
##
## Next, excited state properties for DFT.
%if calculation.dft_exci['nexc'] > 0:
ex
##
## Select our type of excited state based on TDA or no TDA and multiplicity.
## Turbomole doesn't recognise multiplicity if we're unrestricted.
%if calculation.unrestricted_HF:
##
##
${"urpa" if not calculation.dft_exci['TDA'] else "ucis"}
##
##
%else:
##
##
${"rpa" if not calculation.dft_exci['TDA'] else "cis"}${"s" if calculation.dft_exci['multiplicity'] == 1 else "t"}
##
##
%endif
##
## Exit ex.
*
##
## Select symmetry and number of states.
${calculation.dft_exci['symmetry']} ${calculation.dft_exci['nexc']}
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
%if calculation.dispersion['dsp'] is not None:
## Enter dsp menu.
dsp
##
## Set type of dsp (some have weird names).
%if calculation.dispersion['dsp'].upper() == "GD3":
on
%elif calculation.dispersion['dsp'].upper() == "GD2":
old
%elif calculation.dispersion['dsp'].upper() == "GD3BJ":
bj
%elif calculation.dispersion['dsp'].upper() == "GD4":
d4
%else:
${calculation.dispersion['dsp']}
%endif
##
## Turn on three-body params if asked.
%if calculation.dispersion['abc']:
abc
%endif
##
## Quit dsp.
*
%endif
##
## Done.