## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
## Define input for RI (resolution of the identity approximation).
##
## Both the RI-J and RI-HF submenus have a memory and 'on' switch, but both effect the same options.
## For the 'on'switch, this might be a bug because RI-J and RI-HF are activated using different
## options in the define fine ($rij and $rijk respectively), but only $rij is set using the menu.
## Fortunately, these options don't seem to do much except force dscf to crash (correctly forcing use
## of ridft instead).
##
## 
## First, handle RI-J ($rij and $jbas).
%for menu_name, basis_name, ri_type in [("ri", "jbas", "coulomb"), ("rijk", "jkbas", "hartree_fock")]:
%if calculation.method['ri'][ri_type]['calc']:
##
## Enter the menu.
${menu_name}
##
## 'Turn on' (probably just sets a flag for Turbomole to check for and abort dscf if neccessary).
on
##
## Assign basis set.
${basis_name}
##
## If we're using an auto basis set, leave as default.
%if str(calculation.method['ri'][ri_type]['basis_set']) == "auto":
## Auto, enter nothing to use defaults.
*
## Done.
%else:
## Explicit basis set, enter it.
b all ${calculation.method['ri'][ri_type]['basis_set'].to_turbomole()}
*
## Done.
%endif
##
## Done this RI type.
## Exit this submenu
*
%endif
##
##
%endfor
##
## Done.