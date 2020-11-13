## -*- coding: utf-8 -*-
##
<%page args="calculation" />\
##
<%!
	from silico.file.convert.gaussian import Gaussian_input_parser
%>\
##
<%
	# Get our input geometry in .com format.
	input_coords = Gaussian_input_parser(calculation.input_coords.to_format("com"))
%>\
## First, the chk file (if we're using one).
%if calculation.chk_file_name is not None:
%%Chk="${calculation.chk_file_name}"
%endif
##
## Next, the number of processes.
## There are two ways this can be specified; either with CPU_list or num_CPUs.
%if len(calculation.CPU_list) != 0:
%%CPU="${','.join(calculation.CPU_list)}"
%elif calculation.num_CPUs is not None:
%%NProcShared=${calculation.num_CPUs}
%endif
## Now the memory
%if calculation.memory is not None:
%%Mem=${calculation.memory}
%endif
##
## End of Link 0 commands.
##
## Route section.
#p ${calculation.route_section}

## Title card (has no effect on calculation but will crash on 'unusual' characters.
${calculation.safe_name(calculation.descriptive_name)}

## Charge and mult.
${calculation.charge}, ${calculation.multiplicity}
## Geometry
${input_coords.geometry}

## Extended Basis.
%for extended_basis_set in calculation.external_basis_sets:
${extended_basis_set.basis_set}
%endfor

##
## Extended ECP.
%for extended_ECP in calculation.external_ECPs:
${extended_ECP.ECP}
%endfor
##

## Finally, and additional sections that might be in the input file (currently disables).
##%for index, section in enumerate(calculation.input_file.sections[2:]):
##${section}
## Add a blank line (unless this is the last section and it contains a new line)
##%if index != (len(calculation.input_file.sections[2:]) -1) or (len(section) > 0 and section[-1] != "\n"):
##
##%endif
##%endfor