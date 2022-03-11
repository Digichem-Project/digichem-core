## -*- coding: utf-8 -*-
##
<%page args="calculation, write_geom = True" />\
##
<%!
    from silico.file.convert.gaussian import Gaussian_input_parser
%>\
##
## First, the chk file.
%%Chk="${calculation.chk_file_name}"
##
## Next, the rwf file.
%%Rwf="${calculation.rwf_file_name}"
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

## Geometry (if we've been asked to).
%if write_geom:
## Charge and mult.
${calculation.charge}, ${calculation.multiplicity}
## Geometry
<%
    # Get our input geometry in .com format.
    input_coords = Gaussian_input_parser(calculation.input_coords.to_format("com"))
%>\
##
${input_coords.geometry}
##
## This new line is important.
##

##
## External basis set from BSE
##
%if calculation.basis_set['exchange'] is not None:
${calculation.basis_set['exchange'].to_format("gaussian94", calculation.input_coords.elements)}
%endif
##
%endif