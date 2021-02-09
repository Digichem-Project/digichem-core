## -*- coding: utf-8 -*-
<%!
    import shlex
%>\
##
<%page args="program"/>\
##
##
#!/bin/bash
##
## Set the type of parallelization.
## This setting is used by the Turbomole init script to set our path to the correct binaries (SMP, MPI or normal).
%if program.calculation.parallel_mode == "linear":
## Unset PARA_ARCH to use normal binaries.
unset PARA_ARCH
unset PARNODES
%else:
## Set PARA_ARCH appropriately.
export PARA_ARCH=${program.calculation.parallel_mode}
## Set number of nodes to run on.
export PARNODES=${program.calculation.num_CPUs}
%endif
##
## Set scratch dir.
%if program.scratch_base is not None:
##
## Yes scratch.
export TURBOTMPDIR=${shlex.quote(str(program.scratch_base.resolve()))}
%else:
##
## No scratch.
unset TURBOTMPDIR
%endif
##
## Import turbomole init script.
<%include file="wrapper.mako"/>\
##
## Now run the calculation programs we've been given.
%for module in program.calculation.modules:
${module} >> ${shlex.quote(str(program.turbomole_output_path.resolve()))} || { >&2 echo "Failed to execute Turbomole program '${module}'"; exit 1; }
##${module} >> job.last || { >&2 echo "Failed to execute Turbomole program '${module}'"; exit 1; }
%endfor
##
## Done
exit $?