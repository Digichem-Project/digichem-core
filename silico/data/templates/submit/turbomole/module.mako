## -*- coding: utf-8 -*-
<%!
    import shlex
%>\
##
<%page args="program, module"/>\
##
##
<%
    # Setup our module arguemnts.
    # For most, well behaved, arguments we can simply use the command that is returned from module.to_command().
    # However, some more complex modules require additional setup that we can only do here (because, for example,
    # they depend on setup that is performed after the calculation is configured.
    # Currently, this includes mp2prep and NumForce, both of which handle scratch differently to other modules.
    command = module.to_command()
    
    if program.scratch_base is not None and not module.no_scratch:
        if module.name == "mp2prep":
            command += " " + shlex.quote(str(program.scratch_base.resolve()))
        
        elif module.name == "NumForce":
            command += " -scrpath " + shlex.quote(str(program.scratch_base.resolve()))
        

%>\
##
##
#!/bin/bash
##
## Set the type of parallelization.
## This setting is used by the Turbomole init script to set our path to the correct binaries (SMP, MPI or normal).
%if program.calculation.performance['parallel_mode'] == "linear":
## Unset PARA_ARCH to use normal binaries.
unset PARA_ARCH
unset PARNODES
%else:
## Set PARA_ARCH appropriately.
export PARA_ARCH=${program.calculation.performance['parallel_mode']}
## Set number of nodes to run on.
export PARNODES=${program.calculation.performance['num_cpu']}
%endif
##
## IMPORTANT: The NumForce script is buggy and will try to forcably use MPI style parallisation
## if it detects it is inside of a SLURM environment. To prevent this, we unset SLURM_JOB_ID
## if we're using numforce and SMP.
%if program.calculation.performance['parallel_mode'] == "SMP" and module.name == "NumForce":
unset SLURM_JOB_ID
%endif
##
## Set scratch dir.
## IMPORTANT: A bug in NumForce will prevent it running properly if TURBOTMPDIR is set
## (probably because it handles scratch itself).
%if program.scratch_base is not None and module.name != "NumForce" and not module.no_scratch:
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
## If we've been asked to, set the allowed MKL instruction set.
%if program.calculation.performance['intel_AVX2_fix']:
##
export MKL_ENABLE_INSTRUCTIONS=SSE4_2
##
%endif
##
## Now run the calculation module we've been given.
${command} >> ${shlex.quote(str(program.turbomole_output_path.resolve()))} 2>&1 || { >&2 echo "Failed to execute Turbomole module '${module.name}'"; exit 1; }
##
## Done
exit $?
