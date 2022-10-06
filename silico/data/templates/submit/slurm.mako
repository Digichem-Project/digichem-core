## -*- coding: utf-8 -*-
<%!
    from pathlib import Path
    import shlex
    import silico
%>\
##
<%page args="SLURM_target"/>\
##
#!/bin/bash
#SBATCH -J "${SLURM_target.program.calculation.descriptive_name} (${SLURM_target.program.meta['name']})"
##
%if SLURM_target.partition is not None:
#SBATCH -p ${SLURM_target.partition}
%endif
##
%if SLURM_target.num_tasks is not None:
#SBATCH --ntasks=${SLURM_target.num_tasks}
%endif
##
%if SLURM_target.CPUs_per_task is not None:
#SBATCH --cpus-per-task=${SLURM_target.CPUs_per_task}
%endif
##
%if SLURM_target.mem_per_cpu is not None:
#SBATCH --mem-per-cpu=${SLURM_target.mem_per_cpu}
%endif
##
%if SLURM_target.time is not None:
#SBATCH --time="${SLURM_target.time}"
%endif
##
%for option in SLURM_target.options:
#SBATCH --${option}="${SLURM_target.options[option]}"
%endfor
##
## Set slurm stderr and stdout
#SBATCH --output="${SLURM_target.calc_dir.log_file}"
#SBATCH --error="${SLURM_target.calc_dir.log_file}"
##
## exec into resume, so we handle signals from SLURM etc.
exec ${shlex.quote(SLURM_target.silico_command)} resume ${shlex.quote(str(SLURM_target.resume_file_path))}
