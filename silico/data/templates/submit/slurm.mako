## -*- coding: utf-8 -*-
<%!
	from pathlib import Path
%>\
##
#!/bin/bash
#SBATCH -J "${SLURM_target.program.calculation.descriptive_name}"
##
%if SLURM_target.partition is not None:
#SBATCH -p ${SLURM_target.partition}
%endif
##
%if SLURM_target.num_tasks is not None:
#SBATCH -N ${SLURM_target.num_tasks}
%endif
##
%if SLURM_target.CPUs_per_task is not None:
#SBATCH -n ${SLURM_target.CPUs_per_task}
%endif
##
%if SLURM_target.mem_per_CPU is not None:
#SBATCH --mem_per_CPU=${SLURM_target.mem_per_CPU}
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
#SBATCH --output="${Path(SLURM_target.calc_dir.output_directory, 'slurm.out')}"
#SBATCH --error="${Path(SLURM_target.calc_dir.output_directory, 'slurm.out')}"
##
##_silico_resume "${SLURM_target.resume_file_path.relative_to(SLURM_target.calc_dir.input_directory)}"
## exec into resume, so we handle signals from SLURM etc.
exec _silico_resume "${SLURM_target.resume_file_path}"
<%page args="SLURM_target"/>