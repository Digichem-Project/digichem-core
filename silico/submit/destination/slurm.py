# General imports.
from mako.lookup import TemplateLookup
from pathlib import Path
import os
import stat
import subprocess
from subprocess import CalledProcessError

# Silico imports.
import silico
from silico.exception import Submission_error, Silico_exception
from silico.submit import Memory
from silico.submit.destination.resume import Resumable_destination
from silico.config.configurable.option import Option
from silico.misc.base import is_int


class SLURM(Resumable_destination):
    """
    Implementation to allow submission to SLURM, a popular scheduling system.
    """
    
    CLASS_HANDLE = ("SLURM",)
    
    # The name of the script which we pass to sbatch.
    SBATCH_SCRIPT_NAME = "sbatch.submit"
    
    # Configurable options.
    partition = Option(help = "SLURM partition name", type = str, required = True)
    time = Option(help = "Max allowed job time (dd-hh:mm:ss)", default = None, type = str)
    num_tasks = Option(help = "Number of tasks to allocate", default = 1, type = int)
    _CPUs_per_task = Option("CPUs_per_task", help = "Number of CPUs (per task) to allocate. Use 'auto' to use the same number as required for the calculation", default = 'auto', validate = lambda option, configurable, value: value == "auto" or is_int(value))
    _mem_per_CPU = Option("mem_per_CPU", help = "The amount of memory (per CPU) to allocate. Use 'auto' to automatically decide how much memory to allocate based on that required for the calculation. Leave blank to use server defaults", default = "auto", validate = lambda option, configurable, value: value in ["auto", None] or is_int(value))
    options = Option(help = "Additional SLURM options. Any option valid to SLURM can be included here", default = {}, type = dict)
    sbatch_command = Option(help = "The name/path of the sbatch command", default = "sbatch", type = str)
    sinfo_command = Option(help = "The name/path of the sinfo command", default = "sinfo", type = str)
    silico_command = Option(help = "The name/path of the silico command", default = "silico", type = str)
    
    def get_num_nodes(self, idle = False):
        """
        Get the current number of nodes for this partition.
        
        :raises Silico_exception: If the number of nodes could not be determined.
        :param idle: Whether to get the number of idle nodes or total nodes.
        :return: The number of nodes.
        """
        # The signature we'll use to call sinfo.
        sig = [
            self.sinfo_command,
            "-p", self.partition,
            "-o", "%D",
            "-h"
        ]
        
        if idle:
            sig.extend(("-t", "IDLE",))

        try:
            # Call sinfo.
            done = subprocess.run(
                sig,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True,
                check = True,
            )
            
            # The number of nodes will be our output (unless there were none, in which case nothing is outputed).
            num_nodes = int(done.stdout) if done.stdout != "" else 0
        except Exception:
            raise Silico_exception("Failed to retrieve sinfo for partition {}".format(self.partition))
        
        return num_nodes
    
    def get_CPU_info(self):
        """
        Get information on the current allocation of CPUs in this partition.
        
        :return: A tuple of integers specifying the current CPU allocation in the format (allocated, idle, other, total).
        """
        # The signature we'll use to call sinfo.
        sig = [
            self.sinfo_command,
            "-p", self.partition,
            "-o", "%C",
            "-h"
        ]

        try:
            # Call sinfo.
            done = subprocess.run(
                sig,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True,
                check = True,
            )
            
            # Our output is a string of the form alloc/idle/other/total (unless the partition isn't valid, in which case we get nothing).
            cpu_info = done.stdout.split("/")
            
            # This will trigger an error if the format is weird, which is what we want.
            cpu_info = (int(cpu_info[0]),int(cpu_info[1]),int(cpu_info[2]),int(cpu_info[3]))
        except Exception as e:
            raise Silico_exception("Failed to retrieve sinfo for partition {}".format(self.partition)) from e
        
        return cpu_info
    
    @property
    def status(self):
        """
        This method is called to get 'status' about this destination.
        
        For SLURM, we use sinfo to get the number of free nodes for this partition.
        
        :raises Silico_exception: If the number of nodes could not be determined.
        :return: Info string.
        """
        cpu_info = self.get_CPU_info()
        num_nodes = self.get_num_nodes()
        free_nodes = self.get_num_nodes(True)
        return "{} ({:0.0f}%) idle nodes, {} ({:0.0f}%) idle CPUs".format(free_nodes, (free_nodes/num_nodes)*100 if num_nodes != 0 else 0, cpu_info[1], (cpu_info[1]/cpu_info[3])*100 if cpu_info[3] != 0 else 0)
    
    
    ############################
    # Class creation mechanism #
    ############################
    
    class _actual(Resumable_destination._actual):
        """
        Inner class for SLURM.
        """
        
        @property
        def unique_name(self):
            """
            Get a name that is unique for this calculation instance.
            
            SLURM tries to get a unique name based on our allocated SLURM ID, but if we have not yet been submitted to SLURM this method will fallback to a different method.
            """
            if getattr(self, "_unique_name", None) is None:
                try:
                    self._unique_name = os.environ['SLURM_JOB_ID']
                except KeyError:
                    # SLURM_JOB_ID isn't set, means we haven't been submitted to SLURM yet.
                    # Fallback.
                    self._unique_name = super().unique_name
            
            return self._unique_name
    
        @property
        def CPUs_per_task(self):
            """
            Get the number of CPUs to assign for this calculation.
            
            This property will resolve 'auto' to an actual number of processors, use _CPUs_per_task if you do not want this behaviour.
            """
            if self._CPUs_per_task == "auto":
                return self.program.calculation.num_CPUs if self.program.calculation.num_CPUs is not None else 1
            else:
                return self._CPUs_per_task
        
        @CPUs_per_task.setter
        def CPUs_per_task(self, value):
            """
            Set the number of CPUs to assign for this calculation.
            """
            self._CPUs_per_task = value
        
        @property
        def mem_per_CPU(self):
            """
            Get the amount of memory to assign (per CPU).
            
            This property will resolve 'auto' to an actual amount of memory, use _mem_per_CPU if you do not want this behaviour.
            """
            if self._mem_per_CPU is None:
                return None
            elif self._mem_per_CPU == "auto":
                return SLURM_memory(round(float(float(self.program.calculation.memory) * 1.1)) / self.CPUs_per_task)
            else:
                return SLURM_memory(float(self._mem_per_CPU))
    
        @mem_per_CPU.setter
        def mem_per_CPU(self, value):
            """
            Set the amount of memory to assign (per CPU).
            """
            self._mem_per_CPU = value
        
        def write_sbatch_script(self):
            """
            Write, to file, the control script which is passed to sbatch.
            """
            # Get and load our template.
            template_body = TemplateLookup(directories = str(silico.default_template_directory())).get_template("/submit/slurm.mako").render_unicode(SLURM_target = self)
            
            with open(self.sbatch_script_path, "wt") as sbatch_file:
                sbatch_file.write(template_body)
                
            # Make it executable.
            os.chmod(self.sbatch_script_path, os.stat(self.sbatch_script_path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    
        @property
        def sbatch_script_path(self):
            """
            Path to the bash script file that we will write and pass to sbatch.
            """
            return Path(self.calc_dir.input_directory, self.SBATCH_SCRIPT_NAME)
    
        def pre(self):
            """
            Submission method called before pausing.
            """
            # Write the control file.
            self.write_sbatch_script()
            
        def post(self):
            """
            Submission method called after pausing.
            """
            # Call sbatch.
            try:
                subprocess.run(
                    [self.sbatch_command, self.sbatch_script_path],
                    stdout = subprocess.PIPE,
                    stderr = subprocess.STDOUT,
                    check = True,
                    universal_newlines = True,
                )
            except CalledProcessError as e:
                # Something went wrong.
                e.__context__ = None
                
                raise Submission_error(self, "{} did not exit successfully:\n{}".format(self.sbatch_command, e.stdout)) from e
            except FileNotFoundError as e:
                # Couldn't find sbatch.
                e.__context__ = None
                raise Submission_error(self, "unable to locate sbatch executable '{}'".format(self.sbatch_command)) from e
            except Exception as e:
                raise e from None
        
    
class SLURM_memory(Memory):
    
    # SLURM has its own set of units...
    UNITS = {
        'T': 1000000000000,
        'G': 1000000000,
        'M': 1000000,
        'K': 1000,
        }
            
            