# General imports.
import os

# Silico imports.
from silico.interface.urwid.file.browser import File_selector
from silico.config.configurable.option import Option
from silico.result.alignment.base import Alignment, Minimal


class Result_selector(File_selector):
    """
    A file selector for loading calculation result files.
    """

    num_cpu = Option(help = "The number of processes to use in parallel to parse chosen results.", type = int, default = os.cpu_count())
    alignment = Option(help = "The alignment method to use to reorientate molecules.", type = str, choices = Alignment.known_handles(), default = Minimal.CLASS_HANDLE[0])
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, can_choose_folders = True, can_choose_multiple = True, **kwargs)
