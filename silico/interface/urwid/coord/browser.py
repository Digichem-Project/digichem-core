from silico.interface.urwid.file.browser import File_selector
from silico.config.configurable.option import Option


class Coord_selector(File_selector):
    """
    A file selector for loading coordinate files.
    """

    charge = Option(help = "Forcibly set the molecular charge (an integer) of newly loaded files. If blank, the charge will be determined from each loaded file.", type = int, default = None)
    multiplicity = Option(help = "Forcibly set the molecular multiplicity (an integer) of newly loaded files. If blank, the multiplicity will be determined from each loaded file.", type = int, default = None)
    generate_3D = Option(help = "Whether to convert 2D coordinates to 3D by performing a rapid optimisation with molecular mechanics. Even if True, 3D coordinates will only be generated if it can be safely determined that the coordinates are not already in 3D.", type = bool, default = True)
