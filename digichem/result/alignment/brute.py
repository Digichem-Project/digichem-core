# Methodologies for determining molecule linearity.
import math
import numpy

from digichem.result.alignment import Alignment

class Brute_force(Alignment):
    """
    The brute force method for determining molecule linearity. 
    """

    def __init__(self, atoms, *args, charge = None, steps = 180, **kwargs):
        # Number of steps to take in each axis.
        # The total number of iterations will be steps^3
        self.steps = steps
        super().__init__(atoms, *args, charge=charge, **kwargs)
        
    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["Brute force", "Brute"]
        
    def align_axes(self):
        """
        Realign the axes of our coordinate system, so they have the following meaning:
        
        X-axis: The long axis (the kebab skewer), which we define as passing through the pair of atoms with the greatest separation in our set.
        Y-axis: The middle axis, which we define as perpendicular to the the X-axis (obviously) and passing through the atom that is furthest from the X-axis. Note that the Y-axis only has to pass through one atom; there may not be a corresponding atom on the other side (but there will be if the molecule is symmetrical about the X-axis).
        Z-axis: The short axis, defined as perpendicular to both the X and Y-axes (so we have no choice where this goes).
        
        :return: Nothing. The atoms are rearranged in place.
        """
        # Start by creating a copy of our coordinates to experiment with.
        points = [(atom.coords[0], atom.coords[1], atom.coords[2]) for atom in self]

        smallest = math.inf

        for x_angle in numpy.linspace(0, 0.5*math.pi, self.steps):
            for y_angle in numpy.linspace(0, 0.5*math.pi, self.steps):
                for z_angle in numpy.linspace(0, 0.5*math.pi, self.steps):
                    new_points = [
                        self.rotate_coords_XY(
                            self.rotate_coords_XZ(
                                self.rotate_coords_YZ(coord, x_angle), y_angle), z_angle)
                        for coord in points
                    ]
                    #x, y, z = self.get_coordinate_list(new_points)
                    x, y, z = list(map(list, zip(* [(coord[0], coord[1], coord[2]) for coord in new_points])))

                    X = max(x) - min(x)
                    Y = max(y) - min(y)
                    Z = max(z) - min(z)

                    volume = X * Y * Z
                    # print("{}, {}, {}, {}".format(x_angle, y_angle, z_angle, volume))
                    if volume < smallest:
                        smallest = volume

        print("Smallest: {}".format(smallest))

    
    