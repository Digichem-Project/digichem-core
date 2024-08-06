# Alignment procedures based on the force technique
import math
import numpy
from statistics import  mean

from digichem.result.alignment import Alignment

class Pivot(Alignment):
    """
    The pivot method for determining molecule alignment. 
    """

    # Names that uniquely describe this alignment protocol.
    CLASS_HANDLE = ["PIVOT"]

    def align_axes(self):
        """
        Realign the axes of our coordinate system.
        
        :return: Nothing. The atoms are rearranged in place.
        """
        # Begin by translating to centre of coordinates.
        coords = self.get_coordinate_list()
        self.translate((-mean(coords[0]) , -mean(coords[1]), -mean(coords[2])))

        # XY
        print(self.X_length * self.Y_length)
        average = sum(
            (self.get_absolute_theta(atom.coords[1], atom.coords[0]) for atom in self)
        )
        self.rotate_XY(average)
        print(self.X_length * self.Y_length)

        # XZ
        print(self.Z_length * self.Y_length)
        average = sum(
            (self.get_absolute_theta(atom.coords[0], atom.coords[2]) for atom in self)
        )
        self.rotate_XZ(-average)
        print(self.Z_length * self.Y_length)
        print(average)

    
    @classmethod
    def get_absolute_theta(self, opposite, adjacent):
        """
        Determine the angle needed (theta) to rotate a point, always favouring the positive side of the axis.
        
        This function does atan(opposite / adjacent), watching out for div by 0.
        :param opposite: the coordinate of the point along the axis opposite to the angle to calculate. After rotation, the point will have this coord == 0.
        :param adjacent: The coordinate of the point along the axis adjacent to the angle to calculate. After rotation, the point will be on this axis.
        :return: The angle (in radians).
        """
        try:
            angle =  math.atan(opposite / adjacent)

        except (FloatingPointError, ZeroDivisionError):
            # Think it's safe to assume the angle is 90 degrees in this instance.
            angle =  math.pi /2
        
        # TODO: There must be a smarter way of doing this?
        if adjacent > 0:
            return angle

        elif opposite > 0:
            return angle + math.pi

        else:
            return angle - math.pi