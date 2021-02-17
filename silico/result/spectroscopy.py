import math
import numpy
import itertools
import scipy.constants
import scipy.signal

import silico.result.excited_states
from silico.exception.base import Silico_exception

class Spectroscopy_graph():
    """
    Top level class for graphing spectroscopy results (energy vs intensity).
    
    For generating pictures of these graphs, see silico.image.spectroscopy
    """
    
    def __init__(self, coordinates):
        """
        Constructor for Spectroscopy_graph objects
        
        :param coordinates: A list of (energy, intensity) tuples to plot. The units of energy and intensity are irrelevant here (but should be consistent). Note that coordinates with 0 intensity will be removed.
        """
        # We save our coordinates under two properties.
        # Base coordinates are untransformed, coordinates (which is by default is the same as base_coordinates) are transformed.
        self.base_coordinates = [(x, y) for x,y in coordinates if y != 0]
        
    @classmethod
    def from_vibrations(self, vibrations):
        """
        Alternative constructor from a Vibration_list object.
        """
        return self([(vibration.frequency, vibration.intensity) for vibration in vibrations if vibration.intensity != 0])
    
    @property
    def coordinates(self):
        """
        Get the coordinates around which this graph is plotted.
        
        Note that this property may be transformed into different units, use base_coordinates for untransformed variant.
        """
        return self.base_coordinates
    
    def peaks(self, fwhm, resolution = 1, cutoff = 0.01):
        """
        Find peaks in the cumulative graph.
        
        :return: A list of x coords that are peaks.
        """
        coords = self.plot_cumulative_gaussian(fwhm = fwhm, resolution = resolution, cutoff = cutoff)
        y_coords = [y for x,y in coords]
        indexes = scipy.signal.find_peaks(y_coords, height = max(y_coords) * cutoff)[0]
        return [coords[index][0] for index in indexes]
        
    def plot_gaussian(self, fwhm, resolution = 1, cutoff = 0.01):
        """
        Plot a gaussian distribution around a set of coordinates.
        
        :param fwhm: The full-width at half-maximum of the gaussian function.
        :param cutoff: The minimum y value to plot using the gaussian function, as the fraction of the intensity.
        :param resolution: The spacing between points to plot using the gaussian function, in units of the x-axis.
        :return: A list of lists of tuples of (x, y) coordinates plotted by the gaussian function (one list per input coordinate).
        """
        # First, determine our c value.
        c = self.fwhm_to_c(fwhm)
        
        # Next, determine the limits in which we'll plot.
        limits = self.gaussian_limits(c, cutoff, resolution)
        
        # Plot and return.
        return [[(x, self.gaussian(a, b, c, x)) for x in numpy.linspace(*limits)] for b, a in self.base_coordinates]
    
    def plot_cumulative_gaussian(self, *args, **kwargs):
        """
        Plot an additive gaussian distribution around a set of coordinates.
        
        :param *args: The same as for plot_gaussian().
        :return: A single list of tuples of (x, y) coordinates plotted by the gaussian function.
        """
        # First, get our gaussian plots.
        gplots = self.plot_gaussian(*args, **kwargs)
        
        # Sum our y values.
        coords = {}
        for plot in gplots:
            for x, y in plot:
                coords[x] = coords.get(x, 0) + y
                
        # Now return as list of tuples.
        return list(coords.items())    
        
        
    def gaussian_limits(self, c, cutoff, resolution):
        """
        Determine min and max x limits to plot a gaussian function.
        
        :param c: The width of the peak.
        :param cutoff: The minimum y value to plot using the gaussian function, as the fraction of a.
        :param resolution: The spacing between points to plot using the gaussian function, in units of the x-axis.
        :return: A tuple of (minlim, maxlim, num) where minlim is the most negative value, maxlim the most positive and num the integer number of points to plot between them to achieve resolution.
        """
        # Calculate limits for each set of coordinates given to us.
        all_limits = [self.gaussian_x(y, x, c, cutoff * y) for x, y in self.base_coordinates]
        try:
            limits = (min(itertools.chain.from_iterable(all_limits)), max(itertools.chain.from_iterable(all_limits)))
        except ValueError:
            if len(list(itertools.chain.from_iterable(all_limits))) == 0:
                # Nothing to plot.
                raise Silico_exception("'{}' cannot plot spectrum; there are no values".format(type(self).__name__))
            else:
                raise
        
        # Now we need to generate the x values which we'll plot for.
        # This is a little more complicated than it needs to be because there's no simple range() for floats.
        # Calculate the number of points (and round up).
        num_points = round( math.fabs(limits[1] - limits[0]) / resolution)
        
        # Extend our limits so they are a clean multiple of num_points.
        limits = (limits[0], limits[0] + num_points * resolution)
        
        # And now shift our limits so they are still centred.
        xs = [x for x, y in self.base_coordinates]
        centre = (max(xs) - min(xs)) /2 + min(xs)
        limits = (centre - ((limits[1] - limits[0]) / 2) , centre + ((limits[1] - limits[0]) / 2))

        return (limits[0], limits[1], num_points +1)
    
    @classmethod
    def gaussian(self, a, b, c, x):
        """
        An implementation of the gaussian function.
        
        :param a: The maximum height of the peak.
        :param b: The x position of the center of the peak.
        :param c: The width of the peak.
        :param x: The x-value to plot for.
        :return: The corresponding y value.
        """
        return a * math.exp(-( (x - b)**2 / (2 * c**2) ))
    
    @classmethod
    def gaussian_x(self, a, b, c, y):
        """
        An implementation of the gaussian function rearranged to find x.
        
        :param a: The maximum height of the peak.
        :param b: The x position of the center of the peak.
        :param c: The width of the peak.
        :param y: The y-value to plot for.
        :return: A tuple of the two corresponding x values.
        """
        # We have two solutions of the form x = ± d + b
        # Calculate d first.
        try:
            d = math.sqrt( -math.log( y/a ) * 2 * c **2  )
        except ZeroDivisionError:
            # a (max height) is zero; the intensity is zero.
            # We could instead return b?
            raise Silico_exception("'{}' cannot calculate Gaussian function limits; a (intensity) is zero".format(self.__name__))
        
        # Now return our two solutions.
        return (-d + b, d + b)
        
        
    @classmethod
    def fwhm_to_c(self, fwhm):
        """
        Convert a full width at half maximum to the corresponding c-value in the Gaussian function.
        
        :param fwhm: The desired full width at half maximum (in nm).
        :return: The c-value.
        """
        return fwhm / (2 * math.sqrt(2 * math.log(2)))
    
    
class Absorption_emission_graph(Spectroscopy_graph):
    """
    Class for graphing absorption/emission spectra.
    """
    
    def __init__(self, coordinates, *, use_jacobian = True):
        """
        Constructor for Absorption_emission_graph objects.
        
        :param use_jacobian: Whether to use the jacobian transform to scale the y axis.
        """
        # Get our x,y values from our excited states.
        super().__init__(coordinates)
        self.use_jacobian = use_jacobian
        
    @property
    def coordinates(self):
        """
        Get the coordinates around which this graph is plotted.
        
        Note that this property may be transformed into different units, use base_coordinates for untransformed variant.
        """
        return self.nm_coordinates
        
    @classmethod
    def from_excited_states(self, excited_states, **kwargs):
        """
        An alternative constructor that takes a list of excited states as argument.
        
        :param excited_states: An Excited_states_list object to construct from.
        :param **kwargs: Passed to the real constructor.
        """
        return self([(excited_state.energy, excited_state.oscillator_strength if excited_state.oscillator_strength is not None else 0) for excited_state in excited_states], **kwargs)
    
    @property
    def nm_coordinates(self):
        """
        Get the list of coordinates of this graph scaled to nm.
        """
        return [self.energy_to_wavelength(coord, self.use_jacobian) for coord in self.base_coordinates]
    
    @classmethod
    def energy_to_wavelength(self, coord, use_jacobian):
        """
        Convert a pair of x,y coordinates in ev to nm.
        
        :param coord: A tuple of (eV, f).
        :return: A tuple of (nm, i), where i is intensity.
        """
        x, y = coord
        return (silico.result.excited_states.Excited_state.energy_to_wavelength(x), self.jacobian(x, y) if use_jacobian else y)
        
    @classmethod
    def jacobian(self, E, f_E):
        """
        An implementation of the jacobian transform that scales intensity in energy units to intensity in wavelength units.
        
        See J. Phys. Chem. Lett. 2014, 5, 20, 3497 for why this is necessary.
        
        Note that the jacobian transform will maintain the area under the curve regardless of x units (nm or x).
        Sadly, this has the consequence of mangling the intensity units (it becomes tiny; an oscillator strength of 1 at 3 eV becomes 1.163e-12).
        """
        return ((E * scipy.constants.electron_volt)**2  * f_E) / (scipy.constants.Planck * scipy.constants.c)
    
    def plot_gaussian(self, *args, **kwargs):
        """
        Plot a gaussian distribution around our excited state energies.
        
        :param fwhm: The full-width at half-maximum of the gaussian function.
        :param cutoff: The minimum y value to plot using the gaussian function, as the fraction of the intensity.
        :param resolution: The spacing between points to plot using the gaussian function, in units of the x-axis.
        :return: A list of lists of tuples of (x, y) coordinates plotted by the gaussian function (one list per input coordinate).
        """
        # All we need to do over our parent is convert x values from e to wavelength.
        # And scale y values using the jacobian transform.
        return [[self.energy_to_wavelength(coord, self.use_jacobian) for coord in plot] for plot in super().plot_gaussian(*args, **kwargs)]
    
    
    
    
    