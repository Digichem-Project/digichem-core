from silico.image.graph import Graph_image_maker
from matplotlib import ticker
import matplotlib.collections
from silico.exception.base import File_maker_exception
from silico.result.spectroscopy import Spectroscopy_graph,\
    Absorption_emission_graph

class Spectroscopy_graph_maker(Graph_image_maker):
    """
    A graph type for displaying spectroscopy data (UV-vis, IR etc.).
    """
    def __init__(
        self,
        output,
        graph,
        *,
        # Keyword only arguments from here:
        # Gaussian options.
        fwhm = 0.4,
        gaussian_cutoff = 0.001,
        gaussian_resolution = 0.01,
        
        # Options that control what to show in the graph.
        plot_bars = True,
        plot_peaks = False,
        plot_cumulative_peak = True,
        
        # Options controlling the axes limits.
        x_padding = 50,
        peak_cutoff = 0.01,
        x_limits = 'auto',
        y_limits = 'auto',
        **kwargs
    ):
        """
        Constructor for spectroscopy graphs.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param graph: An instance of a silico.result.spectroscopy.Spectroscopy_graph that contains data to plot.
        :param fwhm: The desired full width at half maximum of the peaks that we are going to plot.
        :param gaussian_cutoff: The minimum y value to plot using the gaussian function, as the fraction of the max intensity.
        :param gaussian_resolution: The spacing between points to plot using the gaussian function, in units of the x-axis.
        :param plot_bars: If True, vertical bars are plotted on the graph for each excited state.
        :param plot_peaks: If True, Gaussian peaks are plotted for each excited state.
        :param plot_cumulative_peak: If True, a cumulative graph is plotted from gaussian peaks from all excited states.
        :param x_padding: An amount, in units of the x-axis, to extend beyond the lowest & highest energy peak (so the graph is framed). Note that the graph will never extend beyond 0.
        :param y_cutoff: The minimum y value to include when using auto x limits (the x axis will extend
        :param x_limits_method: String controlling how the x axis limits are set. Options are 'auto' for automatic scaling, showing all plotted peaks. Alternatively, a tuple of (x_min, x_max) which will be used directly as axis limits.
        :param y_limits_method: String controlling how the y axis limits are set. Options are 'auto' for automatic scaling. Alternatively, a tuple of (y_min, y_max) which will be used directly as axis limits.
        
        """
        # Call our image maker parent.
        super().__init__(output, **kwargs)
        
        # Save our graph object (holds coordinates to plot).
        self.graph = graph
        
        # Controls for plotting Gaussian's.
        self.fwhm = fwhm
        self.gaussian_cutoff = gaussian_cutoff
        self.gaussian_resolution = gaussian_resolution
        
        # Options controlling what to plot.
        self.should_plot_bars = plot_bars
        self.should_plot_peaks = plot_peaks
        self.should_plot_cumulative_peak = plot_cumulative_peak
        
        self.init_image_width = 10
        self.init_image_height = 5.26
        
        # Axis scaling.
        self.x_padding = x_padding
        self.peak_cutoff = peak_cutoff
        self.x_limits_method = x_limits
        self.y_limits_method = y_limits
        
    @classmethod
    def transpose(self, coordinates):
        """
        Transpose a list of coordinates from [(x1,y1), (x2,y2)...] to [[x1, x2...], [y1, y2...]]
        """
        return list(map(list, zip(*coordinates)))
    
    def plot(self):
        """
        Plot the contents of our graph.
        """
        # Make sure we actually have something to plot.
        if len(self.graph.coordinates) == 0:
            # We have nothing to plot; we could plot an empty graph (and maybe should?), but for now we'll stop.
            raise File_maker_exception(self, "cannot plot spectrum; there are no values with intensity above 0")
        
        # What we do next depends on what options we've been given.
        if self.should_plot_bars:
            self.plot_columns()
            # If we aren't going to plot any curves, plot some invisible points so we get automatic axis scaling.
            if not self.should_plot_peaks and not self.should_plot_cumulative_peak:
                self.plot_hidden_points()
                
        if self.should_plot_peaks:
            self.plot_lines()
            
        if self.should_plot_cumulative_peak:
            self.plot_cumulative_line()
                        
    def plot_columns(self):
        """
        Plot vertical columns for each plot on the graph we are building.
        """
        columns = [[(x, 0) , (x, y)] for x, y in self.graph.coordinates]
        self.axes.add_collection(matplotlib.collections.LineCollection(columns, colors = (0,0,0)))
        
    def plot_hidden_points(self):
        """
        Plot invisible points to take advantage of matplotlib auto scaling.
        """        
        # Plot. Set size (s) to zero so we don't see our points, but we still get good automatic scaling.
        self.axes.scatter(*self.transpose(self.graph.coordinates), s=0)
        
    def plot_lines(self):
        """
        Plot x and y values as individual lines on the graph we are building.
        """            
        for plot in self.graph.plot_gaussian(self.fwhm, self.gaussian_resolution, self.gaussian_cutoff):
            data = self.transpose(plot)
            self.axes.plot(*data, 'C0-', linewidth = 1)
            
    def plot_cumulative_line(self):
        """
        Plot x and y values as a single line on the graph we are building.
        """
        data = self.transpose(self.graph.plot_cumulative_gaussian(self.fwhm, self.gaussian_resolution, self.gaussian_cutoff))
        self.axes.plot(*data, 'C0-', linewidth = 2)
    
    @property
    def peaks(self):
        """
        A list of peaks (in x units).
        """
        return [x for x,y in self.graph.peaks(self.fwhm, self.gaussian_resolution, self.gaussian_cutoff)]
    
    def selected_peaks(self, decimals = 0, number = None):
        """
        A unique list of sorted peaks, rounded to a given number of decimal points.
        
        :param decimals: The number of decimal points to round to.
        :param number: The number of peaks to return (the most intense peaks will be returned first).
        :returns: A list of ordered peaks as floats (if decimals > 0) or ints.
        """
        if number is not None:
            peaks = self.graph.peaks(self.fwhm, self.gaussian_resolution, self.gaussian_cutoff)
            peaks.sort(key = lambda coord: coord[1])
            peaks = [x for x,y in peaks[-number:]]
        
        else:
            peaks = self.peaks
        
        return sorted(list(set([round(peak, decimals) if decimals > 0 else int(peak) for peak in peaks])))
        
    def auto_x_limits(self):
        """
        Limit the X axis so all plotted peaks are visible.
        """
        # We need to get a list of all peaks that are above our cutoff point.
        # First determine our highest point.
        highest_point = max(self.transpose(self.graph.coordinates)[1])
        
        # Now filter by a fraction of that amount.
        visible_x_values = [x for x, y in self.graph.plot_cumulative_gaussian(self.fwhm, self.gaussian_resolution, self.gaussian_cutoff) if y >= (highest_point * self.peak_cutoff)]
        
        # Now we can just get our min-max, remembering to include our x-padding, and never extending beyond 0 (because negative energy has no meaning in this context(?)).
        self.axes.set_xlim(max(min(visible_x_values) - self.x_padding, 0), max(visible_x_values) + self.x_padding)
        
    def auto_y_limits(self):
        """
        Limit the y axis so all plotted peaks are visible.
        """
        # Use matplotlib autoscaling.
        self.axes.autoscale(enable = True, axis = 'y')
        
        # Clamp to 0 -> pos.
        self.axes.set_ylim(0, self.axes.get_ylim()[1])
            
    def adjust_axes(self):
        """
        Adjust the axes of our graph.
        
        This method is called automatically as part of the make() method.
        """
        super().adjust_axes()
        
        # Call the appropriate scaling method depending on the value of the limits_method.
        for axis in ("x", "y"):
            method = getattr(self, axis + "_limits_method")
            if method == "auto":
                # Call the automatic function.
                getattr(self, "auto_{}_limits".format(axis))()
            elif isinstance(method, tuple) or isinstance(method, list):
                # Manual, explicit limits.
                getattr(self.axes, "set_{}lim".format(axis))(method[0], method[1])
            else:
                # Don't recognise method.
                raise ValueError("Unknown {} limits method '{}'".format(axis, method))
            
        # Use matplotlib to automatically layout our graph.
        self.figure.tight_layout()
        
        # Get a constant scale on our x axis.
        self.constant_scale(0, self.inch_per_x)
        

class Absorption_emission_graph_maker(Spectroscopy_graph_maker):
    """
    A graph for displaying simulated absorptions.
    """
    
    # The key/name of where our options are stored in the main config.
    options_name = None
    
    def __init__(self, output, excited_states, *args, adjust_zero = False, use_jacobian = True, **kwargs):
        """
        Constructor for UV-Vis type absorption graphs.
        
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param excited_states: List of excited states to plot.
        :param adjust_zero: If all the intensities of the given excited states are zero, whether to arbitrarily set the y coords to 1.
        :param peak_cutoff: The minimum oscillator strength a peak must have to be drawn, as a fraction of the tallest peak. Set to 0 for no cutoff.
        :param max_width: The maximum width in pixels of the graph.
        """        
        # Call our parent.
        super().__init__(output, Absorption_emission_graph.from_excited_states(excited_states, adjust_zero = adjust_zero, use_jacobian = use_jacobian), *args, **kwargs)
        
        # The amount of space (in matplotlib inches) to allocate per unit of the x axis.
        self.inch_per_x = 0.0192
        
        # The DPI to save our image with.
        self.output_dpi = 100
        
        # Axes labels.
        self.x_label = "Wavelength /nm"
        if not self.graph.false_intensity:
            self.y_label = "Intensity" if self.graph.use_jacobian else "Oscillator Strength"
        else:
            self.hide_y = True
    
    @classmethod
    def from_options(self, output, *, excited_states, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """    
        return self(
            output,
            excited_states = excited_states,
            **options[self.options_name],
            **kwargs
        )
    
    def adjust_axes(self):
        """
        Adjust the axes of our graph.
        
        This method is called automatically as part of the make() method.
        """
        # Change spacing of tick markers.
        self.axes.xaxis.set_major_locator(ticker.MultipleLocator(50))
        self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(10))

        # Call parent.
        super().adjust_axes()
        
        # Check to see if we've exceeded our max width.
        if self.max_width is not None and (self.figure.get_size_inches()[0] * self.output_dpi) > self.max_width:
            # Adjust our tick labels.
            x_difference = abs(self.axes.get_xlim()[0] - self.axes.get_xlim()[1])
            if x_difference < 1200:
                self.axes.xaxis.set_major_locator(ticker.MultipleLocator(100))
                self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(50))
            else:
                self.axes.xaxis.set_major_locator(ticker.MultipleLocator(200))
                self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(100))
            
            # Update our layout.
            self.figure.tight_layout()
            
class Absorption_graph_maker(Absorption_emission_graph_maker):
    """
    Class for creating absorption (UV-Vis) spectra.
    """
    options_name = "absorption_spectrum"
    
class Emission_graph_maker(Absorption_emission_graph_maker):
    """
    Class for creating emission spectra.
    """
    options_name = "emission_spectrum"

class Frequency_graph_maker(Spectroscopy_graph_maker):
    """
    A graph for displaying vibrational frequencies.
    """
            
    def __init__(self, output, vibrations, *args, fwhm = 80, gaussian_resolution = 1, **kwargs):
        """
        Constructor for frequency graphs.
        
        :param vibrations: List of vibrations that we're going to plot.
        :param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
        :param x_limits_method: String controlling how the x axis limits are set. Options are 'auto' for standard auto scaling, showing all plotted peaks. Alternatively, a tuple of (x_min, x_max) which will be used directly as axis limits.
        :param y_limits_method: String controlling how the y axis limits are set. Options are 'auto' for standard auto scaling. Alternatively, a tuple of (y_min, y_max) which will be used directly as axis limits.
        """
        # Call our parent.
        super().__init__(output, Spectroscopy_graph.from_vibrations(vibrations), *args, fwhm = fwhm, gaussian_resolution = gaussian_resolution, **kwargs)
        
        # Axes labels.
        self.x_label = "Frequency /cm$^{-1}$"
        self.y_label = "Intensity /km mol$^{-1}$"
        
        # Space to allocate per unit of the y axis.
        self.inch_per_x = 0.00252
        
        
    @classmethod
    def from_options(self, output, *, vibrations, options, **kwargs):
        """
        Constructor that takes a dictionary of config like options.
        """        
        return self(
            output,
            vibrations = vibrations,
            **options['IR_spectrum'],
            **kwargs
        )

    def auto_x_limits(self):
        """
        Limit the X axis so all plotted peaks are visible.
        """
        # We use automatic X scaling, except we don't go past 0, we show at least 3500 cm-1 and our scale is reversed (this is typical for IR).
        self.axes.autoscale(enable = True, axis = 'x')
        self.axes.set_xlim(max(self.axes.get_xlim()[1], 3500), 0)
        
    def auto_y_limits(self):
        """
        Limit the y axis so all plotted peaks are visible.
        """
        # We use automatic Y scaling, except we don't go past 0 and we are reversed (this is typical for IR).
        self.axes.autoscale(enable = True, axis = 'y')
        self.axes.set_ylim(self.axes.get_ylim()[1], 0)
        
        
        