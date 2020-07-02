from silico.image.graph import Graph_image_maker
from matplotlib import ticker
import matplotlib.collections
from silico.exception.base import File_maker_exception
from silico.result.spectroscopy import Spectroscopy_graph,\
	Absorption_emission_graph

class Spectroscopy_graph_maker(Graph_image_maker, Spectroscopy_graph):
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
		fwhm = 40,
		gaussian_cutoff: 0.01,
		gaussian_resolution: 0.001,
		
		# Options that control what to show in the graph.
		plot_bars = True,
		plot_peaks = False,
		plot_cumulative_peak = True,
		
		# Options controlling the axes limits.
		x_padding = 40,
		peak_cutoff = 0.01,
		x_limits_method = 'standard',
		y_limits_method = 'standard',
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
		
		# Make sure we actually have something to plot.
		if len(self.graph.coordinates) == 0:
			# We have nothing to plot; we could plot an empty graph (and maybe should?), but for now we'll stop.
			raise File_maker_exception(self, "cannot plot spectrum; there are no values with intensity above 0")
		
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
		self.x_limits_method = x_limits_method
		self.y_limits_method = y_limits_method
		
		
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
		# What we do next depends on what options we've been given.
		if self.should_plot_bars:
			self.plot_columns()
			# If we aren't going to plot any curves, plot some invisible points so we get automatic axis scaling.
			if not self.should_plot_peaks and not self.should_plot_cumulative_peak:
				self.plot_hidden_points()
				
		if self.should_plot_peaks:
			self.plot_peaks()
			
		if self.should_plot_cumulative_peak:
			self.plot_cumulative_peak()
						
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
		
	def auto_x_limits(self):
		"""
		Limit the X axis so all plotted peaks are visible.
		"""
		# We need to get a list of all peaks that are above our cutoff point.
		# First determine our highest point.
		highest_point = max(self.transpose(self.graph.coordinates)[1])
		
		# Now filter by a fraction of that amount.
		visible_x_values = [x for x, y in self.plot_cumulative_gaussian() if y >= (highest_point * self.peak_cutoff)]
		
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
		

class Absorption_graph_maker(Spectroscopy_graph_maker):
	"""
	A graph for displaying simulated absorptions.
	"""
	
	def __init__(self, output, excited_states, max_width = None, *args, **kwargs):
		"""
		Constructor for UV-Vis type absorption graphs.
		
		:param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
		:param excited_states: List of excited states to plot.
		:param peak_cutoff: The minimum oscillator strength a peak must have to be drawn, as a fraction of the tallest peak. Set to 0 for no cutoff.
		:param max_width: The maximum width in pixels of the graph.
		"""		
		# Call our parent.
		super().__init__(output, Absorption_emission_graph.from_excited_states(excited_states), *args, **kwargs)
		
		# The amount of space (in matplotlib inches) to allocate per unit of the x axis.
		self.inch_per_x = 0.0192
		# Our max image width in pixels.
		self.max_width = max_width
		
		# The DPI to save our image with.
		self.output_dpi = 100
		
		# Axes labels.
		self.x_label = "Wavelength /nm"
		self.y_label = "Intensity" if self.graph.use_jacobian else "Oscillator Strength"
		
	@classmethod
	def from_image_options(self, output, *, excited_states, output_base = None, options, **kwargs):
		"""
		An alternative constructor that discards any additional keyword arguments.
		"""
		return self(
			output,
			excited_states = excited_states,
			max_width = options['absorption_graph']['max_width'],
			
			fwhm = options['absorption_graph']['fwhm'],
			gaussian_cutoff = ['absorption_graph']['gaussian_cutoff'],
			gaussian_resolution = ['absorption_graph']['gaussian_resolution'],
			
			x_padding = options['absorption_graph']['x_padding'],
			peak_cutoff = options['absorption_graph']['peak_cutoff'],
			x_limits_method = options['absorption_graph']['x_limits'],
			y_limits_method = options['absorption_graph']['y_limits'],
			
			output_base = output_base,
			dont_modify = options['image']['dont_create_new'],
			use_existing = options['image']['use_existing']
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
		fig_size = self.figure.get_size_inches()
		if (fig_size[0] * self.output_dpi) > self.max_width:
			# We're too big, reduce to max_width.
			fig_size[0] = self.max_width / self.output_dpi
			self.figure.set_size_inches(fig_size)
			
			# Also adjust our tick labels.
			x_difference = abs(self.axes.get_xlim()[0] - self.axes.get_xlim()[1])
			if x_difference < 1200:
				self.axes.xaxis.set_major_locator(ticker.MultipleLocator(100))
				self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(50))
			else:
				self.axes.xaxis.set_major_locator(ticker.MultipleLocator(200))
				self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(100))
			
			# Update layout. Calling this twice is necessary for some reason loool
			self.figure.tight_layout()
			self.figure.tight_layout()
			

class Frequency_graph_maker(Spectroscopy_graph_maker):
	"""
	A graph for displaying vibrational frequencies.
	"""
			
	def __init__(self, output, vibrations, *args, **kwargs):
		"""
		Constructor for frequency graphs.
		
		:param vibrations: List of vibrations that we're going to plot.
		:param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
		:param x_limits_method: String controlling how the x axis limits are set. Options are 'standard' for standard scaling, showing all plotted peaks. Alternatively, a tuple of (x_min, x_max) which will be used directly as axis limits.
		:param y_limits_method: String controlling how the y axis limits are set. Options are 'standard' for automatic scaling. Alternatively, a tuple of (y_min, y_max) which will be used directly as axis limits.
		"""
		# Call our parent.
		super().__init__(output, Spectroscopy_graph.from_vibrations(vibrations), *args, **kwargs)
		
		# Axes labels.
		self.x_label = "Frequency /cm-1"
		self.y_label = "Intensity /km mol-1"
		
		# Space to allocate per unit of the y axis.
		self.inch_per_x = 0.00252
		
		
	@classmethod
	def from_image_options(self, output, *, vibrations, output_base = None, options, **kwargs):
		"""
		An alternative constructor that discards any additional keyword arguments.
		"""
		return self(
			output,
			vibrations = vibrations,
			x_limits_method = options['frequency_graph']['x_limits'],
			y_limits_method = options['frequency_graph']['y_limits'],
			fwhm = options['frequency_graph']['fwhm'],
			output_base = output_base,
			dont_modify = options['image']['dont_create_new'],
			use_existing = options['image']['use_existing']
			
			
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
		
		
		