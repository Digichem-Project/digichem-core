from silico.image.graph import Graph_image_maker
from matplotlib import ticker
import matplotlib.collections
import math
import silico.result.excited_states
from silico.exception.base import File_maker_exception, Silico_exception

class Spectroscopy_graph_maker(Graph_image_maker):
	"""
	A graph type for displaying spectroscopy data (UV-vis, IR etc.).
	"""
	def __init__(
		self,
		output,
		x_values,
		y_values,
		*, # Keyword only arguments from here:
		x_padding = 60,
		fwhm = 40,
		plot_bars = True,
		plot_peaks = False,
		plot_cumulative_peak = True,
		x_limits_method = 'standard',
		y_limits_method = 'standard',
		**kwargs		
	):
		"""
		Constructor for uv-vis type graphs.
		
		:param x_values: List of x_values to plot.
		:param y_values: List of y_values to plot.
		:param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
		:param x_padding: The width that we will add on to either end of our x-axis.
		:param fwhm: The desired full width at half maximum of the peaks that we are going to plot.
		:param plot_bars: If True, vertical bars are plotted on the graph for each excited state.
		:param plot_peaks: If True, Gaussian peaks are plotted for each excited state.
		:parma plot_cumulative_peak: If True, a cumulative graph is plotted from gaussian peaks from all excited states.
		:param x_limits_method: String controlling how the x axis limits are set. Options are 'standard' for standard scaling, showing all plotted peaks. Alternatively, a tuple of (x_min, x_max) which will be used directly as axis limits.
		:param y_limits_method: String controlling how the y axis limits are set. Options are 'standard' for automatic scaling. Alternatively, a tuple of (y_min, y_max) which will be used directly as axis limits.
		
		"""
		super().__init__(output, **kwargs)
		
		self.x_values = x_values
		self.y_values = y_values
		
		# Set some config stuff that will be changeable later.
		# The width (in nm) that we add on to our x axis.
		self.x_padding = x_padding + round(fwhm) + 40
		# The 'c' value in the Gaussian function; decides how wide our peaks are.
		self.fwhm = fwhm
		
		# Options controlling what to plot.
		self.should_plot_bars = plot_bars
		self.should_plot_peaks = plot_peaks
		self.should_plot_cumulative_peak = plot_cumulative_peak
		
		self.init_image_width = 10
		self.init_image_height = 5.26
		
		# Axis scaling.
		self.x_limits_method = x_limits_method
		self.y_limits_method = y_limits_method
	
	@property
	def c_value(self):
		return self.fwhm_to_c(self.fwhm)
	
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
		columns = []
		for coord in zip(self.x_values, self.y_values):
			columns.append([
					(coord[0], 0),
					(coord[0], coord[1])
			])
		self.axes.add_collection(matplotlib.collections.LineCollection(columns, colors = (0,0,0)))
		
	def plot_peaks(self):
		"""
		Plot x and y values as individual peaks on the graph we are building.
		
		:param x_values: The x-values that we are going to plot.
		:param y_list: 2D list of y-values, one list for each excited state.
		"""
		for line_y_values in self.line_y_values:
			# Plot as a line.
			self.axes.plot(self.line_x_values, line_y_values, 'C0-', linewidth = 1)
			
	def plot_cumulative_peak(self):
		"""
		Plot x and y values as a single plot on the graph we are building.
		
		:param x_values: The x-values that we are going to plot.
		:param y_list: 2D list of y-values, one list for each excited state.
		"""
		# Sum all our y-values together.
		line_y_values = self.combined_line_y_values
		# Now plot.
		self.axes.plot(self.line_x_values, line_y_values, 'C0-', linewidth = 2)
		
	def plot_hidden_points(self):
		"""
		Plot invisible points to take advantage of matplotlib auto scaling.
		"""		
		# Plot. Set size (s) to zero so we don't see our points, but we still get good automatic scaling.
		self.axes.scatter(self.x_values, self.y_values, s=0)
	
	@property
	def line_x_values(self):
		"""
		Get the x values that we're going to use to plot our line. This will be the same for all peaks.
		
		:return: The x values as a list.
		"""
		# First decide how many points we're going to plot.
		min_value = max(math.floor(min(self.x_values)) - self.x_padding, 1)
		max_value = math.ceil(max(self.x_values)) + self.x_padding
		
		# Now get our x values, which are just integers from min_value to max_value.
		return list(range(min_value, max_value +1))
	
	@property
	def combined_line_y_values(self):
		"""
		The Y values for our line of all peaks.
		
		Unlike line_y_values, combined_y_values is a 1D list.
		"""
		return list(map(sum, zip(*self.line_y_values)))
	
	@property
	def line_y_values(self):
		"""
		Get the y values to plot a line/peak from a point
		
		:param line_x_values: The x_values (wavelength in nm) to plot our line over.
		:param point_x_value: The x coordinate of the point that we are turning into a peak.
		:param point_y_value: The y coordinate of the point that we are turning into a peak.
		:return: The line_y_values.
		"""
		# A big list of x coordinates that we'll plot our line along. 
		line_x_values = self.line_x_values
		
		# A list (per signal) of y coordinates.
		line_y_list = []
		
		# Go through each (x,y) coordinate (each of which represents that maxima of a peak) and plot a Gaussian curve around it.
		for point_x_value, point_y_value in zip(self.x_values, self.y_values):
			# The (x,y) coordinates of the points that make up the line of this peak.
			line_y_values = []
			
			for line_x_value in line_x_values:
				line_y_values.append(self.gaussian(
					a = point_y_value,
					b = point_x_value,
					c = self.c_value,
					x = line_x_value))
				
			# Add to our list of lists.
			line_y_list.append(line_y_values)
		
		# All done.		
		return line_y_list
	
	def adjust_axes(self):
		"""
		Adjust the axes of our graph.
		
		This method is called automatically as part of the make() method.
		"""
		super().adjust_axes()
		
		# Call the appropriate scaling method depending on the value of the limits_method.
		# X axis scaling depends on what option's been set.
		if self.x_limits_method == "standard":
			self.standard_x_limits()
		elif isinstance(self.x_limits_method, tuple) or isinstance(self.x_limits_method, list):
			self.simple_x_limits(self.x_limits_method[0], self.x_limits_method[1])
		else:
			raise ValueError("Unknown x limits method '{}'".format(self.x_limits_method))
		
		# Y axis scaling depends on what option's been set.
		if self.y_limits_method == "standard":
			self.standard_y_limits()
		elif isinstance(self.y_limits_method, tuple) or isinstance(self.y_limits_method, list):
			self.simple_y_limits(self.y_limits_method[0], self.y_limits_method[1])
		else:
			raise ValueError("Unknown y limits method '{}'".format(self.y_limits_method))
	
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
	def fwhm_to_c(self, fwhm):
		"""
		Convert a full width at half maximum to the corresponding c-value in the Gaussian function.
		
		:param fwhm: The desired full width at half maximum (in nm).
		:return: The c-value.
		"""
		return fwhm / (2 * math.sqrt(2 * math.log(2)))
		

class Absorption_graph_maker(Spectroscopy_graph_maker):
	"""
	A graph for displaying simulated absorptions.
	"""
	
	def __init__(self, output, excited_states, peak_cutoff = 0.05, max_width = None, *args, **kwargs):
		"""
		Constructor for UV-Vis type absorption graphs.
		
		:param output: A path to an output file to write to. The extension of this path is used to determine the format of the file (eg, png, jpeg).
		:param excited_states: List of excited states to plot.
		:param peak_cutoff: The minimum oscillator strength a peak must have to be drawn, as a fraction of the tallest peak. Set to 0 for no cutoff.
		:param max_width: The maximum width in pixels of the graph.
		"""
		# Now go through our excited states and build our data.
		x_values = []
		y_values = []
		
		# Remove states with no oscillator strength (we can't plot them anyway).
		excited_states = [excited_state for excited_state in excited_states if excited_state.oscillator_strength is not None and excited_state.oscillator_strength > 0]
		
		for excited_state in excited_states:
			x_values.append(excited_state.wavelength)
			y_values.append(excited_state.oscillator_strength if excited_state.oscillator_strength is not None else 0)
		
		# Call our parent.
		super().__init__(output, x_values, y_values, *args, **kwargs)
		
		if len(excited_states) == 0:
			# We have nothing to plot; we could plot an empty graph (an maybe should?), but for now we'll stop.
			raise File_maker_exception(self, "cannot plot spectrum; there are no states with oscillator strength above 0")
		
		# The amount of space (in matplotlib inches) to allocate per unit of the x axis.
		self.inch_per_x = 0.0192
		# Our max image width in pixels.
		self.max_width = max_width
		# The min oscillator strength to plot as a percentage.
		self.peak_cutoff = peak_cutoff
		
		# The DPI to save our image with.
		self.output_dpi = 100
		
		# Axes labels.
		self.x_label = "Wavelength /nm"
		self.y_label = "Oscillator Strength"
		
	@classmethod
	def from_image_options(self, output, *, excited_states, output_base = None, options, **kwargs):
		"""
		An alternative constructor that discards any additional keyword arguments.
		"""
		return self(
			output,
			excited_states = excited_states,
			x_limits_method = options['absorption_graph']['x_limits'],
			y_limits_method = options['absorption_graph']['y_limits'],
			max_width = options['absorption_graph']['max_width'],
			peak_cutoff = options['absorption_graph']['peak_cutoff'],
			fwhm = options['absorption_graph']['fwhm'],
			x_padding = options['absorption_graph']['x_padding'],
			output_base = output_base,
			dont_modify = options['image']['dont_create_new'],
			use_existing = options['image']['use_existing']
		)

		
	
# 	@property
# 	def line_y_values(self):
# 		"""
# 		Get the y values to plot a line/peak from a point
# 		
# 		:param line_x_values: The x_values (wavelength in nm) to plot our line over.
# 		:param point_x_value: The x coordinate of the point that we are turning into a peak.
# 		:param point_y_value: The y coordinate of the point that we are turning into a peak.
# 		:return: The line_y_values.
# 		"""
# 		# A big list of x coordinates that we'll plot our line along. 
# 		line_x_values = self.line_x_values
# 		
# 		# A list (per signal) of y coordinates.
# 		line_y_list = []
# 		
# 		# Go through each (x,y) coordinate (each of which represents that maxima of a peak) and plot a Gaussian curve around it.
# 		for point_x_value, point_y_value in zip(self.x_values, self.y_values):
# 			# The (x,y) coordinates of the points that make up the line of this peak.
# 			line_y_values = []
# 			
# 			for line_x_value in line_x_values:
# 				line_y_values.append(self.gaussian(
# 					a = point_y_value,
# 					b = silico.result.excited_states.Excited_state.emission_wavelength_to_energy(point_x_value),
# 					c = self.c_value,
# 					x = silico.result.excited_states.Excited_state.emission_wavelength_to_energy(line_x_value)))
# 				
# 			# Add to our list of lists.
# 			line_y_list.append(line_y_values)
# 		
# 		# All done.		
# 		return line_y_list
	
	def standard_x_limits(self):
		"""
		Limit the X axis so all plotted peaks are visible.
		"""
		# We need to get a list of all peaks that are above our cutoff point.
		# First determine our highest point.
		highest_point = max(self.combined_line_y_values)
		# Now filter by a fraction of that amount.
		visible_x_values = [x for x, y in zip(self.line_x_values, self.combined_line_y_values) if y >= (highest_point * self.peak_cutoff)]
		# Now we can just get our min-max, remembering to include our x-padding.
		self.axes.set_xlim(min(visible_x_values) - self.x_padding, max(visible_x_values) + self.x_padding)
		
	def standard_y_limits(self):
		"""
		Limit the y axis so all plotted peaks are visible.
		"""
		self.axes.autoscale(enable = True, axis = 'y')
		# Clamp to 0 -> pos.
		self.axes.set_ylim(0, self.axes.get_ylim()[1])
		
	def simple_x_limits(self, x_min = 250, x_max = 700):
		"""
		Limit the X axis between between two points.
		
		:param x_min: The start of the y axis.
		:param x_max: The end of the y axis.
		"""
		self.axes.set_xlim(x_min, x_max)
		
	def simple_y_limits(self, y_min = 0, y_max = 1):
		"""
		Limit the X axis between between two points.
		
		:param x_min: The start of the y axis.
		:param x_max: The end of the y axis.
		"""
		self.axes.set_ylim(y_min, y_max)
	
	def adjust_axes(self):
		"""
		Adjust the axes of our graph.
		
		This method is called automatically as part of the make() method.
		"""
		super().adjust_axes()
		
		# Change spacing of tick markers.
		self.axes.xaxis.set_major_locator(ticker.MultipleLocator(50))
		self.axes.xaxis.set_minor_locator(ticker.MultipleLocator(10))

		
		
		# Use matplotlib to automatically layout our graph.
		self.figure.tight_layout()
		
		# Get a constant scale on our x axis.
		self.constant_scale(0, self.inch_per_x)
		
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
		x_values = []
		y_values = []
		for vibration in vibrations:
			x_values.append(vibration.frequency)
			y_values.append(vibration.intensity)
			
		# Call our parent.
		super().__init__(output, x_values, y_values, *args, **kwargs)
		
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

	def standard_x_limits(self):
		"""
		Limit the X axis so all plotted peaks are visible.
		"""
		# We use automatic X scaling, except we don't go past 0, we show at least 3500 cm-1 and our scale is reversed (this is typical for IR).
		self.axes.autoscale(enable = True, axis = 'x')
		self.axes.set_xlim(max(self.axes.get_xlim()[1], 3500), 0)
		
	def standard_y_limits(self):
		"""
		Limit the y axis so all plotted peaks are visible.
		"""
		# We use automatic Y scaling, except we don't go past 0 and we are reversed (this is typical for IR).
		self.axes.autoscale(enable = True, axis = 'y')
		self.axes.set_ylim(self.axes.get_ylim()[1], 0)
		
	def simple_x_limits(self, x_min = 3500, x_max = 0):
		"""
		Limit the X axis between between two points.
		
		:param x_min: The start of the y axis.
		:param x_max: The end of the y axis.
		"""
		self.axes.set_xlim(x_min, x_max)
		
	def simple_y_limits(self, y_min = 1, y_max = 0):
		"""
		Limit the X axis between between two points.
		
		:param x_min: The start of the y axis.
		:param x_max: The end of the y axis.
		"""
		self.axes.set_ylim(y_min, y_max)

	def adjust_axes(self):
		"""
		Adjust the axes of our graph.
		
		This method is called automatically as part of the make() method.
		"""
		super().adjust_axes()
		
		# Adjust our axes limits.
		#self.axes.set_xlim(max(self.axes.get_xlim()[1], 3500), 0)
		#self.axes.set_ylim(self.axes.get_ylim()[1], 0)
		
		# Use matplotlib to automatically layout our graph.
		self.figure.tight_layout()
		
		# Get a constant scale on our x axis.
		self.constant_scale(0, self.inch_per_x)
		
		
		