import textwrap
from itertools import chain


class Node_printer():
	"""
	Class for pretty printing a nested list of options (as used by csubmit -l).
	"""
	
	def __init__(self, options, *, indent = 5, width = 80, space = False):
		"""
		Constructor for Node_printer objects.
		
		:param options: A list of 3-membered tuples of the form (index, string, children).
		:param indent: The amount to indent each level (2 is the min).
		:param width: The maximum line width (set to math.inf for no limit).
		:param space: If True, an extra newline is inserted after each option to space them out more.
		"""
		
		self.options = options
		self.indent = max(indent, 2)
		self.width = width
		self.space = space
		
	def get(self, all = False):
		"""
		
		:param all: Whether to get or not.
		"""
		return "\n".join(self.get_lines(self.options, all = all))
	
	def get_lines(self, options, *, level = 0, all = False, bold = True):
		"""
		
		:param all: Whether to get or not.
		"""
		lines = []
		for index, (option_number, option_body, children) in enumerate(options):
			# First decide if this is the last option at this level or not.
			last_option = (index +1) == len(options)
			
			# Build our string for this option.
			# We start with the option index in square brackets.
			index_string = "[{}] ".format(option_number if option_number is not None else "┐")
			
			# Then we add a space and the option body.
			full_string = "{}{}".format(index_string, option_body)
			
			# Next, we need to wrap if our line is too long.
			# The amount of space available decreases with each indented level.
			# Work out how much space we get.
			#available_space = self.width - (self.indent * level + len(index_string))
			available_space = self.width - (self.indent * level)
			
			# Now wrap and split into lines.
			wrapped_lines = textwrap.wrap(full_string, available_space)
			
			# We actually have slightly less space than this for lines other than the first.
			subsequent_lines = "\n".join(wrapped_lines[1:])
			# Wrap again.
			wrapped_lines = list(chain(wrapped_lines[:1], textwrap.wrap(subsequent_lines, available_space - len(index_string))))
			
			if option_number is None and bold:
				# Make bold.
				# Add control chars.
				if len(wrapped_lines) > 0:
					wrapped_lines[0] = "\033[1m" + wrapped_lines[0]
					wrapped_lines[-1] = wrapped_lines[-1] + "\033[0m"
			
			
			# Add an extra line if we have been asked to.
			if self.space:
				wrapped_lines.append("")
			
			# We now add our indent to each line (which is different for the first line).
			padded_lines = []
			
			for line_number, line in enumerate(wrapped_lines):
				if line_number == 0:
					# If we are the top-most level, then we don't do any padding.
					if level != 0:
						# The first line, add our node 'graphics', which change depending on whether there are options after us.
						indent_string = " ├"  if not last_option else " └"
						# Then pad out to fill width.
						indent_string += "─" * (self.indent -2)
					else:
						indent_string = ""
				else:
					if level != 0:
						# Additional lines; only need to add a continuing down line if there are options after us.
						indent_string = " │"  if not last_option else "  "
						# Then pad out to fill width.
						indent_string += " " * (self.indent -2)
					else:
						indent_string = ""
						
					# If we have children, add another continuation.
					indent_string += " │"  if len(children) > 0 else "  "
					# Add finish padding.
					indent_string += " " * (len(index_string) -2)
					
				# Add indent and add to list.
				padded_lines.append("{}{}".format(indent_string, line))

				
			# If we have children, get them now.
			if len(children) > 0:
				if index == 0 or all or level >= 2:
					child_lines = self.get_lines(children, level = level+1, all = all)
				else:
					# If the all option has not been given, we only print for this first node.
					child_lines = self.get_lines([(None, "(More...)", [])], level = level+1, all = all)
					
				# Pad each appropriately.
				padded_child_lines = []
				for child_line in child_lines:
					if level != 0:
						if not last_option:
							indent_string = " │" + (" " * (self.indent -2))
						else:
							indent_string = " " * self.indent
					else:
						indent_string = ""
					padded_child_lines.append("{}{}".format(indent_string, child_line))
				
			else:
				padded_child_lines = []
						
					
			
			# Done, add the lines.
			lines.extend(padded_lines)
			lines.extend(padded_child_lines)
			
		# All done.
		return lines
				