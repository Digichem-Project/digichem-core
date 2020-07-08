import textwrap
from itertools import chain

class Node():
	"""
	Object containing data that can be printed with Node_printer.
	"""
	
	def __init__(self, ID, name, children = None):
		"""
		Constructor for Node objects.
		
		:param ID: The unique ID of the Node, None is used to indicate an unselectable Node (for grouping purposes).
		:param name: The name of the Node.
		:param children: Optional list of children.
		"""
		self.ID = ID
		self.name = name
		self.children = children if children is not None else []
	
	@classmethod
	def from_list(self, configurable_lists, top_name = "Calculations"):
		"""
		"""
		# First, get a top-level node.
		base_node = self(None, top_name)
		
		# Add.
		for config in configurable_lists[0]:
			self.from_configurable(base_node, config, configurable_lists)
	
	@classmethod
	def from_configurable(self, base_node, config, configs):
		"""
		"""
		# First, turn our config into a Node.
		for group_name in config.GROUP:
			if len(base_node.children) == 0 or base_node.children[-1].name != group_name:
				# Either the base node has no children, or the last child is not our group.
				# Add a new group.
				base_node.children.append(self(None, group_name))
				
			# Set this new node as our base.
			base_node = base_node.children[-1]
			
		# Now add the config to the (current) base node.
		node = self(config.ID, config.NAME if config.GROUP_NAME is not None else config.GROUP_NAME)
		
		# Add our children, using the new node as base.
		if len(configs) > 0:
			for child in config.get_children(configs[1]):
				self.from_configurable(node, child, configs[1:])
			
		# Done, return Node of convenience (although it probably isn't very useful).
		return node
		

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
				