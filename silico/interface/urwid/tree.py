import shlex

import urwid
from urwidtrees.widgets import TreeBox
from urwidtrees.decoration import CollapsibleArrowTree

from silico.interface.urwid.misc import Tab_pile
from urwidtrees.tree import SimpleTree

# Default palette for urwid.
palette = [
	('body', 'black', 'light gray'),
	('boldnode', 'black, bold', 'light gray'),
	('focus', 'light gray', 'dark blue', 'standout'),
	('header', 'light gray,bold', 'dark blue', 'standout'),
	('bars', 'dark blue', 'light gray', ''),
	('arrowtip', 'light blue', 'light gray', ''),
	('connectors', 'light red', 'light gray', ''),
	('calcs', 'light magenta, bold', 'light gray'),
	('good_button', 'black', 'dark green')
]

class Tree_node(urwid.WidgetWrap):
	"""
	Selectable text widget used for the individual nodes in the tree.
	"""
	
	def __init__(self, name, *, ID, attr = "body", none_attr = "boldnode", focus_attr = "focus"):
		"""
		Constructor for Tree_node objects.
		
		:param name: The name (the body) of this node.
		:param ID: The unique ID of this node (can be None).
		:param attr: Attribute to use for this node.
		:param node_attr: Alternative attribute that is used for nodes where ID == None.
		:param focus_attr: Attribute that is used when this node has focus.
		"""
		t = urwid.Text(("{}) ".format(ID) if ID is not None else "") + name)	
		w = urwid.AttrMap(t, (attr if ID is None else none_attr), focus_attr)
		super().__init__(w)
		self.ID = ID
		
	@classmethod
	def from_node(self, node, *args, **kwargs):
		"""
		Get a Tree_node widget from a silico.misc.node_printer.Node object.
		
		The children of the returned node will also be set, if appropriate.
		
		:param node: Top level node to convert.
		:param *args: Passed to the constructor of Tree_node.
		:param **args: Passed to the constructor of Tree_node.
		"""
		children = [self.from_node(child) for child in node.children]
		return (self(node.name, *args, ID = node.ID, **kwargs), children if len(children) > 0 else None)

	def selectable(self):
		return True

	def keypress(self, size, key):
		return key

class Calculation_tree(TreeBox):
	"""
	Modified TreeBox.
	"""
	
	def __init__(self, *args, calcbox, **kwargs):
		"""
		:param calcbox: Editable widget that will be populated with selected options. 
		"""
		super().__init__(*args, **kwargs)
		self.calcbox = calcbox
	
	def keypress(self, size, key):
		"""
		Handler for keypress events.
		"""
		if key in ['enter', ' ']:
			# Toggle expanded state.
			# First, get the focused widget.
			w, focuspos = self.get_focus()
			
			# Check we can expand (and selected is not a leaf with no children).
			if not self._tree._tree.is_leaf(focuspos):
				# Determine whether widget is currently expanded.
				if self._tree.is_collapsed(focuspos):
					super().keypress(size, '+')
				else:
					super().keypress(size, '-')
			else:
				# Leaf node, add selected.
				self.calcbox.edit_text += (" " if self.calcbox.edit_text != "" else "") + self.get_calc_string(focuspos) + " "
		elif key in ['backspace']:
			super().keypress(size, '-')
		elif key in ['delete']:
			# Delete the last calculation.
			self.calcbox.calculations = self.calcbox.calculations[:-1]
		else:
			# Not for us to deal with.
			return super().keypress(size, key)
			
	def get_calc_string(self, pos):
		"""
		Get a m/p/c calculation string from a selected child widget.
		"""
		parent = self._tree.parent_position(pos)
		if parent is not None:
			parentstr = self.get_calc_string(parent)
		else:
			parentstr = None
			
		selfID = self._tree[pos].base_widget.ID
		if selfID is not None:
			if parentstr is not None:
				return parentstr + "/" + str(selfID)
			else:
				return str(selfID)
		else:
			if parentstr is not None:
				return parentstr
			else:
				return None
			
class Calcbox(urwid.Edit):
	"""
	Modified edit widget that offers methods for handling calculation strings.
	"""
	
	@property
	def calculations(self):
		"""
		Get a list of the calculation strings currently entered into this calcbox.
		
		Note that modifications to the returned list are not automatically reflected in this edit widget.
		"""
		return shlex.split(self.edit_text)
	
	@calculations.setter
	def calculations(self, value):
		"""
		Set the text of this edit widget from a list of calculation strings.
		"""
		self.edit_text = "  ".join(value)
		
		

class Calculation_browser(Tab_pile):
	"""
	Urwid box widget that allows selecting calculations.
	"""
	
	def __init__(self, nodes, confirm_action, *, confirm_text = "Confirm", calcbox_attr = "body", calcbox_inner_attr = "calcs", confirm_attr = "good_button"):
		"""
		Constructor for Calculation_browser widgets.
		"""
		# Each browser is made up of 3 main child widgets; the browser itself, a text field (calcbox) and the confirm button.
		self.calcbox = Calcbox((calcbox_attr, 'Selected calculations: '))
		self.confirm = urwid.Button(confirm_text, confirm_action)
		self.browser = Calculation_tree(
			CollapsibleArrowTree(
				SimpleTree(nodes),
				is_collapsed = lambda x: True
			),
			calcbox = self.calcbox
		)
		
		# Call our parent constructor to make the pile.
		super().__init__([
			urwid.LineBox(self.browser),
			(4, urwid.LineBox(urwid.AttrMap(urwid.Filler(self.calcbox), calcbox_inner_attr))),
			(1, urwid.AttrMap(urwid.Filler(urwid.Padding(self.confirm, 'center', 11)), confirm_attr))
		])
	
	@classmethod
	def stop(self, *args, **kwargs):
		"""
		This method is called when the user clicks the confirm button.
		"""
		raise urwid.ExitMainLoop() 
		
	@classmethod
	def run(self, node):
		"""
		Interactively run urwid, using a Calculation_browser as the top-most widget.
		
		:param node: The top-most node used to populate the choices in the calculation browser.
		:return: The selected calculations, as a string.
		"""
		browser = Calculation_browser(Tree_node.from_node(node)[1], confirm_action = self.stop)
	
		header = urwid.AttrMap(urwid.Text('Silico Calculation Browser', align = "center"), 'header')
		footer = urwid.AttrMap(urwid.Text('ENTER: select   DELETE: delete   E: expand all   C: contract all   ctrl-c: quit', align = "center"), 'focus')
		#enclose in a frame
		urwid.MainLoop(
			urwid.Frame(
				urwid.AttrMap(browser, 'body'),
				footer=footer,
				header=header
			),
			palette,
	
		).run()  # go
	
		return browser.calcbox.edit_text
	
