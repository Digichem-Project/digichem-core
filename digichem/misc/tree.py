# Copyright (C) 2013  Patrick Totzke <patricktotzke@gmail.com>
# This file is released under the GNU GPL, version 3 or a later revision.

# Modified from the above source

import urwid
from urwidtrees.tree import SimpleTree
from urwidtrees.widgets import TreeBox
from urwidtrees.decoration import ArrowTree  # for Decoration
from urwidtrees.decoration import CollapsibleArrowTree
import shlex


# define some colours
palette = [
	('body', 'black', 'light gray'),
	('boldnode', 'black, bold', 'light gray'),
	('focus', 'light gray', 'dark blue', 'standout'),
	('header', 'light gray,bold', 'dark blue', 'standout'),
	('bars', 'dark blue', 'light gray', ''),
	('arrowtip', 'light blue', 'light gray', ''),
	('connectors', 'light red', 'light gray', ''),
	('calcs', 'light magenta, bold', 'light gray'),
	('button', 'black', 'dark green')
]

# We use selectable Text widgets for our example..

class LeeBox(TreeBox):
	"""
	Tmp custom TreeBox implementation that supports some additional key presses.
	"""
	
	def __init__(self, *args, calcbox, **kwargs):
		"""
		:param calcbox: Editable widget that will be populated with selected options. 
		"""
		super().__init__(*args, **kwargs)
		self.calcbox = calcbox
	
	def keypress(self, size, key):
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
				# Lead node, add selected.
				self.calcbox.edit_text += (" " if self.calcbox.edit_text != "" else "") + self.get_calc_string(focuspos) + " "
		elif key in ['backspace']:
			super().keypress(size, '-')
		elif key in ['delete']:
			# Delete the last calculation.
			# Split on each calc string.
			calc_strings = shlex.split(self.calcbox.edit_text)
			# Remove the last (if we can).
			calc_strings = calc_strings[:-1]
			# Rebuild.
			self.calcbox.edit_text = "  ".join(calc_strings)
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

class FocusableText(urwid.WidgetWrap):
	"""Selectable Text used for nodes in our example"""
	def __init__(self, name, *, ID):
		t = urwid.Text(("{}) ".format(ID) if ID is not None else "") + name)	
		w = urwid.AttrMap(t, ('body' if ID is None else 'boldnode'), 'focus')
		urwid.WidgetWrap.__init__(self, w)
		self.ID = ID

	def selectable(self):
		return True

	def keypress(self, size, key):
		return key

class Tabpile(urwid.Pile):
	
	def keypress(self, size, key):
		if key == 'tab':
			#print(self.focus_position)
			self.focus_position = self.focus_position +1 if self.focus_position+1 < len(self.contents) else self.focus_position
		elif key == 'shift tab':
			self.focus_position = self.focus_position -1 if self.focus_position > 0 else self.focus_position
		else:
			return super().keypress(size, key)

def node_list_to_tree(node):
	"""
	"""
	children = [node_list_to_tree(child) for child in node.children]
	return (FocusableText(node.name, ID = node.ID), children if len(children) > 0 else None)

def construct_example_tree(node):
	# define a list of tree structures to be passed on to SimpleTree
	forrest = [node_list_to_tree(node)]

	# stick out test tree into a SimpleTree and return
	return SimpleTree(forrest)

def stop(*args, **kwargs):
	raise urwid.ExitMainLoop()
	
def run(node):
	# get example tree
	stree = SimpleTree([node_list_to_tree(node)])
	stree = SimpleTree(node_list_to_tree(node)[1])
	# Here, we add some decoration by wrapping the tree using ArrowTree.
	atree = CollapsibleArrowTree(stree,
						is_collapsed = lambda x: stree.depth(x) > -1,
					  # customize at will..
					  # arrow_hbar_char=u'\u2550',
					  # arrow_vbar_char=u'\u2551',
					  #  arrow_tip_char=u'\u25B7 ',
					  #  indent=6
					  # arrow_connector_tchar=u'\u2560',
					  # arrow_connector_lchar=u'\u255A',
					  )

	# put the into a treebox
	#treebox = urwid.AttrMap(TreeBox(atree), 'body')
	#treebox = TreeBox(atree)
	
	# An editable textbox where selected calcs will be placed.
	calcbox = urwid.Edit(('body', 'Selected calculations: '))
	
	treebox = LeeBox(atree, calcbox = calcbox)
	
	# Container for top-level widgets.
	#root = urwid.Pile([treebox, (2, urwid.Filler(calcbox))])
	#root = urwid.Pile([
	root = Tabpile([
		urwid.LineBox(treebox),
		(4, urwid.LineBox(urwid.AttrMap(urwid.Filler(calcbox), 'calcs'))),
		(1, urwid.AttrMap(urwid.Filler(urwid.Padding(urwid.Button("Submit", stop), 'center', 10)), 'button'))
	])
	#root = urwid.Pile([treebox, (4, urwid.LineBox(urwid.Filler(calcbox)))])
	
	
	#rootwidget = urwid.AttrMap(treebox, 'body')
	#add a text footer
	header = urwid.AttrMap(urwid.Text('Silico Calculation Browser', align = "center"), 'header')
	footer = urwid.AttrMap(urwid.Text('ENTER: select   DELETE: delete   E: expand all   C: contract all   ctrl-c: quit', align = "center"), 'focus')
	#enclose in a frame
	urwid.MainLoop(
		urwid.Frame(
			urwid.AttrMap(root, 'body'),
			footer=footer,
			header=header
		),
		palette,

	).run()  # go

	return calcbox.edit_text
	
	
	
	