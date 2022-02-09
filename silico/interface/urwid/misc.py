import urwid.numedit

class Tab_pile(urwid.Pile):
    """
    An enhanced Pile class that allows changing focus between children with the tab and shift-tab keys.
    """
    
    def keypress(self, size, key):
        """
        Handler for keypress events.
        """
        key = super().keypress(size, key)
        
        # Add support for traversing between our children.
        if key == 'tab':
            pass
            # Move down (if we can).
            next_pos = 1
            changed_focus = False
            while self.focus_position +next_pos < len(self.contents):
                # Check if it can be selected.
                next_focus = self.contents[self.focus_position +next_pos][0]
                if next_focus.selectable():
                    self.focus_position = self.focus_position +next_pos
                    changed_focus = True
                    break
                 
                else:
                    next_pos += 1
                     
            if not changed_focus:
                return key
            
        elif key == 'shift tab':
            # Move up (if we can).
            prev_pos = 1
            changed_focus = False
            while self.focus_position -prev_pos >= 0:
                # Check if it can be selected.
                prev_focus = self.contents[self.focus_position -prev_pos][0]
                if prev_focus.selectable():
                    self.focus_position = self.focus_position -prev_pos
                    changed_focus = True
                    break
                
                else:
                    prev_pos += 1
                    
            if not changed_focus:
                return key

        else:
            return key
        
        
class Blank(urwid.Pile):
    """
    A placeholder widget that takes up no space and does nothing.
    """
    
    
    def __init__(self):
        """
        """
        super().__init__([])
        
class IntEditZero(urwid.numedit.IntegerEdit):
    """
    An int edit widget that allows specifying zeroes.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.trimLeadingZeros = False
        
class FloatEditZero(urwid.numedit.FloatEdit):
    """
    An int edit widget that allows specifying zeroes.
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.trimLeadingZeros = False