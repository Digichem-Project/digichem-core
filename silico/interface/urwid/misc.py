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
            # Move down (if we can).
            if self.focus_position +1 < len(self.contents):
                self.focus_position = self.focus_position +1
            
            else:
                # Continue.
                return key
            
        elif key == 'shift tab':
            # Move up (if we can).
            if self.focus_position > 0:
                self.focus_position = self.focus_position -1
            
            else :
                # Continue;
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