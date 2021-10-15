from silico.interface.urwid.misc import Tab_pile
import urwid

class Control_wrapper(Tab_pile):
    """
    Abstract class for classes that wrap widgets with forwards/backwards controls.
    """
    
    def __init__(self, body, top, cancel_callback = None):
        """
        Constructor for Control_wrapper objects.
        
        :param body: The inner widget to wrap.
        :param top: The topmost widget to use for navigation.
        """
        self.body = body
        self.top = top
        self.cancel_callback = self.convert_callback(cancel_callback)
        
        controls = urwid.Columns(self.get_controls())
        
        super().__init__([
            self.body,
            ('pack', controls)
        ])

    def convert_callback(self, callback):
        """
        Convert certain values which can be interpreted as a function into a real function.
        
        The interpretation is as follows:
            - int: A number of steps to call back() on the top widget.
            - otherwise: A functiojn to call.
            
        :returns: The possibly new callback function.
        """
        def back_callback():
            self.top.back(callback)
            # Stop going back any further.
            return False
        
        if isinstance(callback, int):
            return back_callback
        
        else:
            return callback
        
    def keypress(self, size, key):
        """
        Handler for keypress events.
        """
        # Allow children to intercept first.
        key = super().keypress(size, key)
        
        if key == "esc":
            self.cancel()
        
        return key
        
    def cancel(self):
        """
        Function called when the back button is pressed.
        """
        self.retval = False
        if self.cancel_callback is not None:
            if self.cancel_callback() is not False:
                self.top.back()
        
    def get_controls(self):
        raise NotImplementedError("Implement in subclass")

class Confirm_or_cancel(Control_wrapper):
    """
    A widget wrapper that provides a cancel and a confirm button.
    """
    
    def __init__(self, body, top, cancel_callback = None, submit_callback = None):
        """
        Constructor for Confirm_or_cancel objects.
        
        :param body: The inner widget to display.
        :param top: The topmost widget to use for navigation.
        :param cancel_callback: A function that will be called when the cancel button is pressed. If this function returns False, the top-most widget will not be swapped back.
        :param submit_callback: A function that will be called when the submit button is pressed. If this function returns False, the top-most widget will not be swapped back.
        """
        self.submit_callback = self.convert_callback(submit_callback)
        
        # A value indicating which of our two buttons was pressed.
        self.retval = None
        
        super().__init__(body, top, cancel_callback = cancel_callback)

    def get_controls(self):
        """
        Get the buttons used to control this widget.
        """
        return [
            urwid.AttrMap(urwid.Button("Back", lambda button: self.cancel()), "button", "button--focus"),
            urwid.AttrMap(urwid.Button("Confirm", lambda button: self.submit()), "button--good", "button--good--focus")
        ]

    def submit(self):
        """
        Function called when the submit button is pressed.
        """
        self.retval = True
        if self.submit_callback is not None:
            if self.submit_callback() is not False:
                self.top.back()


class Confirm(Control_wrapper):
    """
    A widget wrapper that provides a single button.
    """
    
    def __init__(self, body, top, submit_callback = None):
        super().__init__(body, top, cancel_callback = submit_callback)
    
    def get_controls(self):
        """
        Get the buttons used to control this widget.
        """
        return [
            urwid.AttrMap(urwid.Button("Confirm", lambda button: self.cancel()), "button--good", "button--good--focus")
        ]
    