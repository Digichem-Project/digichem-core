from silico.interface.urwid.setedit.configurable import Configurable_browser


class Swappable_browser(Configurable_browser):
    """
    A settings browser designed for Swappable widgets.
    """
    
    def __init__(self, view):
        super().__init__(view.top, view, on_change_callback = view.on_settings_change)


