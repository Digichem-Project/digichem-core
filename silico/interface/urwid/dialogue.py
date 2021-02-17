import urwid

class Confirm(urwid.Overlay):
    """
    A confirmation box, designed to be used as a pop-up. Presents the user with a message and two buttons.
    """
    
    def __init__(self, title, message, confirm_action, top, title_attr = "title", message_attr = "dialogue"):
        """
        Constructor for Confirm dialogue boxes.
        
        :param title: Text markup to use as the title.
        :param message: Text markup to use as the message body.
        :param confirm_action: Callback to execute when the user chooses confirm.
        """
        self.tOp = top
        self.confirm_action = confirm_action
        
        # Buttons.
        cancel = urwid.Button("Cancel", self.cancel)
        confirm = urwid.Button("Confirm", self.confirm)
        
        # Get our pile.
        pile = urwid.Pile([
            ('pack', urwid.Divider()),
            ('pack', urwid.Text((title_attr, title), align = "center")),
            ('pack', urwid.Divider()),
            ('pack', urwid.Padding(urwid.Text((message_attr, message), align = "center"), align = "center", width = "pack", left = 1, right = 1)),
            ('pack', urwid.Divider()),
            ('pack', urwid.Columns([
                urwid.AttrMap(urwid.Padding(cancel, align = "center", width = 10), 'bad_button'),
                urwid.AttrMap(urwid.Padding(confirm, align = "center", width = 11), 'good_button')
            ]))
        ])
        
        # Call parent.
        super().__init__(
            #urwid.LineBox(urwid.Filler(pile)),
            urwid.AttrMap(urwid.LineBox(pile), message_attr),
            self.tOp.original_widget,
            align = "center",
            width = ('relative', 80),
            #width = 'pack',
            valign = "middle",
            height = 'pack',
            #height = ('relative', 80)
        )
        
    def cancel(self, *args, **kwargs):
        """
        Function called when the 'cancel' button is chosen.
        """
        self.tOp.back()
        
    def confirm(self, *args, **kwargs):
        """
        Function called when the 'confirm' button is chosen.
        """
        self.confirm_action()
        self.tOp.back()