
class Invalid_method_parent(Exception):
    """
    A calculation was created with a program that can not run it, or a program was created with a destination that cannot run it.
    """
    
    def __init__(self, method_part, parent):
        """
        """
        self.method_part = method_part
        self.parent = parent
        
        super().__init__("The {} '{}' is not compatible with {} '{}'".format(method_part.TYPE, method_part.name, parent.TYPE, parent.name))