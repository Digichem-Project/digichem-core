import shlex


class ID_splitter(shlex.shlex):
    """
    A custom parser used for splitting strings which identify certain configurables.
    
    An identifier is something that can be used to identify a particular configurable loader.
    
    Typically they are either an integer index or a string describing a path to a particular loader.
    """
    
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, posix = True)
        # We only have one splitting character.
        self.whitespace = "/"
        self.whitespace_split = True
    
