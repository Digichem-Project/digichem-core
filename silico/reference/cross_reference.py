

class Captions():
    """
    A class for setting and retrieving cross-references.
    """
    
    def __init__(self):
        self.database = {}
                
    def cross(self, caption_type, name):
        """
        Retrieve a cross-reference to a caption.
        
        If the caption has not previously been referenced, it will first be created.
        """
        try:
            return self.database[caption_type].index(name) +1
        except (KeyError, ValueError):
            self.register(caption_type, name)
            return self.database[caption_type].index(name) +1
        
    def register(self, caption_type, name):
        """
        Insert a new caption into the database.
        """
        if caption_type not in self.database:
            self.database[caption_type] = []
        
        # If this caption has already been registered, issue a warning.
        if name in self.database[caption_type]:
            raise ValueError("The caption '{}'/'{}' has already been registered ({})".format(caption_type, name, self.database[caption_type].index(name)))
        
        else:
            self.database[caption_type].append(name)
            
    def __call__(self, caption_type, name):
        return self.cross(caption_type, name)