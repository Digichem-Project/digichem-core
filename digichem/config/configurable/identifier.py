import csv

import silico.log
from silico.misc.base import is_int


class Identifier():
    """
    A class that identifies a configurable. This could be:
     - A single integer (which would reference a loader index).
     - A string (which would reference a loader tag).
     - A list of strings (which would reference a loader path).
    """
    
    def __init__(self, identifier):
        # Use csv to split our string.
        parser = csv.reader([identifier], delimiter = ":")
        
        # First, determine if a Namespace has been defined by a double colon.
        split = next(parser)
        
        # A namespace is defined by an initial double colon ('::').
        # We can check for this by looking if the second item is empty.
        if len(split) > 2 and split[1] == "":
            self.namespace = split[0].strip()
            split = split[2:]
        
        else:
            self.namespace = None
        
        # Remove excess whitespace, add namespace (if given) and remove empty parts.
        self.raw_tag_list = []
        for tag_part in split:
            # Remove surrounding whitespace.
            # NOTE: This is nearly always the right thing to do, but currently there's no way to turn this off.
            # This means we cannot select tags that start or end with whitespace.
            tag_part = tag_part.strip()
            
            # Skip if empty.
            if tag_part == "":
                silico.log.get_logger().warn("ignoring blank section of method tag list")
                continue
            
            # Add
            self.raw_tag_list.append(tag_part)
            
        # Add a list with namespace.
        if self.namespace is not None:
            self.tag_list = [self.namespace + " " + tag_part for tag_part in self.raw_tag_list]
        
        else:
            # NOTE: This list is the same object as raw_tag_list, where as above it is a separate list...
            self.tag_list = self.raw_tag_list
        
        # Now get our actual value.
        # Turn single numbers into int.
        if len(self.tag_list) == 1 and is_int(self.tag_list[0]):
            self.value = int(self.tag_list[0])
        
        else:
            self.value = self.tag_list
    
    
    def __str__(self):
        if not isinstance(self.value, list):
            return str(self.value)
        
        elif self.namespace is not None:
            return self.namespace + ":: " + " : ".join(self.raw_tag_list)
        
        else:
            return " : ".join(self.raw_tag_list)
