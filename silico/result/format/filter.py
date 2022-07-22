from silico.misc.base import is_number, to_number
from builtins import isinstance
from silico.result.base import Result_object


class Result_filter():
    """
    Class for filtering data form a result set object.
    """
    
    def __init__(self, filter_string):
        """
        """
        # TODO: This splitting is fine for now but is very naive;
        # we offer no support for semi-colons in filter names (escaping etc).
        self.filters = filter_string.split(":")
        
    def filter(self, result):
        """
        Get some data from a result set object.
        """
        # The part of the result set we are currently looking at.
        cur_item = result
        # The data we've filtered.
        data = None
        
        for filter in self.filters:
            if is_number(filter):
                num_filter = to_number(filter)
                
            else:
                num_filter = filter
            
            # If the current item is a list and we have an int, use as an index.
            if isinstance(cur_item, list) and isinstance(num_filter, int):
                child = cur_item[num_filter]
                # If this child object is not a Result_object but the parent is, use the dumped version instead.
                # (because we have to call dump() at some point, and we can't at any point in the child).
                if not isinstance(child, Result_object) and isinstance(cur_item, Result_object):
                    cur_item = cur_item.dump()[num_filter]
                
                else:
                    cur_item = child
                
            # Else, if the filter is the name of an attribute of current item, use that.
            # Likewise to above, if this child item does not have a dump() function, get the dumped equivalent instead.
            elif hasattr(cur_item, filter) and isinstance(getattr(cur_item, filter), Result_object):
                cur_item = getattr(cur_item, filter)
            
            # Else, if our current item has a find method, use that.
            elif hasattr(cur_item, "find"):
                cur_item = cur_item.find(filter)

            # Else, assume filter is a key for a dict.
            # If we are already a dict, just get the item
            elif isinstance(cur_item, dict):
                # If the filter looks like an int, make it a real one.
                # NOTE: This is mostly necessary for supporting emission results,
                # which are stored in a dict
                cur_item = cur_item[num_filter]
            
            # Otherwise, dump the current item to a dict and continue.
            else:
                cur_item = cur_item.dump()[num_filter]
        
        if isinstance(cur_item, Result_object):
            data = cur_item.dump()
        else:
            data = cur_item
        
        return data