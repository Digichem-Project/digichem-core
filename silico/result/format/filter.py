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
        # The header/name of the current item.
        header = self.filters
        
        for filter in self.filters:
            if is_number(filter):
                num_filter = to_number(filter)
                
            else:
                num_filter = filter
            
            # If the current item is a list (or has a 'values' item) and we have an int, use as an index.
            if (isinstance(cur_item, list) or isinstance(cur_item, dict) and "values" in cur_item) and isinstance(num_filter, int):
                try:
                    child = cur_item[num_filter] if isinstance(cur_item, list) else cur_item["values"][num_filter]
                    # If this child object is not a Result_object but the parent is, use the dumped version instead.
                    # (because we have to call dump() at some point, and we can't at any point in the child).
                    if not isinstance(child, Result_object) and isinstance(cur_item, Result_object):
                        cur_item = cur_item.dump()[num_filter]
                    
                    else:
                        cur_item = child
                
                except IndexError:
                    raise IndexError("filter '{}' is out of range".format(num_filter)) from None
                
            # Else, if the filter is the name of an attribute of current item, use that.
            # Likewise to above, if this child item does not have a dump() function, get the dumped equivalent instead.
            elif hasattr(cur_item, filter) and isinstance(getattr(cur_item, filter), Result_object):
                cur_item = getattr(cur_item, filter)
            
            # Else, if the filter is the name of a key in either cur_item or cur_item.dump(), use that.
            elif isinstance(cur_item, dict):
                # If the filter looks like an int, make it a real one.
                # NOTE: This is mostly necessary for supporting emission results,
                # which are stored in a dict
                cur_item = cur_item[num_filter]
                
            elif num_filter in cur_item.dump():
                cur_item = cur_item.dump()[num_filter]
            
            # Else, if our current item has a find method, use that.
            elif hasattr(cur_item, "find"):
                cur_item = cur_item.find(filter)
                
            else:
                # Panic, we don't know what to filter by.
                raise ValueError("Unable to filter by '{}'".format(filter))
        
        if isinstance(cur_item, Result_object):
            cur_item = cur_item.dump()
            
        # Add header.
        if len(header) > 0:
            data = {":".join(header): cur_item}
        
        else:
            data = cur_item
        
        return data