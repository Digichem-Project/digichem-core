
def andjoin(listlike):
    """
    Joins a list of string-like variables with commans, except for the last item which is joined with 'and'.
    
    For example:
    >>> andjoin(["item1"])
    "item1"
    >>> andjoin(["item1", "item2"])
    "item1 and item2"
    >>> andjoin(["item1", "item2", "item3"])
    "item1, item2 and item3"
    """
    # Now combine.
    joined = ", ".join([str(item) for item in listlike[:-1]])
    
    if len(listlike) > 1:
        joined += " and " + str(listlike[-1])
    elif len(listlike) == 1:
        joined = listlike[0]
    
    return joined