import collections

def tail(file, num_lines = 20):
    """
    Return the last n lines of a file.
    """
    last_lines = collections.deque(maxlen = num_lines)
    
    # We'll read through the file from the top.
    # This is probably inefficient for huge files but is easy to implement and I don't care just now.
    for line in file:
        last_lines.append(line)
        
    return list(last_lines)