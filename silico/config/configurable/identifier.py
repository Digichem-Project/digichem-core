import shlex

# This splitting mechanism is now deprecated because it interferes with parsing of tag lists.
#
# Splitting with shlex strips out quote marks, normally this is a good thing but it makes it hard to
# specify tag lists which need to be quoted,
# eg:
#
# Single Node SLURM/Gaussian 16/[Opt, B3LYP, Gas-Phase, '6-31G(d,p)']
#
# Becomes a list that looks like:
#
# ["Single Node SLURM", "Gaussian 16", "[Opt, B3LYP, Gas-Phase, 6-31G(d,p)]"]
#
# When this gets parsed by yaml, the final tag path will then get parsed incorrectly into:
#
# ["Single Node SLURM", "Gaussian 16", ["Opt", "B3LYP", "Gas-Phase", "6-31G(d", "p)"]]
#
# The solution is to simply use split("/"), this is less intelligent but preserves whitespace etc.
# To make this work, the '/' character is now disallowed in tag names (which probably makes sense anyway).
 

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
    
