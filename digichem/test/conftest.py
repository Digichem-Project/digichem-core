import numpy

# Set numpy errors (not sure why this isn't the default...)
numpy.seterr(invalid = 'raise', divide = 'raise')
