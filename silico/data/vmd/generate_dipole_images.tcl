# This script expects 9 arguments, which are:
# cube_file    A .cube file to load molecular data from.
set cube_file [lindex $argv 0]
# tcl_common    Path to our common library to source.
set tcl_common [lindex $argv 1]
# A string describing the rendering style to use; 'silico' or 'gaussian'.
set rendering_style [lindex $argv 2]
# A string of x,y,z translations to perform.
set translations [lindex $argv 3]
# rotations        A string list of rotations to perform.
set rotations [lindex $argv 4]
# dipole_start    {x y z} coordinates of the start of the dipole.
set dipole_start [lindex $argv 5]
# dipole_end    {x y z} coordinates of the end of the dipole.
set dipole_end [lindex $argv 6]
# x0y0z0        A file name to write one of the output images to (in png format).
set x0y0z0 [lindex $argv 7]
# x90y0z0        A file name to write one of the output images to (in png format).
set x90y0z0 [lindex $argv 8]
# x0y90z0        A file name to write one of the output images to (in png format).
set x0y90z0 [lindex $argv 9]
# x45y45z45        A file name to write one of the output images to (in png format).
set x45y45z45 [lindex $argv 10]

# Load our common library.
source $tcl_common

# Set our visual style
use_style $rendering_style

# Load our molecule and keep track of its numerical handle.
set mol_handle [molecule new $cube_file]

# Rotate as we've been told.
#rotate_molecule $mol_handle $translations $rotations

# Make our molecule transparent so we can better see our dipole arrow.
material change Opacity Transparent 0.5
molecule modmaterial 0 $mol_handle Transparent

# Use standard display settings.
standard_molecule_style 0 $mol_handle

# Draw our dipole moment.
draw color red
draw arrow $dipole_start $dipole_end

# And save our pictures.
render_images $rotations $x0y0z0 $x90y0z0 $x0y90z0 $x45y45z45 

# All done.
exit
