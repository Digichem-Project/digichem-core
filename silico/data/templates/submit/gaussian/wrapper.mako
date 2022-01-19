## -*- coding: utf-8 -*-
##
<%page args="program" />\
##
#!/bin/bash
##
##
## First, set the env variables gaussian needs.
export GAUSS_SCRDIR="${program.calculation.scratch_directory}"
export ${program.root_environ_name}="${program.root.parent}"
##
## Now we need to load the init script provided by gaussian. Kinda weird but fun way to do init (they even provide a c-shell version for crazy kids (was anyone still using csh in 2016?!)).
. "${program.init_file}" || { >&2 echo "Failed to load Gaussian init file"; exit 1; }
##
##
## Now run Gaussian!
## From the manual it appears that Gaussian does some interpretation of the arguments passed to it, in some cases appending a .gjf suffix (it's not clear when, most likely in instance where the name is lacking a dot (.) character). This should be fine as we always use a .com filename, but it would probably be safer to have Gaussian read from stdin (so it can't get confused).
## We use resolve() here to make the paths absolute (because we are being executed from a different working directory by necessity).
"${program.executable}" "${program.com_file_path.resolve()}" "${program.log_file_path.resolve()}"
exit $?