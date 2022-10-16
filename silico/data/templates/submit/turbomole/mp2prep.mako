## -*- coding: utf-8 -*-
<%!
    import shlex
%>\
##
<%page args="program"/>\
##
##
#!/bin/bash
##
##
## Import turbomole init script.
<%include file="wrapper.mako"/>\
##
## Now run mp2prep
## We specify up to two options:
## - -e or -g for energy or gradient respectively.
## - Scratch dir if we're using one.
${program.mp2prep_executable} ${"-g" if program.calculation.properties['opt']['calc'] or program.calculation.properties['freq']['calc'] else "-e"} ${shlex.quote(str(program.scratch_base.resolve())) if program.scratch_base is not None else ""} || { >&2 echo "Failed to execute mp2prep"; exit 1; }
##
## Done
exit $?