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
## Now run the commands we've been given.
${program.define_executable} < ${shlex.quote(str(program.define_input_path.resolve()))} > ${shlex.quote(str(program.define_output_path.resolve()))} || { >&2 echo "Failed to execute define"; exit 1; }
##
## Done
exit $?