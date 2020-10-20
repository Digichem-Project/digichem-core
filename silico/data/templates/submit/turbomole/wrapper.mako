## -*- coding: utf-8 -*-
##
<%!
	import shlex
%>\
##
<%page args="program" />\
##
##
## First, set the env variables turbomole needs.
export TURBODIR=${program.root}
##
## Now load init script.
. ${shlex.quote(str(program.init_file))} || { >&2 echo "Failed to load Turbomole init file"; exit 1; }
##
##