#!/bin/bash
#
#
# Helper script for tab completion.
#
# TODO: Improve by someone who knows better bash than I...

_silico_complete() {
	
	if [ "$COMP_CWORD" -eq 1 ]; then
		# First argument is a sub command name.
		COMPREPLY=($(compgen -W "config convert interactive report result status submit -v -h" "${COMP_WORDS[1]}"))
	else
	    local IFS=$'\n'
	    files=($(compgen -f "${COMP_WORDS[$COMP_CWORD]}"))
		COMPREPLY=("${files[@]}")
	fi
}

complete -o filenames -F _silico_complete silico
