#!/bin/bash
#
#
# Helper script for tab completion.
#

_silico_complete() {
	#echo "${#COMP_WORDS[@]}"
	
	if [ "$COMP_CWORD" -eq 1 ]; then
		# First argument is a sub command name.
		COMPREPLY=($(compgen -W "submit convert report result status" "${COMP_WORDS[1]}"))
	else
		COMPREPLY=($(compgen -f "${COMP_WORDS[$COMP_CWORD]}"))
	fi
}

complete -F _silico_complete silico

# _silico_dev_complete() {
# 	
# }
# 
# complete -F _silico_complete silico-dev