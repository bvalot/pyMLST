#!/usr/bin/env bash

_wgmlst_completions()
{
	if [ "${#COMP_WORDS[@]}" != "2" ]; then
		return
	fi

	COMPREPLY=($(compgen -W "add_strain create_db" "${COMP_WORDS[1]}"))
}

complete -F _wgmlst_completions wgMLST
