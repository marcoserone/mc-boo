#!/bin/bash

for T in {4..10}; do
	sed s/Param_it/$T/g   mc-bash-I.sh > mc-bash-fixed_$T.sh
	sed s/Param_it/$T/g   fixed_template.txt > fix_external_iteration_$T.txt
 done
