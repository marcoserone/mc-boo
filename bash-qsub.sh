#!/bin/bash


force=$1


start=4
for irun in $(seq $start 10)
do

	if [ $irun -eq $start ]; then
		job=`qsub mc-bash-fixed_$irun.sh`
		echo $job
	else
		job_next=`qsub -W depend=afterany:$job mc-bash-fixed_$irun.sh`
		job=$job_next
	fi

done
