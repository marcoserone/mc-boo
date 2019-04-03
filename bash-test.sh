#!/bin/bash



mkdir fit_external_runs
cd fit_external_runs

start=1
for irun in $(seq $start 5)
do

	file="test_$irun.txt"
	echo " Get["mc-boo.m"]; " > $file
	echo "ParallelTable[ fitExternalWrapper[100,101,123,4,5,50,3000,1001,3/2,11/10,10^(-3),1,temp/2,1 + $irun/10,1/5,"bin_scan",1/10,1/10,0,0] {temp,5,8} ]" >> $file
	echo " Exit[]; " >> $file
	cat $file

done

