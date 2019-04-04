#!/bin/bash




start=4
for irun in $(seq $start 6)
do

	file="split_$irun.txt"
	echo " Get[\"mc-boo.m\"]; " > $file
	echo "fixedExternalWrapperSplit[$irun,100,101,123,4,6,5,200,101,1,11/10,10^(-3),1,1,\"bblb\",1/10,1/10,0,0]" >> $file
	echo " Exit[]; " >> $file
	cat $file

done

for irun in $(seq $start 6)
do

	file="split_$irun.txt"
	/opt/sissa/sissa-mathematica113/root/bin/math -noprompt -run "<<$file" > Output_split_$irun.txt

done
