#!/bin/bash
#PBS -q extra
#PBS -T extra_sw
#PBS -l nodes=1:ppn=20
#PBS -l walltime=24:00:00
#PBS -mn
#PBS -joe


module load testing
module load mathematica/10.2

math -noprompt -run "<<fix_external_iteration_104.txt" > Output_fix_external_iteration_104.txt

