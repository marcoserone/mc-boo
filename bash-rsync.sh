#!/bin/bash
rsync -avuz $prune-empty-dirs $include='Res*.txt' $include='Log*' $include='info_*' $include='*/' $exclude='*' uluviano@frontend1.hpc.sissa.it:~/uluviano datsync
