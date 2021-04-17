#!/bin/bash

for dir in */; do
    cd "$dir"
	sbatch jobscript.sh
	cd ..
done
squeue -u dvolkov
