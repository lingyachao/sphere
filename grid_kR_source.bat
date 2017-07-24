#!/bin/bash
matlab -r "startup; quit;"

for kR_arg in 1.5 2.5 3.5
do
    for source_arg in 2 2.5 3 3.5
    do
	for D22min_arg in 0.1 0.6 1.1
	do
	    bsub -q medium -c 24:00 -o output.txt matlab -nodisplay -r "fn_grid_kR_source($kR_arg, $source_arg, $D22min_arg)"
	done
    done
done