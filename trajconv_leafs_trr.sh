#!/bin/bash

root_dir="/run/media/argha/disk3/pure_dppc/dt20"

dirs=(298 302 306 310 314)

subdirs=(r1 )

for dir in "${dirs[@]}"; do
    for subdir in "${subdirs[@]}"; do
        
        current_dir="${root_dir}/${dir}/${subdir}"

        cd "$current_dir"

        gmx_mpi convert-tpr -s production.tpr -n index_analysis.ndx -o leaf1.tpr <<< "1"
        
        gmx_mpi convert-tpr -s production.tpr -n index_analysis.ndx -o leaf2.tpr <<< "2"

        gmx_mpi trjconv -f production.trr -s production.tpr -n index_analysis.ndx -o leaf1_${dir}_${subdir}_20fs.trr -pbc mol -ur compact -dt 1 <<< "1"
        
        gmx_mpi trjconv -f production.trr -s production.tpr -n index_analysis.ndx -o leaf2_${dir}_${subdir}_20fs.trr -pbc mol -ur compact -dt 1 <<< "2"

    done
done
