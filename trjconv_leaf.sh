#!/bin/bash
root_dir="/run/user/1000/gvfs/smb-share:server=10.43.10.64,share=project/argha/20fs"
ndx_file="${root_dir}/index_analysis.ndx"
dirs=(298 302 306 310 314)
#subdirs=(r1)
subdirs=(r11 r12 r13 r14 r15 r16 r17 r18 r19 r20)
count=0
log_file="${root_dir}/gmx_processing.log"
echo "Starting GROMACS processing at $(date)" > "$log_file"

for dir in "${dirs[@]}"; do
    for subdir in "${subdirs[@]}"; do
        current_dir="${root_dir}/${dir}/${subdir}"
        (
            echo "Processing $dir/$subdir at $(date)" >> "$log_file"
            cp "$ndx_file" "$current_dir"
            cd "$current_dir"
            echo "Running for $dir/$subdir: gmx_mpi convert-tpr for leaf1" >> "$log_file"
            gmx_mpi convert-tpr -s production.tpr -n index_analysis.ndx -o leaf1.tpr <<< "1"
            echo "Running for $dir/$subdir: gmx_mpi convert-tpr for leaf2" >> "$log_file"
            gmx_mpi convert-tpr -s production.tpr -n index_analysis.ndx -o leaf2.tpr <<< "2"
            echo "Running for $dir/$subdir: gmx_mpi trjconv for leaf1" >> "$log_file"
            gmx_mpi trjconv -f production.trr -s production.tpr -n index_analysis.ndx -o leaf1_${dir}_${subdir}_20fs.trr -pbc mol -ur compact -dt 1 <<< "1"
            echo "Running for $dir/$subdir: gmx_mpi trjconv for leaf2" >> "$log_file"
            gmx_mpi trjconv -f production.trr -s production.tpr -n index_analysis.ndx -o leaf2_${dir}_${subdir}_20fs.trr -pbc mol -ur compact -dt 1 <<< "2"
            echo "Finished processing $dir/$subdir at $(date)" >> "$log_file"
        ) &
	let count+=1
        # Wait for a batch of 64 tasks to finish before launching the next batch
        if [ $count -eq 10 ]; then
            wait
            count=0
        fi
    done
done
wait
echo "Completed GROMACS processing at $(date)" >> "$log_file"
