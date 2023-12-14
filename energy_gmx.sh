#!/bin/bash
root_dir="/run/user/1000/gvfs/smb-share:server=nas-nano.local,share=project/argha/dapc"
#dirs=(298)
dirs=(298 302 306 310 314)
#subdirs=(replica1)
subdirs=(r1 r2 r3 r4 r5 r6 r7 r8 r9 r10)
log_file="${root_dir}/gmx_processing_energy.log"
energy_dir="${root_dir}/energy_mixed"
mkdir -p "$energy_dir"

echo "Starting GROMACS processing at $(date)" > "$log_file"


for dir in "${dirs[@]}"; do
    for subdir in "${subdirs[@]}"; do
        current_dir="${root_dir}/${dir}/${subdir}"
        (
            echo "Processing $dir/$subdir at $(date)" >> "$log_file"
            cd "$current_dir"
            echo "Running for $dir/$subdir: gmx_mpi energy" >> "$log_file"
            gmx_mpi energy -f production.edr -s production.tpr -o "${energy_dir}/energy_mixed_${dir}_${subdir}.xvg" <<< "10"
            echo "Finished processing $dir/$subdir at $(date)" >> "$log_file"
        ) &
        if (( ++count == 1 )); then
            wait
            count=0
        fi
    done
done
wait
echo "Completed GROMACS ENERGY processing at $(date)" >> "$log_file"

