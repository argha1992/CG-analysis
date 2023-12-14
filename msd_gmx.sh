#!/bin/bash
root_dir="/run/user/1000/gvfs/smb-share:server=nas-nano.local,share=project/argha/mixed_m2"
dirs=(298 302 306 310 314)
subdirs=(r1 r2 r3 r4 r5)
log_file="${root_dir}/gmx_processing_msd.log"
msd_dir="${root_dir}/msd"
mkdir -p "$msd_dir"

echo "Starting GROMACS processing at $(date)" > "$log_file"

for dir in "${dirs[@]}"; do
    for subdir in "${subdirs[@]}"; do
        current_dir="${root_dir}/${dir}/${subdir}"
        (
            echo "Processing $dir/$subdir at $(date)" >> "$log_file"
            cd "$current_dir"

            for leaf in 1 2; do
                for time in 3 4; do
                    start_time=$((time * 1000000))
                    end_time=$(((time + 1) * 1000000))
                    time_interval="${time}$((${time} + 1))"

                    echo "Running for $dir/$subdir: gmx_mpi msd for leaf${leaf} of ${time}-${time+1}us" >> "$log_file"
                    gmx_mpi msd -f leaf${leaf}_${dir}_${subdir}_dt100.xtc -s leaf${leaf}.tpr -lateral z -dt 1000 -o "${msd_dir}/msd_leaf${leaf}_${dir}_${subdir}_${time_interval}.xvg" -b $start_time -e $end_time <<< "2"
                done
            done

            echo "Finished processing $dir/$subdir at $(date)" >> "$log_file"
        ) &
        if (( ++count == 1 )); then
            wait
            count=0
        fi
    done
done
wait
echo "Completed GROMACS MSD processing at $(date)" >> "$log_file"

