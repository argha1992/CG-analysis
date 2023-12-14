#!/bin/bash

SCRIPT_PATH="/run/media/argha/disk3/pure_dppc/dt20/correlated_flow1.py"
DIR="/run/media/argha/disk3/pure_dppc/dt20"
LOG_FILE="/run/media/argha/disk3/pure_dppc/dt20/execution.log"

# Counter for concurrent processes
counter=0

# Initialize the log file (overwrites if it exists)
echo "Starting the correlation flow processing:" > $LOG_FILE

for temp in 298 302 306 310 314; do
    for replica in r1 r2 r3 r4 r5 r6 r7 r8 r9 r10; do
        for leaf in leaf1 leaf2; do
	    tpr_file="${DIR}/${temp}/${replica}/${leaf}.tpr"
	    trr_file="${DIR}/${temp}/${replica}/${leaf}_${temp}_${replica}_20fs.trr"

            if [[ -e "$tpr_file" ]] && [[ -e "$trr_file" ]]; then
                echo "Processing: $tpr_file $trr_file" >> $LOG_FILE
                mpirun python3 $SCRIPT_PATH $tpr_file $trr_file &
                let counter+=1
                # If 8 processes are running, wait
                if [[ $counter -eq 1 ]]; then
                    wait
                    counter=0
                fi
            fi
        done
    done
done

# Wait for any remaining processes to finish
wait

echo "Finished all processes." >> $LOG_FILE
