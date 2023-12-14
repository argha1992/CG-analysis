#!/bin/bash

SCRIPT_BASE_DIR="/run/user/1000/gvfs/smb-share:server=10.43.10.64,share=project/argha/20fs/"

# Define the path to your Python script
SCRIPT_PATH="${SCRIPT_BASE_DIR}/msdtime_avg_grid_mpi.py"

# Arrays of directories and subdirectories
declare -a TEMP_DIRS=("298" "302" "306" "310" "314")
declare -a RUN_DIRS=("r11" "r12" "r13" "r14" "r15" "r16" "r17" "r18" "r19" "r20")

# Array of .tpr files
declare -a TPR_FILES=("leaf1.tpr" "leaf2.tpr")

# Starting directory
START_DIR=$(pwd)

# Loop through each temperature directory
for TEMP_DIR in "${TEMP_DIRS[@]}"; do

    # Loop through each run directory within the current temperature directory
    for RUN_DIR in "${RUN_DIRS[@]}"; do

        # Loop through the .tpr files, generate the corresponding .trr filename, and run the Python script
        for TPR_FILE in "${TPR_FILES[@]}"; do
            # Extract base name from .tpr (e.g., "leaf1" from "leaf1.tpr")
            BASE_NAME=$(basename $TPR_FILE .tpr)
            
            # Construct the corresponding .trr filename
            TRR_FILE="${BASE_NAME}_${TEMP_DIR}_${RUN_DIR}_20fs.trr"
            
            # Absolute paths for the TPR and TRR files
            ABSOLUTE_TPR_PATH="${START_DIR}/${TEMP_DIR}/${RUN_DIR}/${TPR_FILE}"
            ABSOLUTE_TRR_PATH="${START_DIR}/${TEMP_DIR}/${RUN_DIR}/${TRR_FILE}"
            
            mpirun python3 $SCRIPT_PATH -t $ABSOLUTE_TPR_PATH -r $ABSOLUTE_TRR_PATH >> ${START_DIR}/script_execution_2.log 2>&1

        done

    done

done
