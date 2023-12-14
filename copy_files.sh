#!/bin/bash

# Parent directories
parent_dirs=("298" "302" "306" "310" "314")

# Subdirectories from r1 to r10
sub_dirs=("r1" "r2" "r3" "r4" "r5" "r6" "r7" "r8" "r9" "r10")

# Remote user and host
#remote_user_host="phyarmi@parampravega.iisc.ac.in"
remote_user_host="phyarmi@beena.physics.iisc.ac.in"
# Loop through the parent directories
for parent in "${parent_dirs[@]}"; do
    # Check if parent directory exists locally, if not, create it
    if [ ! -d "$parent" ]; then
        mkdir -p "$parent"
    fi

    # Loop through the subdirectories
    for sub in "${sub_dirs[@]}"; do
        # Construct the local and remote paths
        local_path="$parent/$sub/leaf*"
        #remote_path="/scratch/phyarmi/project/dppc_20fs/$parent/$sub"
        remote_path="/home/phyarmi/project/analysis/$parent/$sub"

        # Check if subdirectory exists locally, if not, create it
        if [ ! -d "$parent/$sub" ]; then
            mkdir -p "$parent/$sub"
        fi

        # Ensure the directory exists on the remote server
        ssh $remote_user_host "mkdir -p $remote_path"

        # Execute the scp command
        scp $local_path "$remote_user_host:$remote_path"

        # Check if the scp command was successful
        if [ $? -eq 0 ]; then
            echo "Successfully copied files from $local_path to $remote_path"
        else
            echo "Error copying files from $local_path to $remote_path"
        fi
    done
done
