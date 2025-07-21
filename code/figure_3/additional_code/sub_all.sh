#!/bin/bash

# List of jobs to submit (replace these with actual job scripts or commands)
jobs=(
	## Commnands from  Download all here
)

# Function to check if the number of running or queued jobs is less than or equal to 25
check_jobs() {
    job_count=$(squeue -u $USER -t R,PD,ST,CD --noheader | wc -l)
    echo $job_count
}

# Loop through each job in the list
for job in "${jobs[@]}"
do
    # Wait until the number of jobs is less than or equal to 25
    while [ $(check_jobs) -gt 22 ]; do
        echo "Job count is more than 22. Waiting..."
        sleep 30  # Wait for 30 seconds before checking again
    done
    
    echo "Submitting job: $job..."
    eval $job  # Submit the job from the list
done

