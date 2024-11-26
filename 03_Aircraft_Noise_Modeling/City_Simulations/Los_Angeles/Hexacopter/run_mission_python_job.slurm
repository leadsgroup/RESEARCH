#!/bin/bash
#SBATCH --job-name=hexacopter_mission_run             # Job name
#SBATCH --output=output_%j.log            # Output log file (%j will be replaced by the job ID)
#SBATCH --error=error_%j.log              # Error log file (%j will be replaced by the job ID)
#SBATCH --nodes=1                         # Request nodes
#SBATCH --ntasks=1                        # Number of tasks (processes)
#SBATCH --cpus-per-task=2                 # Number of CPU cores per task
#SBATCH --partition=primary               # Partition/queue to submit the job
#SBATCH --mail-type=ALL                   # Send email on job start, end and fail
#SBATCH --mail-user=aidanrm2@illinois.edu # Email for notifications

# Load any necessary modules
module load /opt/modulefiles/anaconda3/2024.10  # This has the latest version of python

# Activate the virtual environment
source $HOME/scratch/leads_research/rcaideenv/bin/activate

# Run the Python script
python3 HC_Operations_LA.py

# Deactivate the virtual environment after completion
deactivate
