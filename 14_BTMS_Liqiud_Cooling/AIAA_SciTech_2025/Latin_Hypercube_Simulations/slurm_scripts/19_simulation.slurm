#!/bin/bash
#SBATCH --job-name=19_simulation             # Job name
#SBATCH --output=outputs/output_19.log      # Output log file
#SBATCH --error=outputs/error_19.log        # Error log file
#SBATCH --nodes=1                         # Request one node
#SBATCH --ntasks=1                        # Number of tasks (processes)
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --partition=primary               # Partition/queue to submit the job

# Load the necessary module
module load /opt/modulefiles/anaconda3/2024.10

# Activate the virtual environment
source /home/sshekar2/scratch/leads_research/.rcaideenv/bin/activate

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURMD_NODENAME"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Number of CPUs: $SLURM_CPUS_ON_NODE"

# Run the simulation with parameters
python3 -u  BTMS_Battery_Degradation_study.py 1907.5886407088653 2543.5701738359317 0.5605806921146017 --storage_dir /home/sshekar2/storage/optimization_LHS_11_21/case_19 19

# Deactivate the virtual environment after the script finishes
deactivate
