#!/bin/bash
#SBATCH --job-name=1_simulation             # Job name
#SBATCH --output=outputs/output_1.log      # Output log file
#SBATCH --error=outputs/error_1.log        # Error log file
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
python3 -u  BTMS_Battery_Degradation_study.py 1325.2529971975935 4338.670323923622 0.5143522581019287 --storage_dir /home/sshekar2/storage/optimization_LHS_11_21/case_1 1

# Deactivate the virtual environment after the script finishes
deactivate
