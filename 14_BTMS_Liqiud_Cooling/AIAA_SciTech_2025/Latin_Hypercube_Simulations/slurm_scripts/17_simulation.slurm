#!/bin/bash
#SBATCH --job-name=17_simulation             # Job name
#SBATCH --output=outputs/output_17.log      # Output log file
#SBATCH --error=outputs/error_17.log        # Error log file
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
python3 -u  BTMS_Battery_Degradation_study.py 9487.312767846126 7744.092631252577 0.5923920361376593 --storage_dir /home/sshekar2/storage/optimization_LHS_11_21/case_17 17

# Deactivate the virtual environment after the script finishes
deactivate
