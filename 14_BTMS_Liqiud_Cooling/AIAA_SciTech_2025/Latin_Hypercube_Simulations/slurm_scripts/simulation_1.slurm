#!/bin/bash
#SBATCH --job-name=simulation_1             # Job name
#SBATCH --output=output_1.log      # Output log file
#SBATCH --error=error_1.log        # Error log file
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
python3 -u  BTMS_Battery_Degradation_study.py 8661.078982576946 3997.6385947829126 0.34139015426109237 --storage_dir /home/sshekar2/storage/optimization_LHS_11_21/case_1

# Deactivate the virtual environment after the script finishes
deactivate
