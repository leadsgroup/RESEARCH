# ----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import os
import subprocess
from pathlib import Path
import time
import latin_hypercube_sampler  # Assuming this contains `generate_lhs_samples`



def main():
    # Define parameters
    total_no_sims = 2  # Total number of simulations to run
    parallel_sims = 2  # Number of simulations to run in parallel

    variable_limits = [(1000, 12500), (500, 9000), (0.1, 0.7)]
    variable_names  = ['Power_HAS', 'Power_HEX', 'Dim_RES']

    storage_dir  = '/home/sshekar2/storage/optimization_LHS_11_21'
    automation_optimizations(total_no_sims, parallel_sims, variable_limits,variable_names,storage_dir)
    
    return

# ----------------------------------------------------------------------
#   Automation Optimization Script
# ---------------------------------------------------------------------
def automation_optimizations(total_simulations, parallel_simulations, variable_limits,variable_names,storage_dir):

    # Generate Latin Hypercube Samples
    lhs_samples = latin_hypercube_sampler.generate_lhs_samples(variable_limits, total_simulations,variable_names)
    print(f"Generated {total_simulations} Latin Hypercube samples.")

    # Create SLURM job scripts directory
    script_dir = Path("slurm_scripts")
    script_dir.mkdir(exist_ok=True)

    # Write SLURM scripts for all simulations
    job_files = []
    for sim_id, params in enumerate(lhs_samples):
        sim_storage_dir = f"{storage_dir}/case_{sim_id}"
        # Create the SLURM job file with the subdirectory path as a string
        job_file = create_slurm_job(sim_id, params, script_dir, sim_storage_dir)
        job_files.append(job_file)

    # Step 4: Submit jobs in parallel batches
    for i in range(0, total_simulations, parallel_simulations):
        batch_jobs = job_files[i:i + parallel_simulations]
        print(f"Submitting batch {i // parallel_simulations + 1}: {len(batch_jobs)} jobs")

        for job_file in batch_jobs:
            submit_slurm_job(job_file)

        # Wait for the current batch to finish
        wait_for_jobs_to_complete()

    print("All simulations completed.")

# ----------------------------------------------------------------------
#   Helper Functions
# ---------------------------------------------------------------------
def create_slurm_job(sim_id, sim_params, script_dir,sim_storage_dir):
    """
    Creates a SLURM job file for a given simulation.
    """
    job_name = f"simulation_{sim_id}"
    job_file = Path(script_dir) / f"{job_name}.slurm"

    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}             # Job name
#SBATCH --output=output_{sim_id}.log      # Output log file
#SBATCH --error=error_{sim_id}.log        # Error log file
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
python3 -u  BTMS_Battery_Degradation_study.py {sim_params[0]} {sim_params[1]} {sim_params[2]} --storage_dir {sim_storage_dir}

# Deactivate the virtual environment after the script finishes
deactivate
"""

    with open(job_file, "w") as f:
        f.write(slurm_script)
    
    return job_file

def submit_slurm_job(job_file):
    """
    Submits a SLURM job file.
    """
    subprocess.run(["sbatch", str(job_file)], check=True)

def wait_for_jobs_to_complete():
    """
    Waits for all SLURM jobs to complete before proceeding.
    """
    print("Waiting for jobs to complete...")
    while True:
        result = subprocess.run(["squeue", "-u", os.getlogin()], capture_output=True, text=True)
        if len(result.stdout.strip().splitlines()) <= 1:  # No running jobs except the header
            break
        time.sleep(30)  # Check every 30 seconds

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()
