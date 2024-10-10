{\rtf1\ansi\ansicpg1252\cocoartf2818
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 AndaleMono;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red255\green255\blue255;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\cssrgb\c100000\c100000\c100000\c0;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs24 \cf2 \cb3 \CocoaLigature0 #!/bin/bash\
#SBATCH --job-name= Aidan_Test_Run_1 # Job name\
#SBATCH --output=output_%j.log            # Output log file (%j will be replaced by the job ID)\
#SBATCH --error=error_%j.log              # Error log file (%j will be replaced by the job ID)\
#SBATCH --ntasks=1                        # Number of tasks (processes)\
#SBATCH --cpus-per-task=4                # Number of CPU cores per task\
#SBATCH --mem=16G                          # Memory per node\
#SBATCH --time=01:00:00                   # Time limit hrs:min:sec\
#SBATCH --partition=secondary             # Partition/queue to submit the job\
#SBATCH --mail-type=ALL                   # Send email on job start, end and fail\
#SBATCH --mail-user=aidanrm2@illinois.edu # Email for notifications\
\
# Load any necessary modules\
module load anaconda3/2024.06-Jun # This has the latest version of python\
\
# Activate the virtual environment\
source $HOME/scratch/RCAIDE_Develop_Virtual_Environment/bin/activate\
\
# Run the Python script\
python3 Curved_Flight_Test.py\
\
# Deactivate the virtual environment after completion\
deactivate\
~           }