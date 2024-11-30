#!/bin/bash

# Remote server details
REMOTE_USER="sshekar2"
REMOTE_HOST="malamute.aerospace.illinois.edu"
REMOTE_BASE_PATH="/home/sshekar2/storage/optimization_LHS_11_21"

# Local destination directory
LOCAL_BASE_PATH="/Users/sai/Documents/Research/LEADS_Research/RESEARCH/14_BTMS_Liqiud_Cooling/AIAA_SciTech_2025/Optimization_Results_Analysis"

# Loop through case_0 to case_74
for i in $(seq 0 74); do
  REMOTE_FILE="$REMOTE_BASE_PATH/case_$i/Raw_Data/e_Twin_Otter_nmc_case_${i}.0.xlsx"
  LOCAL_FILE="$LOCAL_BASE_PATH/e_Twin_Otter_nmc_case_${i}.0.xlsx"
  
  # Copy the file using scp
  scp ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_FILE} ${LOCAL_FILE}
  
  # Check if the copy was successful
  if [[ $? -eq 0 ]]; then
    echo "Copied case_$i successfully."
  else
    echo "Failed to copy case_$i." >&2
  fi
done
