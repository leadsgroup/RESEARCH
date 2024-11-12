import os
import pandas as pd
import numpy as np
from RCAIDE.load import load as load_results
from RCAIDE.save import save as save_results 


def initialize_raw_data(csv_files, raw_data={}):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    for file_path in csv_files:
        base_name = os.path.basename(file_path)
        c_rate_str = base_name.replace('Battery_Simulation_C_rate_', '').replace('.csv', '').replace('_', '.')
        full_path = os.path.join(current_dir, file_path)
        df = pd.read_csv(full_path)
        
        # Initialize dictionary for the current C-rate if it doesn't exist
        raw_data[c_rate_str] = {}
        
        # Loop through each unique initial temperature and store the data
        for temp in df['Initial_Temperature'].unique():
            temp_df = df[df['Initial_Temperature'] == temp]
            
            # Prepare data lists for discharge, voltage, and temperature
            discharge = temp_df['Capacity'].tolist()
            voltage = temp_df['Voltage'].tolist()
            temperature = temp_df['Mean_Temperature'].tolist()
            
            # Store the data under the C-rate and initial temperature key
            raw_data[c_rate_str][f'Initial Temperature {temp}'] = {
                'discharge': discharge,
                'voltage': voltage,
                'temperature': temperature
            }
    return raw_data

def main():
  
    # List of CSV files
    csv_files = [
    'Battery_Simulation_C_rate_0_1.csv',
    'Battery_Simulation_C_rate_0_5.csv',
    'Battery_Simulation_C_rate_1.csv',
    'Battery_Simulation_C_rate_2.csv',
    'Battery_Simulation_C_rate_4.csv',
    'Battery_Simulation_C_rate_6.csv',
    'Battery_Simulation_C_rate_8.csv',
    'Battery_Simulation_C_rate_10.csv',
    'Battery_Simulation_C_rate_12.csv',
    'Battery_Simulation_C_rate_14.csv',
    'Battery_Simulation_C_rate_16.csv',
    'Battery_Simulation_C_rate_18.csv',
    'Battery_Simulation_C_rate_20.csv'
]
    
    # Initialize the raw_data dictionary
    raw_data = initialize_raw_data(csv_files)
    save_results(raw_data, os.path.join(os.path.dirname(__file__), 'lfp_raw_data.res'))



if __name__ == "__main__":
    main()
