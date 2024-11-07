import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


def main():
    # Read the data from the CSV file
    file_paths = ['HPPC_Panasonic_NCR18650GA_MZTC_Channel_5_Wb_1.CSV',
                 'HPPC_Panasonic_NCR18650GA_MZTC_Channel_6_Wb_1.CSV',
                 'HPPC_Panasonic_NCR18650GA_MZTC_Channel_7_Wb_1.CSV',
                 'HPPC_Panasonic_NCR18650GA_MZTC_Channel_8_Wb_1.CSV']
    
    # Create full paths for each file
    file_paths = [os.path.join(os.path.dirname(__file__), f) for f in file_paths]
    df = read_data(file_paths)
    
    plot_data(df)
    return

def read_data(file_paths):
    df_list = []
    for file_path in file_paths:
        df = pd.read_csv(file_path)
        # Get the last Discharge Capacity value from Step 3
        step3_last_capacity = df[df['Step Index'] == 3]['Discharge Capacity (Ah)'].iloc[-1]
        # Filter data where Step Index > 7
        filtered_df = df[df['Step Index'] > 7].copy()
        # Reset Discharge Capacity relative to the last value from Step 3
        filtered_df['Discharge Capacity (Ah)'] = filtered_df['Discharge Capacity (Ah)'] - step3_last_capacity
        df_list.append(filtered_df)
    return df_list

def plot_data(df_list):
    num_plots = len(df_list)
    variables = [
        ('Voltage (V)', 'Voltage vs Time'),
        ('Current (A)', 'Current vs Time'),
        ('Aux_Temperature_2 (C)', 'Temperature vs Time'),
        ('Capacity (Ah)', 'Capacity vs Time'),
        ('Discharge Capacity (Ah)', 'Discharge Capacity vs Time')
    ]
    
    for y_var, title_suffix in variables:
        fig, axes = plt.subplots(int(num_plots/2), 2, figsize=(8, 6))
        axes = axes.flatten()
        
        for i, (df, ax) in enumerate(zip(df_list, axes)):
            ax.plot(df['Test Time (s)']/60, df[y_var], label=f'Channel {i+5}')
            ax.set_xlabel('Time (min)')
            ax.set_ylabel(y_var)
            ax.legend()
            ax.set_title(f'Channel {i+5} {title_suffix}')
            ax.grid()
        
        plt.tight_layout()


if __name__ == "__main__":
    main()
    plt.show()