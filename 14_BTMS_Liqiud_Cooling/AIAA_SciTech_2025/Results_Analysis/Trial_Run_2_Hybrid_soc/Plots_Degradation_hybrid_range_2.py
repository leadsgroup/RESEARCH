import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def main(): 
    file_name = 'e_Twin_Otter_hybrid_range.xlsx'
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    df = read_excel_file(file_path)
    variables = [
        ('Cycle Day [power_bus]', 'Cell SOC (%) [power_bus]', 'LFP SOC vs Time'),
        ('Cycle Day [cruise_bus]', 'Cell SOC (%) [cruise_bus]', 'NMC SOC vs Time'),
        # ('Cycle Day', 'Resistance Growth', 'Resistance Growth vs Cycle Day'),
        # ('Cycle Day', 'Capacity Fade', 'Capacity Fade vs Cycle Day'),
        # ('Cycle Day', 'Cell Temperature (C)', 'Cell Temperature vs Cycle Day'),
        # ('Time (min)', 'Range (nmi)', 'Range vs Time'),
        # ('Time (min)', 'Module Heat Generation (W)', 'Module Heat Generation vs Time')
    ]
    plot_data(df, variables)
    return

def read_excel_file(file_path):
    # Read the Excel file
    excel_file = pd.ExcelFile(file_path)
    sheet_names = excel_file.sheet_names
    print("\nAvailable sheets in the Excel file:")
    
    # Dictionary to store concatenated data
    combined_data = {}
    last_time = 0  # Keep track of the last time value
    
    # Iterate through each sheet
    for sheet in sheet_names:
        print(f"- {sheet}")
        # Read the sheet into a DataFrame
        df = pd.read_excel(excel_file, sheet_name=sheet)
        
        # Add the last_time to the current sheet's time values
        if 'Time (min)' in df.columns:
            df['Time (min) [power_bus])'] = df['Time (min) [power_bus]'] + last_time
            # Update last_time for the next sheet
            last_time = df['Time (min)'].iloc[-1]
        
        # For each column in the current sheet
        for column in df.columns:
            # If column name already exists in combined_data, concatenate
            if column in combined_data:
                combined_data[column] = pd.concat([combined_data[column], df[column]], ignore_index=True)
            # If it's a new column name, add it to the dictionary
            else:
                combined_data[column] = df[column]
    
    # Convert the dictionary to a DataFrame
    result_df = pd.DataFrame(combined_data)
    return result_df

def plot_data(df, variables):
    for x_var, y_var, title_suffix in variables:
        plt.figure(figsize=(10, 6))
        
        plt.plot(df[x_var], df[y_var], label=y_var)
        
        plt.xlabel(x_var)
        plt.ylabel(y_var)
        plt.title(title_suffix)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
       

if __name__ == "__main__":
    main()
    plt.show()