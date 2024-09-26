import pandas as pd
import matplotlib.pyplot as plt


def main():
    file_paths = [
    'collins1_Channel_9_Wb_1.CSV',
    'collins1_Channel_9_Wb_2.CSV',
    'collins1_Channel_9_Wb_3.CSV',
    'collins1_Channel_9_Wb_4.CSV'
]
    plot_current_voltage(file_paths)
    return

def plot_current_voltage(file_paths):
    """
    Function to plot Current (A) and Voltage (V) against Test Time (s) from multiple CSV files.
    
    Args:
    file_paths (list): List of file paths for the CSV files.
    
    Returns:
    None: Displays the plot.
    """
    # Read all the CSV files into individual DataFrames
    dfs = [pd.read_csv(file) for file in file_paths]
    
    # Combine all the dataframes into a single dataframe
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Convert 'Test Time (s)' to numeric to avoid issues
    combined_df['Test Time (s)'] = pd.to_numeric(combined_df['Test Time (s)'], errors='coerce')
    
    # Create the plot
    plt.figure(figsize=(12, 6))
    
    # Plot Current
    plt.plot(combined_df['Test Time (s)'], combined_df['Current (A)'], label='Current (A)', color='blue')
    
    # Plot Voltage
    plt.plot(combined_df['Test Time (s)'], combined_df['Voltage (V)'], label='Voltage (V)', color='red')
    
    # Add labels and title
    plt.title('Current and Voltage vs Time')
    plt.xlabel('Test Time (s)')
    plt.ylabel('Current (A) / Voltage (V)')
    plt.legend()
    return

# ----------------------------------------------------------------------        
#   Call Main
# ----------------------------------------------------------------------    

if __name__ == '__main__':
    main()
    plt.show()

   
