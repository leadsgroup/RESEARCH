import pandas as pd
import matplotlib.pyplot as plt


def main():
    file_paths = [
    'collins1_Channel_9_Wb_1.CSV',
    'collins1_Channel_9_Wb_2.CSV',
    'collins1_Channel_9_Wb_3.CSV',
    'collins1_Channel_9_Wb_4.CSV'
]
    current, voltage, time =  plot_current_voltage(file_paths)
    #plot_hppc_only(current, voltage, time)
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
    #plt.figure(figsize=(12, 6))
    
    # Plot Current
    df.plot(combined_df['Test Time (s)'], combined_df['Current (A)'], label='Current (A)', color='blue')
    
    # Plot Voltage
    df.plot(combined_df['Test Time (s)'], combined_df['Voltage (V)'], label='Voltage (V)', color='red')
    
    # Add labels and title
    #plt.title('Current and Voltage vs Time')
    #plt.xlabel('Test Time (s)')
    #plt.ylabel('Current (A) / Voltage (V)')
    #plt.legend()
    return combined_df['Current (A)'], combined_df['Voltage (V)'],combined_df['Test Time (s)']

def plot_hppc_only(current, voltage, time):
    """
    Function to identify and plot only the time intervals in which HPPC tests are defined.
    
    Args:
    
    
    Returns:
    None: Displays the plot with HPPC test intervals.
    """
    
    # Identify HPPC tests by detecting significant changes in current
    current_diff = combined_df['Voltage (V)'].diff() 
    
    # Define a threshold for significant changes in current to identify HPPC tests
    hppc_threshold = current_diff.abs().mean() * 80 # Adjust this value based on your dataset
    
    
    # Create a boolean mask that indicates periods of significant current change
    hppc_mask = current_diff > hppc_threshold 

    # Add a buffer around these events (e.g., include a few time steps before and after to capture the full pulse)
    combined_df['HPPC'] = hppc_mask.astype(int).rolling(window=5, min_periods=1).max()

    
    # Filter out only the HPPC test periods
    hppc_tests = combined_df[combined_df['HPPC'] == 1]
    
    # Plot only the HPPC test intervals if any are found
    if not hppc_tests.empty:
        plt.figure(figsize=(12, 6))
        
        # Plot Current during HPPC test intervals
        plt.plot(hppc_tests['Test Time (s)'], hppc_tests['Current (A)'], label='HPPC Current (A)', color='blue')
        
        # Plot Voltage during HPPC test intervals
        plt.plot(hppc_tests['Test Time (s)'], hppc_tests['Voltage (V)'], label='HPPC Voltage (V)', color='red')
        
        # Add labels and title
        plt.title('HPPC Current and Voltage vs Time (Filtered HPPC Intervals)')
        plt.xlabel('Test Time (s)')
        plt.ylabel('Current (A) / Voltage (V)')
        plt.legend()
        plt.grid(True)
        
        # Show the plot
        plt.show()
    else:
        print("No HPPC test intervals identified based on the current threshold.")


# ----------------------------------------------------------------------        
#   Call Main
# ----------------------------------------------------------------------    

if __name__ == '__main__':
    main()
    plt.show()

   
