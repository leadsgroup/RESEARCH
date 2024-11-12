import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
import RCAIDE
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D


def main():
    # Load raw data and create the interpolant once
    interpolantnd, interpolantlin = create_interpolant()

    # Input parameters for voltage retrieval
    C_rate_input = 1.5          # Desired C-rate
    temperature_input = 25      # Desired initial temperature
    discharge_capacity_input = 0  # Desired discharge capacity

    # Get the interpolated voltage
    voltagend,voltagelin = get_voltage(interpolantnd, interpolantlin, C_rate_input, temperature_input, discharge_capacity_input)
    print(f"Interpolated Voltage Nearest : {voltagend}")
    print(f"Interpolated Voltage Linear : {voltagelin}")


    discharge_capacities = np.linspace(0, 2.6, 50)  # Replace with your actual data range
    temperatures = np.linspace(-20, 60, 20)       # Replace with your actual data range
    c_rates =  np.linspace(0.1,20, 200)         # Replace with your actual C-rates
    plot_voltage_surface(interpolantnd, discharge_capacities, temperatures, c_rates)
    plot_voltage_surface(interpolantlin, discharge_capacities, temperatures, c_rates)
    
    return

# Function to load data and create an interpolant
def create_interpolant():
    # Define the relative path to your data file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(current_dir, 'LFP_raw_data.res')

    # Load the raw_data using RCAIDE.load()
    raw_data = RCAIDE.load(full_path)

    # Prepare lists for the data needed for interpolation
    c_rates = []
    temperatures = []
    discharge_capacities = []
    voltages = []

    # Iterate through the data structure to populate the lists
    for c_rate_key, temp_data in raw_data.items():
        c_rate = float(c_rate_key)  # Convert C-rate to float

        for temp_key, data in temp_data.items():
            initial_temp = int(temp_key.split()[-1])  # Extract temperature as an integer

            # Extend lists with the discharge, voltage, and other data
            discharge_capacities.extend(data['discharge'])
            voltages.extend(data['voltage'])
            c_rates.extend([c_rate] * len(data['discharge']))
            temperatures.extend([initial_temp] * len(data['discharge']))

    # Convert lists to numpy arrays
    points = np.array([c_rates, temperatures, discharge_capacities]).T
    values = np.array(voltages)

    # Create the interpolant
    interpolantnd = NearestNDInterpolator(points, values)
    interpolantlin = LinearNDInterpolator(points,values)
    return interpolantnd, interpolantlin

# Define a function to get voltage using the created interpolant
def get_voltage(interpolantnd, interpolantlin, C_rate, temperature, discharge_capacity):
    voltagend = interpolantnd(C_rate, temperature, discharge_capacity)
    voltagelin = interpolantlin(C_rate, temperature, discharge_capacity)
    return voltagend,voltagelin

def plot_voltage_surface(interpolant, discharge_capacities, temperatures, c_rates):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create a dense grid for plotting the surface
    discharge_plot, temp_plot, c_rate_plot = np.meshgrid(discharge_capacities, temperatures, c_rates, indexing='ij')
    voltage_plot = interpolant(c_rate_plot, temp_plot, discharge_plot)

    # Plot slices through different C-rates for visualization clarity
    for c_rate in [0.1, 1, 2, 5,6]:  # Slice through selected C-rates
        voltage_slice = interpolant(c_rate, temp_plot[:, :, 0], discharge_plot[:, :, 0])
        # Add label to the surface for the legend
        surf = ax.plot_surface(discharge_plot[:, :, 0], temp_plot[:, :, 0], voltage_slice, alpha=0.7)
        # Set the label for the legend
        surf.set_label(f'C-rate={c_rate}')

    ax.set_xlabel('Discharge Capacity')
    ax.set_ylabel('Initial Temperature')
    ax.set_zlabel('Voltage')
    ax.set_title('Voltage Surface Plot Across Discharge Capacity, Temperature, and C-rate')
    ax.legend()  # This will now work correctly
    
    return


if __name__ == "__main__":
    main()
    plt.show()

