import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator,griddata,make_interp_spline
from mpl_toolkits.mplot3d import Axes3D
import os
def load_data(file_path, sheet_name):
    """
    Load data from an Excel sheet.

    :param file_path: Path to the Excel file
    :param sheet_name: Sheet name to load data from
    :return: Pandas DataFrame
    """
    return pd.read_excel(file_path, sheet_name=sheet_name)

def create_surface_plot(data):
    """
    Create a 3D surface plot for HAS power, HEX power, and Cycle Day, and highlight the maxima.

    :param data: Pandas DataFrame containing HAS_power, HEX_power, and Cycle Day
    """
    # Extract relevant columns
    x = data['HAS_power']
    y = data['HEX_power']
    z = data['Cycle Day']

    # Define the grid for interpolation
    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate Z values for the grid using LinearNDInterpolator
    interpolator = LinearNDInterpolator((x, y), z)
    zi = interpolator(xi, yi)

    # Find the maximum value and its coordinates
    max_value = np.nanmax(zi)  # Handle NaN values
    max_index = np.unravel_index(np.nanargmax(zi), zi.shape)  # Get indices of max value
    max_coords = (xi[max_index], yi[max_index])  # Get x, y coordinates

    # Create the 3D surface plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    surf = ax.plot_surface(xi, yi, zi, cmap='viridis', edgecolor='none', alpha=0.8)

    # Highlight the maximum point
    ax.scatter(max_coords[0], max_coords[1], max_value, color='red', s=100, label='Max Cycle Day')
    ax.text(max_coords[0], max_coords[1], max_value + 5,  # Offset text slightly above the point
            f'Max: {max_value:.2f}\n({max_coords[0]:.2f}, {max_coords[1]:.2f})',
            color='red', fontsize=10)

    # Add a color bar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
    cbar.set_label('Cycle Day', fontsize=12)

    # Add labels and title
    ax.set_title('Surface Plot: HAS Power vs HEX Power vs Cycle Day', fontsize=16)
    ax.set_xlabel('HAS Power', fontsize=12)
    ax.set_ylabel('HEX Power', fontsize=12)
    ax.set_zlabel('Cycle Day', fontsize=12)

    # Add legend
    ax.legend()

    # Adjust layout and show the plot
    plt.tight_layout()

def weight_plot(data):
    """
    Plot Cycle Day vs TMS Weight with smooth interpolation.

    :param data: Pandas DataFrame containing 'TMS Weight' and 'Cycle Day' columns
    """
    # Extract relevant columns
    x = data['TMS Weight'].values
    x1 = data['HAS_power'].values
    x2= data['HEX_power'].values
    x3= data['RES_dimensions'].values
    y = data['Cycle Day'].values

    # Sort the data by x for proper interpolation
    sorted_indices1 = np.argsort(x1)
    sorted_indices2 = np.argsort(x2)
    sorted_indices3 = np.argsort(x3)
   # x = x[sorted_indices]
    x1 = x1[sorted_indices1]
    x2 = x2[sorted_indices2]
    x3 = x3[sorted_indices3]
    y1 = y[sorted_indices1]
    y2= y[sorted_indices2]
    y3 = y[sorted_indices3]
    #Create a dense grid of points for smooth interpolation
    x_smooth_1 = np.linspace(x1.min(), x1.max(), 300)
    x_smooth_2 = np.linspace(x2.min(), x2.max(), 300)
    x_smooth_3 = np.linspace(x3.min(), x3.max(), 300)

     # Use cubic spline interpolation to generate smooth y values
    spline_1 = make_interp_spline(x1, y, k=1)  
    y_smooth_1 = spline_1(x_smooth_1)
    spline_2 = make_interp_spline(x2, y, k=1)  
    y_smooth_2 = spline_2(x_smooth_2)
    spline_3 = make_interp_spline(x3, y, k=1)  
    y_smooth_3 = spline_3(x_smooth_3)

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot the original data points
    #ax.scatter(x, y, color='red', label='Original Data', alpha=0.7)


    # Plot the smooth interpolation
    #ax.plot(y,x, linestyle='-', color='blue')
    ax.plot(x_smooth_1,y_smooth_1,linestyle='-')
    ax.plot(x_smooth_2,y_smooth_2,linestyle='-')
    ax.plot(x_smooth_3,y_smooth_3,linestyle='-')

    # Add labels, title, and legend
    ax.set_title('Cycle Day vs TMS Weight)', fontsize=16, weight='bold')
    ax.set_xlabel('Cycle Day', fontsize=14)
    ax.set_ylabel('All Vairables', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)

    # Adjust layout and show the plot
    plt.tight_layout()

def create_3d_surface(data):
    """
    Create a 3D surface plot to show Cycle Day as a function of HAS Power, HEX Power, and RES Dimension.

    :param data: Pandas DataFrame containing HAS_power, HEX_power, RES_dimensions, and Cycle Day
    """
    # Extract relevant columns
    x = data['HAS_power']
    y = data['HEX_power']
    z = data['RES_dimensions']
    cycle_day = data['Cycle Day']

    # Define grid for interpolation
    xi = np.linspace(x.min(), x.max(), 50)
    yi = np.linspace(y.min(), y.max(), 50)
    zi = np.linspace(z.min(), z.max(), 50)
    xi, yi, zi = np.meshgrid(xi, yi, zi)

    # Interpolate Cycle Day values onto the grid
    cycle_day_grid = griddata(
        points=(x, y, z),
        values=cycle_day,
        xi=(xi, yi, zi),
        method='linear'
    )

    # Create 3D surface plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    surf = ax.plot_surface(
        xi[:, :, 0], yi[:, :, 0], cycle_day_grid[:, :, 0],
        cmap='viridis', edgecolor='none', alpha=0.8
    )

    # Add color bar to represent Cycle Day
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
    cbar.set_label('Cycle Day', fontsize=12)

    # Add labels and title
    ax.set_title('3D Surface Plot: Cycle Day as a Function of Inputs', fontsize=16)
    ax.set_xlabel('HAS Power', fontsize=12)
    ax.set_ylabel('HEX Power', fontsize=12)
    ax.set_zlabel('Cycle Day', fontsize=12)

    # Show plot
    plt.tight_layout()



def main():
    """
    Main function to load data, process it, and plot the surface.
    """
    # Load the data
    file_name = 'consolidated_exit_conditions.xlsx'
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    sheet_name = 'Sheet1'
    data = load_data(file_path, sheet_name)

    # Add a condition for marking exit condition
    if 'Exit Condition' in data.columns:
        data['Color'] = data['Exit Condition'].apply(
            lambda condition: 'red' if 'Temperature' in condition else 'blue'
        )

    # Create the surface plot
    weight_plot(data)
    #create_surface_plot(data)
    #create_3d_surface(data)

# Entry point for the script
if __name__ == "__main__":
    main()
    plt.show()
