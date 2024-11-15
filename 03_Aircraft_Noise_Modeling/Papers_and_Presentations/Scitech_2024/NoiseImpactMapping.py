import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
# Read the census data
gdf_census = gpd.read_file('combined_tracts_data_all.geojson')


# Define the folder path for saving plots
folder_path = 'Plots'

# Check if the folder exists, and if not, create it
if not os.path.isdir(folder_path):
    os.mkdir(folder_path)

# Columns to plot
column_names = ['L_Amax', 'L_AeqT', 'L_AeqT_24hr', 'SEL', 'L_dn', 'L_Aeq_jetliner']

# Plot choropleth maps for each data column
for data_column in column_names:
    plt.figure(figsize=(12, 12))
    plt.rcParams['axes.linewidth'] = 2.0
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.major.pad': 1,
        'ytick.major.pad': 0,
    }
    plt.rcParams.update(parameters)
    
    # Plot the GeoDataFrame
    ax = gdf_census.plot(
        column=data_column,  # Use the data column
        cmap='plasma',  # Color map
        edgecolor='white',
        linewidth=0.15,
        alpha=0.75,
        legend=True,
        legend_kwds={
            "label": f'{data_column} (dB)',  # Dynamic label for the legend
            "orientation": "vertical",
            "shrink": 0.75,  # Adjust the color bar size
            "ticks": np.linspace(45, 90, num=6)  # Specify fewer bins
        },
        vmin=45,  # Minimum value for color scaling
        vmax=90,  # Maximum value for color scaling
        ax=None
    )

    # Set map limits and labels
    ax.set_ylim(33.6, 34.4)
    ax.set_xlim(-119, -117)
    plt.title(f'{data_column}')
    plt.ylabel('Latitude ($\degree$)')
    plt.xlabel('Longitude ($\degree$)')

    # Save the plot as an image
    plot_file_name = os.path.join(folder_path, f'{data_column}.png')
    plt.savefig(plot_file_name)  # Save with higher DPI
    plt.close()