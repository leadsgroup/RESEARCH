import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import pandas as pd
import os
import numpy as np

def impact_plotting(gdf_file, column_names,min_lon,min_lat,max_lon,max_lat):
    gdf_census = gdf_file

    # Define the folder path for saving plots
    folder_path = 'Plots'

    # Check if the folder exists, and if not, create it
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)

    # Create a Polygon object for the bounding box
    bbox = Polygon([(min_lon, min_lat), 
                    (max_lon, min_lat), 
                    (max_lon, max_lat), 
                    (min_lon, max_lat), 
                    (min_lon, min_lat)])

    # filter the GeoDataFrame based on bounds
    gdf_census = gdf_census[gdf_census.geometry.intersects(bbox)] #can change form within to intersects

    # Plot choropleth maps for each data column
    for data_column in column_names:
        plt.figure(figsize=(12, 12))
        plt.rcParams['axes.linewidth'] = 2.0
        plt.rcParams["font.family"] = "Times New Roman"
        parameters = {
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'axes.titlesize': 0,
            'xtick.major.pad': 1,
            'ytick.major.pad': 0,
        }
        
        # Plot the GeoDataFrame
        ax = gdf_census.plot(
            column=data_column,  # Use the data column
            cmap='coolwarm',  # Color map
            edgecolor='white',
            linewidth=0.15,
            alpha=0.75,
            legend=True,
            legend_kwds={
                "label": '(dB)',  # label for the legend
                "orientation": "vertical",
                "shrink": 0.75,  # Adjust the color bar size
                "ticks": np.linspace(45,90,6) # Specify fewer bins
            },
            vmin=45,  # Minimum value for color scaling
            vmax=90,  # Maximum value for color scaling
            ax=None
        )

        # Set map limits and labels
        ax.set_ylim(min_lat,max_lat)
        ax.set_xlim(min_lon,max_lon)
        ax.set_aspect('auto')
        plt.ylabel('Latitude ($\degree$)')
        plt.xlabel('Longitude ($\degree$)')

        # Save the plot as an image
        plot_file_name = os.path.join(folder_path, f'{data_column}.png')
        plt.savefig(plot_file_name)  # Save with higher DPI
        plt.close()

def census_plotting(gdf_file, column_names,min_lon,min_lat,max_lon,max_lat):
    gdf_census = gdf_file

    # Define the folder path for saving plots
    folder_path = 'Plots'

    # Check if the folder exists, and if not, create it
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)

    # Create a Polygon object for the bounding box
    bbox = Polygon([(min_lon, min_lat), 
                    (max_lon, min_lat), 
                    (max_lon, max_lat), 
                    (min_lon, max_lat), 
                    (min_lon, min_lat)])

    # filter the GeoDataFrame based on bounds
    gdf_census = gdf_census[gdf_census.geometry.intersects(bbox)] #can change form within to intersects

    # Plot choropleth maps for each data column
    for data_column in column_names:
        plt.figure(figsize=(12, 12))
        plt.rcParams['axes.linewidth'] = 2.0
        plt.rcParams["font.family"] = "Times New Roman"
        parameters = {
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'axes.titlesize': 0,
            'xtick.major.pad': 1,
            'ytick.major.pad': 0,
        }
        
        # Plot the GeoDataFrame
        ax = gdf_census.plot(
            column=data_column,  # Use the data column
            cmap='coolwarm',  # Color map
            edgecolor='white',
            linewidth=0.15,
            alpha=0.75,
            legend=True,
            legend_kwds={
                "label": '(dB)',  # label for the legend
                "orientation": "vertical",
                "shrink": 0.75  # Adjust the color bar size
            },
            ax=None
        )

        # Set map limits and labels
        ax.set_ylim(min_lat,max_lat)
        ax.set_xlim(min_lon,max_lon)
        ax.set_aspect('auto')
        plt.ylabel('Latitude ($\degree$)')
        plt.xlabel('Longitude ($\degree$)')

        # Save the plot as an image
        plot_file_name = os.path.join(folder_path, f'{data_column}.png')
        plt.savefig(plot_file_name)  # Save with higher DPI
        plt.close()