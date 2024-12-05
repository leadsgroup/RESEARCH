import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
import shapely
import os
import matplotlib.pyplot as plt



race_dictionary = {
    "Total": "B03002001",
    "White": "B03002003",
    "Black or African American": "B03002004",
    "American Indian and Alaska Native": "B03002005",
    "Asian": "B03002006",
    "Native Hawaiian and Other Pacific Islander": "B03002007",
    "Hispanic or Latino": "B03002012"
}

income_dictionary = {
    "Total": "B19001001",
    "Less than $10,000": "B19001002",
    "$10,000 to $14,999": "B19001003",
    "$15,000 to $19,999": "B19001004",
    "$20,000 to $24,999": "B19001005",
    "$25,000 to $29,999": "B19001006",
    "$30,000 to $34,999": "B19001007",
    "$35,000 to $39,999": "B19001008",
    "$40,000 to $44,999": "B19001009",
    "$45,000 to $49,999": "B19001010",
    "$50,000 to $59,999": "B19001011",
    "$60,000 to $74,999": "B19001012",
    "$75,000 to $99,999": "B19001013",
    "$100,000 to $124,999": "B19001014",
    "$125,000 to $149,999": "B19001015",
    "$150,000 to $199,999": "B19001016",
    "$200,000 or More": "B19001017"
}



def community_annoyance(gdf_data, noise_sensitive_structures, sensitivity_levels):
    # Read the census data (GeoDataFrame)
    gdf_census = gdf_data

    total_population = gdf_census['B03002001']  # Total population in census tract
    #convert to projected coordinates to get accurate area
    gdf_census_projected = gdf_census.to_crs(epsg=3857)
    tract_area = (gdf_census_projected.geometry.area)*(3.861e-7)
    
    population_density = (total_population / tract_area) # Population density (per unit area)

    temp_gdf = gdf_census.copy()

    # Calculate net noise-sensitive structure influence for each structure
    for structure in noise_sensitive_structures:
        # Net influence based on sensitivity level
        temp_gdf[f'{structure}s_net'] = (gdf_census[structure] * sensitivity_levels[structure])/gdf_census[structure]

    # Summation of all net influences
    temp_gdf['summation'] = temp_gdf[[f'{structure}s_net' for structure in noise_sensitive_structures]].sum(axis=1)

    # Community annoyance calculation

    gdf_census['L_dnlog(PD)'] = (gdf_census['L_dn'] * (np.log10(population_density)))
    gdf_census['CA'] = (gdf_census['L_dn'] * (np.log10(population_density)) * np.log(np.exp(1)+temp_gdf['summation']))

    return gdf_census


def structure_analysis(gdf_file,column_names,noise_type,noise_threshold,min_lon,min_lat,max_lon,max_lat):
    gdf_census = gdf_file

    # Define the folder path for saving plots
    folder_path = 'Bar_Plots'

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

    gdf_census = gdf_census[gdf_census[noise_type] >= noise_threshold]



    values = gdf_census[column_names].sum()
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
    plt.rcParams.update(parameters)
    #  Bar plot
    plt.bar(column_names, values, width = 0.5)
    plot_file_name = os.path.join(folder_path, f'{noise_threshold}dB_Plot.png')
    plt.savefig(plot_file_name) 
    plt.close()



def community_analysis(gdf_file,column_names,noise_type,noise_threshold,min_lon,min_lat,max_lon,max_lat):
    gdf_census = gdf_file

    # Define the folder path for saving plots
    folder_path = 'Bar_Plots'

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

    gdf_census = gdf_census[gdf_census[noise_type] >= noise_threshold]



    values = gdf_census[column_names].sum()/gdf_census['B03002001'].sum()*100
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
    plt.rcParams.update(parameters)
    #  Bar plot
    plt.bar(column_names, values, width = 0.5)
    plot_file_name = os.path.join(folder_path, f'{noise_threshold}dB_Plot.png')
    plt.savefig(plot_file_name) 
    plt.close()
