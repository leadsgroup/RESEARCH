import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon
import os
import matplotlib.pyplot as plt



race_dictionary = {
    "Total": "B03002001",
    "White": "B03002003",
    "Black": "B03002004",
    "Native American": "B03002005",
    "Asian": "B03002006",
    "Pacific Islander": "B03002007",
    "Hispanic": "B03002012"
}


race_dictionary_flipped = dict([(value, key) for key, value in race_dictionary.items()])

income_dictionary = {
    "B19001001": "Total Households",
    "B19001002": "<$10k",
    "B19001003": "$10k to $15k",
    "B19001004": "$15k to $20k",
    "B19001005": "$20k to $25k",
    "B19001006": "$25k to $30k",
    "B19001007": "$30k to $35k",
    "B19001008": "$35k to $40k",
    "B19001009": "$40k to $45k",
    "B19001010": "$45k to $50k",
    "B19001011": "$50k to $60k",
    "B19001012": "$60k to $75k",
    "B19001013": "$75k to $100k",
    "B19001014": "$100k to $125k",
    "B19001015": "$125k to $150",
    "B19001016": "$150k to $200k",
    "B19001017": "$200k+ "
}


income_aggregation2 = {
    '<50K': [
        "B19001002",
        "B19001003",
        "B19001004",
        "B19001005",
        "B19001006",
        "B19001007",
        "B19001008",
        "B19001009",
        "B19001010",
    ],
    '50k to 100k': [
        "B19001011",
        "B19001012",
        "B19001013",
    ],
    '100k to 150k': [
        "B19001014",
        "B19001015",
    ],
    '150k to 200k': [
        "B19001016",
    ],
    '200k+': [
        "B19001017",
    ],
}


income_dictionary_flipped = dict([(value, key) for key, value in income_dictionary.items()])


#Census Tract 469, Riverside, CA large area negative value

def community_annoyance(gdf_data, noise_sensitive_structures, sensitivity_levels):
    # Read the census data (GeoDataFrame)
    gdf_census = gdf_data
    race_list = list(race_dictionary.values())
    # income_list = list(income_dictionary.keys())

    income_columns =['<50K','50k to 100k','100k to 150k','150k to 200k','200k+']
    for income in income_columns:
        gdf_data[income]= gdf_data[income_aggregation2[income]].sum(axis=1)

    #convert to projected coordinates to get accurate area
    gdf_census_projected = gdf_census.to_crs(epsg=3857)
    tract_area = (gdf_census_projected.geometry.area)*(3.861e-7)
    gdf_census['area sq mi'] = tract_area
    gdf_census['total_structures'] = gdf_census[noise_sensitive_structures].sum(axis=1)


    for col in race_list:
        total_population = gdf_census[col]  # Total population in census tract
        
        population_density = (total_population / tract_area) # Population density (per unit area)

        temp_gdf = gdf_census.copy()

        # # Calculate net noise-sensitive structure influence for each structure
        # for structure in noise_sensitive_structures:
        #     # Net influence based on sensitivity level
        #     temp_gdf[f'{structure}s_net'] = (gdf_census[structure] * sensitivity_levels[structure])/gdf_census[structure]

        # # Summation of all net influences
        # temp_gdf['summation'] = temp_gdf[[f'{structure}s_net' for structure in noise_sensitive_structures]].sum(axis=1)

        # Initialize a column for the net influence of noise-sensitive structures
        temp_gdf['net_influence'] = 0

        # Loop through each noise-sensitive structure
        for structure in noise_sensitive_structures:
            # Calculate weighted influence for the structure
            weighted_influence = gdf_census[structure] * sensitivity_levels[structure]
            
            # Add the weighted influence to the net influence column
            temp_gdf['net_influence'] += weighted_influence

        # Normalize by the total structures in each tract
        temp_gdf['s_net'] = 1+ np.nan_to_num(temp_gdf['net_influence'] / gdf_census['total_structures'])
        # gdf_census['s_net'] = 1+ np.nan_to_num(temp_gdf['net_influence'] / gdf_census['total_structures'])
        # gdf_census['lastterm'] = (np.nan_to_num(np.log(1+gdf_census['total_structures']/tract_area))+1)


        gdf_census[f'{race_dictionary_flipped[col]}_LogPOP'] = (gdf_census['L_dn'] * (np.log10(population_density)))

        gdf_census[f'{race_dictionary_flipped[col]}_CA'] = (gdf_census['L_dn'] * np.nan_to_num(np.log10(population_density)) * temp_gdf['s_net']*(np.nan_to_num(np.log(1+gdf_census['total_structures']/tract_area))+1))



    
    for col in income_columns:
        total_population = gdf_census[col]  # Total population in census tract
        
        population_density = (total_population / tract_area) # Population density (per unit area)

        temp_gdf = gdf_census.copy()
        # gdf_census['total_structures'] = gdf_census[noise_sensitive_structures].sum(axis=1)

        # Initialize a column for the net influence of noise-sensitive structures
        temp_gdf['net_influence'] = 0

        # Loop through each noise-sensitive structure
        for structure in noise_sensitive_structures:
            # Calculate weighted influence for the structure
            weighted_influence = gdf_census[structure] * sensitivity_levels[structure]
            
            # Add the weighted influence to the net influence column
            temp_gdf['net_influence'] += weighted_influence

        # Normalize by the total structures in each tract
        temp_gdf['s_net'] = 1+ np.nan_to_num(temp_gdf['net_influence'] / gdf_census['total_structures'])

        gdf_census[f'{col}_LogPOP'] = (gdf_census['L_dn'] * (np.log10(population_density)))
        
        gdf_census[f'{col}_CA'] = (gdf_census['L_dn'] * np.nan_to_num(np.log10(population_density)) * temp_gdf['s_net']*(np.nan_to_num(np.log(1+gdf_census['total_structures']/tract_area))+1))

    
        # # Calculate net noise-sensitive structure influence for each structure
        # for structure in noise_sensitive_structures:
        #     # Net influence based on sensitivity level
        #     temp_gdf[f'{structure}s_net'] = (gdf_census[structure] * sensitivity_levels[structure])/gdf_census[structure]

        # # Summation of all net influences
        # temp_gdf['summation'] = temp_gdf[[f'{structure}s_net' for structure in noise_sensitive_structures]].sum(axis=1)

        # # Community annoyance calculation
        # # gdf_census[f'{col}'] = (gdf_census['L_dn'] * (np.log10(population_density)))
        # gdf_census[f'{col}'] = (gdf_census['L_dn'] * (np.log10(population_density)) * np.log(np.exp(1)+temp_gdf['summation']))


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
