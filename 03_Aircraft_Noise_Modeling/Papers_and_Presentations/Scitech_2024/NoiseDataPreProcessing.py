import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from rtree import index
import time
import numpy as np

ti=time.time()

# Read census tracts and noise data
census_tracts = gpd.read_file('combined_tracts_data_all.geojson')



def noise_preprocess(file_list):
    for file in file_list:
        # Load noise data
        noise_data = pd.read_csv(file)

        # Extract noise columns (excluding the first two columns for Latitude and Longitude)
        noise_columns = noise_data.columns.values.tolist()[2:]
        
        # Dictionary to store minimum noise levels for each column (used if no points are found within a tract)
        threshold_values = {column: noise_data[column].min() for column in noise_columns}

        # Create a spatial index for noise points
        spatial_index = index.Index()
        for idx, noise_point in noise_data.iterrows():
            spatial_index.insert(idx, (noise_point['Longitude'] - 360, noise_point['Latitude'], 
                                       noise_point['Longitude'] - 360, noise_point['Latitude']))

        # Initialize an empty DataFrame for the results per census tract
        noise_per_tract = pd.DataFrame(columns=['geometry'] + noise_columns)

        # Process each census tract
        for idx, tract in census_tracts.iterrows():
            tract_polygon = tract['geometry']

            # Get noise points within the current tract
            possible_matches_idx = list(spatial_index.intersection(tract_polygon.bounds))
            possible_matches = noise_data.iloc[possible_matches_idx]

            # Filter points within the tract and calculate the average noise for each column
            for column in noise_columns:
                noise_within_tract = possible_matches[
                    possible_matches.apply(lambda row: tract_polygon.contains(Point(row['Longitude'] - 360, row['Latitude'])), axis=1)]
                
                # Calculate the average noise level or use the threshold value if no points are found
                if not noise_within_tract.empty:
                    noise_level = noise_within_tract[column].mean()
                else:
                    noise_level = threshold_values[column]
                
                # Populate the average noise level and geometry for the current tract
                noise_per_tract.loc[idx, column] = noise_level
            noise_per_tract.loc[idx, 'geometry'] = tract_polygon

        # Convert to GeoDataFrame and save as GeoJSON
        noise_per_tract_gdf = gpd.GeoDataFrame(noise_per_tract, geometry='geometry')
        noise_per_tract_gdf.to_file(f"{file}_noise.geojson", driver="GeoJSON")

def sensitive_area_preprocess(file_list):
    demo_per_tract = census_tracts

    for file in file_list:
        # Read each CSV file
        loc_data = pd.read_csv(file+'.csv')

        # Create spatial index for the current file's points
        spatial_index = index.Index()
        for idx, point in loc_data.iterrows():
            spatial_index.insert(idx, (point['Longitude'], point['Latitude'], point['Longitude'], point['Latitude']))

        # Initialize a temporary DataFrame for the current file
        temp_df = pd.DataFrame(columns=['geometry'] + [file])
        
        # Iterate through each census tract
        for idx, tract in demo_per_tract.iterrows():
            tract_polygon = tract['geometry']

            # Check if the tract geometry is valid
            if tract_polygon.is_valid:
                # Get points within the current tract using spatial indexing
                possible_matches_idx = list(spatial_index.intersection(tract_polygon.bounds))
                possible_matches = loc_data.iloc[possible_matches_idx]

                # Filter points that are within the current tract
                points_within_tract = possible_matches[possible_matches.apply(
                    lambda row: tract_polygon.contains(Point(row['Longitude'], row['Latitude'])), axis=1)]

                # Store the count of points in the tract and its centroid location
                temp_df.loc[idx, 'geometry'] = tract_polygon
                temp_df.loc[idx, file] = len(points_within_tract)
        demo_per_tract = pd.merge(demo_per_tract, temp_df, on=['geometry'], how='left')
    
    # Save the combined results to one Excel file
    demo_per_tract.to_file('sensitive_results_tract.geojson', driver="GeoJSON") 

sensitive_list = ['Churches', 'Schools_Colleges_and_Universities']
sensitive_area_preprocess(sensitive_list)


file_name= ['TR_1000ft_LA_10min_All','SR_1000ft_LA_10min_All']
noise_preprocess(file_name)

to = time.time()

print((to-ti)/60)
