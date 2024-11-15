import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from rtree import index
import time
import numpy as np

ti=time.time()

def noise_preprocess(file_list, census_tracts_file):
    # Load the existing census tracts GeoDataFrame
    census_tracts = gpd.read_file(census_tracts_file)

    for file in file_list:
        # Load noise data
        noise_data = pd.read_csv(f'Noise_Data/{file}.csv')

        # Extract noise columns (excluding the first two columns for Latitude and Longitude)
        noise_columns = noise_data.columns.values.tolist()[2:]

        # Dictionary to store minimum noise levels for each column (used if no points are found within a tract)
        threshold_values = {column: noise_data[column].min() for column in noise_columns}

        # Convert noise data to GeoDataFrame with adjusted longitude
        noise_data['geometry'] = noise_data.apply(
            lambda row: Point(row['Longitude'] - 360, row['Latitude']), axis=1
        )
        noise_gdf = gpd.GeoDataFrame(noise_data, geometry='geometry', crs=census_tracts.crs)

        # Perform a spatial join to associate noise points with census tracts
        joined = gpd.sjoin(noise_gdf, census_tracts, how='inner', predicate='within')

        # Group by census tract and calculate the mean noise levels for each column
        for column in noise_columns:
            mean_noise = joined.groupby('index_right')[column].mean()

            # Add a new column for the calculated mean noise levels to the census tracts GeoDataFrame
            census_tracts[column] = census_tracts.index.map(mean_noise).fillna(threshold_values[column])

    # Save the updated GeoDataFrame back to a GeoJSON file
    census_tracts.to_file(census_tracts_file, driver="GeoJSON")


def sensitive_area_preprocess(file_list, census_tracts_file):
    # Load the existing census tracts GeoDataFrame
    demo_per_tract = gpd.read_file(census_tracts_file)

    for file in file_list:
        # Read each CSV file
        loc_data = pd.read_csv(f'Raw_Data/{file}.csv')

        # Convert points to a GeoDataFrame
        loc_data['geometry'] = loc_data.apply(
            lambda row: Point(row['Longitude'], row['Latitude']), axis=1
        )
        loc_gdf = gpd.GeoDataFrame(loc_data, geometry='geometry', crs=demo_per_tract.crs)

        # Perform a spatial join to associate points with census tracts
        joined = gpd.sjoin(loc_gdf, demo_per_tract, how='inner', predicate='within')

        # Count points within each census tract
        point_counts = joined.groupby('index_right').size()

        # Add a new column for the point counts to the census tracts GeoDataFrame
        demo_per_tract[file] = demo_per_tract.index.map(point_counts).fillna(0).astype(int)

    # Save the updated GeoDataFrame back to a GeoJSON file
    demo_per_tract.to_file(census_tracts_file, driver="GeoJSON")

sensitive_list = ['Churches', 'Schools_Colleges_and_Universities']
sensitive_area_preprocess(sensitive_list,'combined_tracts_data_all.geojson')


file_name= ['TR_1000ft_LA_10min_All']
noise_preprocess(file_name,'combined_tracts_data_all.geojson')

to = time.time()

print((to-ti)/60)
