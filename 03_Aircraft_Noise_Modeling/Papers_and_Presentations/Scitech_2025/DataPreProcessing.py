import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import shapely
from rtree import index
import numpy as np
import os

race_dictionary = {
    "B03002001": "Total",
    "B03002003": "White",
    "B03002004": "Black",
    "B03002005": "Native American",
    "B03002006": "Asian",
    "B03002007": "Pacific Islander",
    "B03002012": "Hispanic"
}

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


# def noise_preprocess(noise_file, census_tracts_file):
#     # Load the existing census tracts GeoDataFrame
#     census_tracts = census_tracts_file
#     try:
#         census_tracts = census_tracts.drop(['index_right'], axis=1)
#     except:
#         pass
#     # Load noise data
#     noise_data = noise_file

def noise_preprocess(noise_file, census_tracts_file):
    # Load the existing census tracts GeoDataFrame
    census_tracts = census_tracts_file

    # Load noise data
    noise_data = noise_file

    # Extract noise columns (excluding Latitude and Longitude)
    noise_columns = noise_data.columns[2:]

    # Compute minimum noise values for each column
    threshold_values = {col: noise_data[col].min() for col in noise_columns}
    
    # Convert points to a GeoDataFrame
    noise_data['geometry'] = noise_data.apply(
        lambda row: Point(row['Longitude'] - 360, row['Latitude']), axis=1
    )
    noise_gdf = gpd.GeoDataFrame(noise_data, geometry='geometry', crs=census_tracts.crs)

    # Perform spatial join to associate noise points with census tracts
    joined = gpd.sjoin(noise_gdf, census_tracts, how='inner', predicate='within')

    # Calculate mean noise levels for each census tract
    mean_noise = joined.groupby('index_right')[noise_columns].mean()

    # Merge the mean noise levels back into the census tracts GeoDataFrame
    census_tracts = census_tracts.join(mean_noise)

    # Fill missing values with threshold values
    for col in noise_columns:
        census_tracts[col] = census_tracts[col].fillna(threshold_values[col])

    return census_tracts


def sensitive_area_preprocess(file_list, census_tracts_file,file_names):
    # Load the existing census tracts GeoDataFrame
    demo_per_tract = census_tracts_file

    try:
        demo_per_tract = demo_per_tract.drop(['index_right'], axis=1)
    except:
        pass

    i=0
    for file in file_list:
        loc_data = file
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
        demo_per_tract[file_names[i]] = demo_per_tract.index.map(point_counts).fillna(0).astype(int)
        i+=1

    return demo_per_tract


def save_file(file,save_name):
    file.to_file(save_name, driver="GeoJSON")

def convert_csv(file,save_name):
    df = file.copy()
    df = df.drop(columns=['Latitude'])
    df = df.drop(columns=['Longitude'])
    df = df.drop(columns=['name_y'])
    df = df.drop(columns=['geoid_y'])
    df = df.loc[:, ~df.columns.str.contains('Error', case=False)]
    # df_proj = df.to_crs(epsg=3857)
    # df['Area']  = (df_proj.geometry.area)*(3.861e-7)
    df = df.drop(columns=['geometry'])

    combined_dict = {**race_dictionary, **income_dictionary}
    df = df.rename(columns=combined_dict)

    df.to_csv(save_name, index=False)

