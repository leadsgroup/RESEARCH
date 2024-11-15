import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import pandas as pd

# Read the noise data and census data
gdf_census = gpd.read_file('combined_tracts_data_all.geojson')

merged_gdf = gdf_census


# Filter the census data based on the bounding box
file_list = ['combined_results_per_tract.xlsx']
for file in file_list:
    data = pd.read_excel(file)

    # Create GeoDataFrames for noise data
    geometry = [Point(xy) for xy in zip(data['Longitude'], data['Latitude'])]
    data_gdf = gpd.GeoDataFrame(data, geometry=geometry, crs='EPSG:4326')

    # Perform a spatial join to combine the census geometries with the noise data
    merged_gdf = gpd.sjoin(merged_gdf, data_gdf, how='left', op='contains')

merged_gdf.to_file('processed.geojson', driver="GeoJSON") 