import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon

def community_annoyance(gdf_data, noise_sensitive_structures, sensitivity_levels):
    # Read the census data (GeoDataFrame)
    gdf_census = gpd.read_file(gdf_data)

    # Define the bounding box coordinates (optional, can change from within to intersects)
    min_lon, max_lon = -119, -117.3
    min_lat, max_lat = 33.6, 34.4

    # Create a Polygon object for the bounding box
    bbox = Polygon([(min_lon, min_lat), 
                    (max_lon, min_lat), 
                    (max_lon, max_lat), 
                    (min_lon, max_lat), 
                    (min_lon, min_lat)])

    # Filter census tracts within the bounding box (can also use `within` instead of `intersects`)
    gdf_census = gdf_census[gdf_census.geometry.intersects(bbox)]

    # Compute population density (vectorized)
    total_population = gdf_census['B03002001']  # Total population in census tract
    tract_area = gdf_census['geometry'].area  # Area of the census tract
    print(tract_area[0:20])
    population_density = total_population / tract_area  # Population density (per unit area)

    print(population_density)

    temp_gdf = gdf_census.copy()
    noise_sens_struct_sum = 0

    for structure in noise_sensitive_structures:
        noise_sens_struct_sum = gdf_census[structure]
    
        noise_sens_struct_density = noise_sens_struct_sum/tract_area

        log_noise_sens_struct = np.nan_to_num(np.log(noise_sens_struct_density),nan=0)

        temp_gdf[structure+'s_net'] = log_noise_sens_struct *  sensitivity_levels[structure]

    temp_gdf['summation'] = sum(temp_gdf[structure + 's_net'] for structure in noise_sensitive_structures)

    print(temp_gdf['summation'][:20])
    density_calc = np.nan_to_num(np.log(population_density),nan=0)
    print(density_calc[:20])
    c_a = gdf_census['L_dn'] * np.nan_to_num(np.log(population_density),nan=0) * temp_gdf['summation']

    gdf_census['Community Annoyance'] = c_a

    # Save the updated GeoDataFrame to a GeoJSON file
    # gdf_census.to_file('CA'+gdf_data, driver='GeoJSON')

    return c_a  


sensitivity_levels = {
"Schools_Colleges_and_Universities": .6,
"Churches": .7
}

sensitive_list = ['Churches', 'Schools_Colleges_and_Universities']
file = 'TRcombined_tracts_data_all.geojson'

test = community_annoyance(file,sensitive_list,sensitivity_levels)

# print(np.where(test>0))