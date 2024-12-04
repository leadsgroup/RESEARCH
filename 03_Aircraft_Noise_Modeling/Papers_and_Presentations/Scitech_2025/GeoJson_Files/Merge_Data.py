import geopandas as gpd
import pandas as pd

# # Load the data
# gdf_censusla1 = gpd.read_file('Race/LosAngeles.geojson')
# gdf_censusla2 = gpd.read_file('Income/LosAngeles.geojson')


# gdf_censusoc1 = gpd.read_file('Race/Orange.geojson')
# gdf_censusoc2 = gpd.read_file('Income/Orange.geojson')


# gdf_censussb1 = gpd.read_file('Race/SanB.geojson')
# gdf_censussb2 = gpd.read_file('Income/SanB.geojson')

# gdf_censusriv1 = gpd.read_file('Race/Riverside.geojson')
# gdf_censusriv2 = gpd.read_file('Income/Riverside.geojson')

# gdf1  = gpd.sjoin(gdf_censusla1,gdf_censusla2)
# gdf2  = gpd.sjoin(gdf_censusoc1,gdf_censusoc2)
# gdf3  = gpd.sjoin(gdf_censussb1,gdf_censussb2)
# gdf4  = gpd.sjoin(gdf_censusriv1,gdf_censusriv2)

# # gdf1  = gdf_censusla1.merge(gdf_censusla2, on='geoid')
# # gdf2  = gdf_censusoc1.merge(gdf_censusoc2, on='geoid')
# # gdf3  = gdf_censussb1.merge(gdf_censussb2, on='geoid')
# # gdf4  = gdf_censusriv1.merge(gdf_censusriv2, on='geoid')


# # Concatenate the GeoDataFrames
# combined_gdf =(pd.concat([gdf1, gdf2, gdf3, gdf4], 
#                                           ignore_index=True))
# combined_gdf = combined_gdf.drop(columns=['index_right'])
# combined_gdf.to_file('LA_Area_Test.geojson', driver="GeoJSON")

files = [
    "Race/LosAngeles.geojson",
    "Income/LosAngeles.geojson",
    "Race/Orange.geojson",
    "Income/Orange.geojson",
    "Race/SanB.geojson",
    "Income/SanB.geojson",
    "Race/Riverside.geojson",
    "Income/Riverside.geojson",
]


gdfs = [gpd.read_file(file) for file in files]
common_crs = gdfs[0].crs

merged_gdfs = []
for i in range(0, len(gdfs), 2):  # Process in pairs (race, income)
    gdf_race = gdfs[i]
    gdf_income = gdfs[i + 1]
    
    # Merge the race and income data based on a common attribute, e.g., 'GEOID'
    merged_gdf = gdf_race.merge(gdf_income, on="name", how="inner") 
    merged_gdfs.append(merged_gdf)

# Concatenate all merged GeoDataFrames into one
if 'geometry' in merged_gdf.columns:
    merged_gdf = gpd.GeoDataFrame(merged_gdf, geometry='geometry', crs=common_crs)
else:
    # If geometry is missing, reassign it from one of the original dataframes
    merged_gdf = gpd.GeoDataFrame(merged_gdf, geometry=gdf_race.geometry, crs=common_crs)

final_gdf = pd.concat(merged_gdfs, ignore_index=True)



# Save the final combined GeoDataFrame to a new GeoJSON file
final_gdf.to_file('LA_Area_Tract.geojson', driver="GeoJSON")


