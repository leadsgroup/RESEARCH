import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import pandas as pd
import os
import numpy as np
# Read the census data
gdf_census = gpd.read_file('combined_tracts_data_all.geojson')

# Define the bounding box coordinates
min_lon, max_lon = -119, -117
min_lat, max_lat = 33, 35

# Filter the census data based on bounding box
bbox_polygon = Polygon([(min_lon, min_lat), (max_lon, min_lat), (max_lon, max_lat), (min_lon, max_lat)])
gdf_census_inside_bbox = gdf_census[gdf_census['geometry'].intersects(bbox_polygon)]

types = ['SR','TR']
timesteps = ['10']

# Define the folder path for saving plots
folder_path = 'newstuff'

# Check if the folder exists, and if not, create it
if not os.path.isdir(folder_path):
    os.mkdir(folder_path)
    
file_name_noise = ['TR_1000ft_LA_10min_Allper_tract.xlsx','SR_1000ft_LA_10min_Allper_tract.xlsx']
for file in file_name_noise:

    noise_data = pd.read_excel(file)

    column_names = noise_data.columns.values.tolist()

    # Create a GeoDataFrame from the filtered noise data
    geometry = [Point(xy) for xy in zip(noise_data['Longitude'], noise_data['Latitude'])]

    data_dict = {col_name: noise_data[col_name] for col_name in column_names}
    noise_gdf = gpd.GeoDataFrame(data_dict, geometry=geometry, crs='EPSG:4326')

    # Perform a spatial join to combine the filtered census geometries with the filtered noise data
    merged_gdf_inside_bbox = gpd.sjoin(gdf_census_inside_bbox, noise_gdf.set_geometry('geometry'), how='left', op='contains')

    # Plot chloropleth maps for each data type
    for data_column in column_names[2:]:
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
        ax = merged_gdf_inside_bbox.plot(
            column=data_column,  # Use the data column
            cmap='plasma',
            edgecolor='white',
            linewidth=0.15,
            alpha=0.75,
            legend=True,
            legend_kwds={
                "label": r'$L_{AeqT}$' + " (dB)",
                "orientation": "vertical",
                "shrink": 0.75,  # Adjust the color bar size
                "ticks": np.linspace(45, 90, num=6) # Specify fewer bins
            },
            ax=None,
            vmin=45,
            vmax=90
        )

        ax.set_ylim(33.6, 34.4)

        plt.title('')
        plt.ylabel('Latitude $\degree$')
        plt.xlabel('Longitude $\degree$')

        plot_file_name = os.path.join(folder_path, f'{data_column}.png')
        plt.savefig(plot_file_name)
        plt.close()


# file_name= ['combined_results_per_tract.xlsx']
# census_tracts = gpd.read_file('filtered_race_data_lacounty (1).geojson')
# for file in file_name:
#     demo_data = pd.read_excel(file)

#     column_names = demo_data.columns.values.tolist()

#     # Create a GeoDataFrame from the filtered noise data
#     geometry = [Point(xy) for xy in zip(demo_data['Longitude'], demo_data["Latitude"])]

#     # tract_polygon = gdf_census_inside_bbox['geometry'].area

#     # data_dict = {col_name: demo_data[col_name] / tract_polygon if col_name not in ['Longitude', 'Latitude'] else demo_data[col_name] 
#     # for col_name in column_names[3:]}
#     data_dict = {col_name: demo_data[col_name] for col_name in column_names[3:]}
    
#     demo_gdf = gpd.GeoDataFrame(data_dict, geometry=geometry, crs='EPSG:4326')


#     # Perform a spatial join to combine the filtered census geometries with the filtered noise data
#     merged_gdf_inside_bbox = gpd.sjoin(gdf_census_inside_bbox, demo_gdf.set_geometry('geometry'), how='left', op='contains')

#     # Plot chloropleth maps for each data type
#     for data_column in column_names[3:]:
#         plt.figure(figsize=(12, 12))
#         plt.rcParams['axes.linewidth'] = 2.0
#         plt.rcParams["font.family"] = "Times New Roman"
#         parameters = {
#             'axes.labelsize': 16,
#             'xtick.labelsize': 14,
#             'ytick.labelsize': 14,
#             'axes.titlesize': 16,
#             'xtick.major.pad': 1,
#             'ytick.major.pad': 0,
#         }
#         plt.rcParams.update(parameters)

#         ax = merged_gdf_inside_bbox.plot(
#             column=data_column,  # Use the data column
#             cmap='plasma',
#             edgecolor='white',
#             linewidth=0.15,
#             alpha=0.75,
#             legend=True,
#             legend_kwds={
#                 "label": str(data_column) + " (dB)",
#                 "orientation": "vertical",
#                 "shrink": 0.75,  # Adjust the color bar size
#                 # "ticks": np.linspace(45, 90, num=6) # Specify fewer bins
#             },
#             ax=None
#             # vmin=45,
#             # vmax=90
#         )

#         ax.set_ylim(33.6, 34.4)

#         plt.title('')
#         plt.ylabel('Latitude $\degree$')
#         plt.xlabel('Longitude $\degree$')

#         plot_file_name = os.path.join(folder_path, f'{data_column}.png')
#         plt.savefig(plot_file_name)
#         plt.close()
