# Python imports
import matplotlib.pyplot as plt  
import sys 
import numpy as np 
import json 
import geopandas as gpd
from scipy.stats import gaussian_kde
from shapely.geometry import Point
from matplotlib.colors import Normalize

# Function to calculate KDE and plot the contour map
def kde_density_plot(ax, coords, boundary_gdf, title):
    kde = gaussian_kde(coords)
    xmin, ymin, xmax, ymax = boundary_gdf.total_bounds
    x = np.linspace(xmin, xmax, 100)
    y = np.linspace(ymin, ymax, 100)
    X, Y = np.meshgrid(x, y)
    xy = np.vstack([X.ravel(), Y.ravel()])
    
    Z = np.reshape(kde(xy).T, X.shape)

    # Create a mask for grid points inside the boundary
    mask = np.zeros(Z.shape, dtype=bool)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            point = Point(X[i, j], Y[i, j])
            if boundary_gdf.geometry.contains(point).any():
                mask[i, j] = True
    Z_masked = np.ma.masked_where(~mask, Z)

    # Plot LA County boundary and KDE
    boundary_gdf.plot(ax=ax, color='lightgrey', edgecolor='black')
    contour = ax.contourf(X, Y, Z_masked, levels=20, cmap='viridis', alpha=0.7, 
                          norm=Normalize(vmin=Z.min(), vmax=Z.max()))
    plt.colorbar(contour, ax=ax, label='Density')
    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect('equal')

# Main function
def main(): 
# replace username in file paths with your username
    boundary_gdf = gpd.read_file('/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles_Metroplex/data race la county/Edited_CensusData_LACounty_filtered.geojson')

    # Sensitive area datasets
    data_files = {
        'Schools': '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles_Metroplex/various other data/Schools_Colleges_and_Universities_-4260208272874444542.geojson',
        'Hospitals': '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles_Metroplex/various other data/Hospitals_and_Medical_Centers.geojson',
        'Churches': '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles_Metroplex/various other data/Churches.geojson'
    }

    # Plot KDE for each sensitive area
    fig, axs = plt.subplots(1, 3, figsize=(18, 12))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, wspace=0.5) 
    for i, (title, file) in enumerate(data_files.items()):
        sensitive_gdf = gpd.read_file(file)
        coords = np.vstack([sensitive_gdf.geometry.x, sensitive_gdf.geometry.y])
        kde_density_plot(axs[i], coords, boundary_gdf, f'LA County {title} Density Contour')

    plt.show()

    # Plot KDE for all combined sensitive areas
    fig, ax = plt.subplots(figsize=(12, 12))
    all_coords = np.hstack([np.vstack([gpd.read_file(file).geometry.x, gpd.read_file(file).geometry.y]) 
                            for file in data_files.values()])
    kde_density_plot(ax, all_coords, boundary_gdf, 'LA County Combined Density Contour for Schools, Churches, and Hospitals/Medical Centers')
    plt.show()

if __name__ == '__main__': 
    main()
