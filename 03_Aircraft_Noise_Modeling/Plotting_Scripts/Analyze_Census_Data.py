# Python imports
import matplotlib.pyplot as plt  
import numpy as np 
import geopandas as gpd
from scipy.stats import gaussian_kde
from shapely.geometry import Point
from matplotlib.colors import Normalize

# Load combined race data
combined_gdf_race = gpd.read_file(
    '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles/data counties combined/3combined_data_race.geojson'
)

# Create 1x3 subplots for race-related plots
fig, axs = plt.subplots(1, 3, figsize=(18, 6))
ax1, ax2, ax3 = axs

# Define the color scale range
vmin, vmax = 0, 100

# Plot 1: Percent Hispanic or Latino
combined_gdf_race.plot(
    column='Percent Hispanic or Latino',
    cmap='RdYlGn_r',
    edgecolor='white',
    linewidth=0.0,
    alpha=0.75,
    ax=ax1,
    legend=True,
    vmin=vmin,
    vmax=vmax
)
ax1.axis("off")
ax1.set_title("Percent Hispanic or Latino")

# Plot 2: Percent Non-Hispanic Black
combined_gdf_race.plot(
    column="Percent Non Hispanic Black or African American alone",
    cmap='RdYlGn_r',
    edgecolor='white',
    linewidth=0.0,
    alpha=0.75,
    ax=ax2,
    legend=True,
    vmin=vmin,
    vmax=vmax
)
ax2.axis("off")
ax2.set_title("Percent Non Hispanic Black")

# Plot 3: Percent Non-Hispanic White
combined_gdf_race.plot(
    column="Percent Non Hispanic White alone",
    cmap='RdYlGn_r',
    edgecolor='white',
    linewidth=0.0,
    alpha=0.75,
    ax=ax3,
    legend=True,
    vmin=vmin,
    vmax=vmax
)
ax3.axis("off")
ax3.set_title("Percent Non Hispanic White")

plt.tight_layout()
plt.show()

# Load combined income data
combined_gdf_income = gpd.read_file(
    '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles/data counties combined/3combined_data_income.geojson'
)

# Columns to sum for income greater than $30,000
columns_to_sum2 = [
    '$30,000 to $34,999', '$35,000 to $39,999', '$40,000 to $44,999', '$45,000 to $49,999',
    '$50,000 to $59,999', '$60,000 to $74,999', '$75,000 to $99,999', '$100,000 to $124,999',
    '$125,000 to $149,999', '$150,000 to $199,999', '$200,000 or more'
]

# Sum the specified columns and create a new column called 'TotH'
combined_gdf_income['TotH'] = combined_gdf_income[columns_to_sum2].sum(axis=1)

# Create 1x2 subplots for income-related plots
fig, axs = plt.subplots(1, 2, figsize=(15, 6))
ax1, ax2 = axs

# Plot 1: Population with income less than $30,000
combined_gdf_income.plot(
    column='TotL',
    cmap='plasma',
    edgecolor='white',
    linewidth=0.0,
    alpha=0.75,
    ax=ax1,
    legend=True
)
ax1.axis("off")
ax1.set_title("Population with Total Income Less Than $30,000")

# Plot 2: Population with income $30,000 or more
combined_gdf_income.plot(
    column="TotH",
    cmap='plasma',
    edgecolor='white',
    linewidth=0.0,
    alpha=0.75,
    ax=ax2,
    legend=True
)
ax2.axis("off")
ax2.set_title("Population with Total Income $30,000 or More")

plt.tight_layout()
plt.show()

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
    boundary_mask = boundary_gdf.geometry.unary_union
    mask = np.array([[boundary_mask.contains(Point(X[i, j], Y[i, j])) for j in range(Z.shape[1])] for i in range(Z.shape[0])])
    Z_masked = np.ma.masked_where(~mask, Z)

    # Plot LA County boundary and KDE
    boundary_gdf.plot(ax=ax, color='lightgrey', edgecolor='black')
    contour = ax.contourf(X, Y, Z_masked, levels=20, cmap='viridis', alpha=0.7, norm=Normalize(vmin=Z.min(), vmax=Z.max()))
    plt.colorbar(contour, ax=ax, label='Density')
    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect('equal')

# Main function for KDE plots
def main():
    # Load LA County boundary data
    boundary_gdf = gpd.read_file(
        '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles/data race la county/Edited_CensusData_LACounty_filtered.geojson'
    )

    # Sensitive area datasets
    data_files = {
        'Schools': '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles/various other data/Schools_Colleges_and_Universities_-4260208272874444542.geojson',
        'Hospitals': '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles/various other data/Hospitals_and_Medical_Centers.geojson',
        'Churches': '/Users/avacipriani/Desktop/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/Cencus_Data/Los_Angeles/various other data/Churches.geojson'
    }

    # Plot KDE for each sensitive area
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    for i, (title, file) in enumerate(data_files.items()):
        sensitive_gdf = gpd.read_file(file)
        coords = np.vstack([sensitive_gdf.geometry.x, sensitive_gdf.geometry.y])
        kde_density_plot(axs[i], coords, boundary_gdf, f'LA County {title} Density Contour')

    plt.tight_layout()
    plt.show()

    # Plot KDE for all combined sensitive areas
    fig, ax = plt.subplots(figsize=(12, 8))
    all_coords = np.hstack([np.vstack([gpd.read_file(file).geometry.x, gpd.read_file(file).geometry.y]) for file in data_files.values()])
    kde_density_plot(ax, all_coords, boundary_gdf, 'LA County Combined Density Contour for Schools, Churches, and Hospitals/Medical Centers')
    plt.show()

# Run the main function
if __name__ == '__main__':
    main()
