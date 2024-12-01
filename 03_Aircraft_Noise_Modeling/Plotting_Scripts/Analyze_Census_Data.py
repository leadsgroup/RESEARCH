# Python imports
from RCAIDE.Framework.Core import  Data
import matplotlib.pyplot as plt  
import numpy as np 
import geopandas as gpd
from scipy.stats import gaussian_kde
from shapely.geometry import Point
from matplotlib.colors import Normalize
import seaborn as sns

import os  
import sys 
import pandas as pd
import numpy as  np 
import matplotlib.cm as cm

local_path_1 =  os.path.split(sys.path[0])[0] 
sys.path.append( os.path.join(local_path_1, 'Papers_and_Presentations'))


def main(): 
    #generate_population_plots()
    
    generate_census_analysis_plots()
    
    return  

def generate_population_plots():

    ospath                = os.path.abspath(__file__)
    separator             = os.path.sep
    relative_path         = os.path.dirname(ospath) + separator + '..' + 'Cencus_Data'  + separator + 'Los_Angeles' + separator  
    
    # Load LA County boundary data
    boundary_gdf = gpd.read_file(
        relative_path + 'data race la county'  + separator + 'Edited_CensusData_LACounty_filtered.geojson'
    )

    # Sensitive area datasets
    data_files = {
        'Schools': relative_path + 'various other data'  + separator + 'Schools_Colleges_and_Universities_-4260208272874444542.geojson',
        'Hospitals': relative_path + 'various other data'  + separator + 'Hospitals_and_Medical_Centers.geojson',
        'Churches': relative_path +  'various other data'  + separator + 'Churches.geojson'
    }

    # Plot KDE for each sensitive area
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    for i, (title, file) in enumerate(data_files.items()):
        sensitive_gdf = gpd.read_file(file)
        coords = np.vstack([sensitive_gdf.geometry.x, sensitive_gdf.geometry.y])
        kde_density_plot(axs[i], coords, boundary_gdf, f'LA County {title} Density Contour')

    plt.tight_layout()
    
    # Plot KDE for all combined sensitive areas
    fig, ax = plt.subplots(figsize=(12, 8))
    all_coords = np.hstack([np.vstack([gpd.read_file(file).geometry.x, gpd.read_file(file).geometry.y]) for file in data_files.values()])
    kde_density_plot(ax, all_coords, boundary_gdf, 'LA County Combined Density Contour for Schools, Churches, and Hospitals/Medical Centers')
    
        
    # Load combined race data
    combined_gdf_race = gpd.read_file(
        relative_path + 'data counties combined'  + separator + '3combined_data_race.geojson'
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
        relative_path + 'data counties combined'  + separator + '3combined_data_income.geojson'
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
    return 


def generate_census_analysis_plots():

    # ----------------------------------------------------------------------------------------------------------------------
    # FILE IMPORTS 
    # ----------------------------------------------------------------------------------------------------------------------            
    ospath                = os.path.abspath(__file__)
    separator             = os.path.sep 
    routes_filepath       = sys.path[-1] + separator  + 'Scitech_2025' + separator  + 'LA_Census_Data.csv' 
    cencus_data           = pd.read_csv(routes_filepath)
    
    income_ranges = ['<50', '50-100', '100-150', '150-200', '>200']
    income_tag    =  ['Perc_50k', 'Perc_100k','Perc_150k','Perc_200k','Perc_200kplus']
    pp = define_plot_parameters()
    #colors = cm.tab10(6)
    #colors = cm.tab20(6)
    #colors = cm.tab20c(6)
    for i in  range(len(income_ranges)):
        fig_name = 'income_range_' + str(i+1)  
        plt.figure(figsize=(4, 6))
        sns.regplot(data=cencus_data, x="Perc_Asian"          , scatter=False, y= income_tag[i], line_kws={"color": pp.colors[0]}, label = 'Asian') 
        sns.regplot(data=cencus_data, x="Perc_Black"          , scatter=False, y= income_tag[i], line_kws={"color": pp.colors[1]}, label = 'Black') 
        sns.regplot(data=cencus_data, x="Perc_Hispanic"       , scatter=False, y= income_tag[i], line_kws={"color": pp.colors[2]}, label = 'Hispanic') 
        sns.regplot(data=cencus_data, x="Perc_AmericanIndian" , scatter=False, y= income_tag[i], line_kws={"color": pp.colors[3]}, label = 'Native American') 
        sns.regplot(data=cencus_data, x="Perc_PacificIslander", scatter=False, y= income_tag[i], line_kws={"color": pp.colors[4]}, label = 'Pacific Islander')
        sns.regplot(data=cencus_data, x="Perc_White"          , scatter=False, y= income_tag[i], line_kws={"color": pp.colors[5]}, label = 'White' )  
        plt.ylabel('% Population Income Between ' + income_ranges[i] + 'k per Tract') 
        plt.xlabel('% Population of Given Demographic') 
        plt.tight_layout()
        
    plt.figure(figsize=(10, 1))
    plt.legend(fontsize = pp.legend_font_size,bbox_to_anchor=(0.05, 0.95),loc='upper center')          

    
    return 

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


# ------------------------------------------------------------------ 
# Define plot parameters 
# ------------------------------------------------------------------  
def define_plot_parameters(): 

    plt.rcParams['axes.linewidth'] = 2.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 12,
                  'legend.fontsize': 12,
                  'xtick.labelsize': 12,
                  'ytick.labelsize': 12,
                  'axes.titlesize': 12}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 2
    plot_parameters.line_styles      = ['-','--','--',':','--']
    plot_parameters.figure_width     = 10
    plot_parameters.figure_height    = 7
    plot_parameters.marker_size      = 10
    plot_parameters.legend_font_size = 20
    plot_parameters.alpha_val        = 0.25
    plot_parameters.root_color       = 'grey'
    plot_parameters.plot_grid        = True   

    plot_parameters.colors           = cm.tab10(np.linspace(0, 1, 10)) 

    plot_parameters.colors_1         = ['black','darkmagenta','mediumblue','darkgreen','darkgoldenrod','darkred']   
    plot_parameters.colors_2         = ['grey','orchid','darkcyan','green','orange','red']        
    plot_parameters.markers          = ['s','o','v','P','p','^','D','X','*']   
    
    return plot_parameters

# Run the main function
if __name__ == '__main__':
    main()
    plt.show()