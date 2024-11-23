import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
import pandas as pd

# Styling function
def plot_style():
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {
        'axes.labelsize': 6,
        'xtick.labelsize': 6,
        'ytick.labelsize': 3,
        'figure.dpi': 600
    }
    plt.rcParams.update(parameters)

# Apply styling
plot_style()

# Get Set2 colormap
cmap = get_cmap("Set2")
colors = [cmap(i) for i in range(5)]

def emissionsCO2_plot(title, time, CO2, color, filename):
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot CO2 emissions for each fuel type
    labels = ["Jet A1", "Ethane", "Methane", "Propane", "LNG"]
    for i in range(len(CO2)):
        ax.plot(time, CO2[i], label=labels[i], color=color[i], linewidth=1)

    # Labels and ticks
    ax.set_ylabel("CO2 Emissions (CO2e kg)", fontsize=5.5, labelpad=0)
    ax.set_xlabel("Time (minutes)", fontsize=5.5, labelpad=0)
    
    # Set more ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(50))  # Major ticks every 50 seconds
    ax.tick_params(axis='x', which='major', labelsize=4.5)  # Font size for major x-axis ticks
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))   # Minor ticks every 10 seconds
    
    ax.yaxis.set_major_locator(plt.MultipleLocator(1e4))  # Major ticks every 10,000
    ax.tick_params(axis='y', which='major', labelsize=4.5)  # Font size for major y-axis ticks
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.2e4))  # Minor ticks every 2000
    
    # Major and minor gridlines
    ax.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.7)
    ax.grid(True, which='minor', linestyle='-', linewidth=0.3, color='gray', alpha=0.3)
    
    # Ensure gridlines are behind the lines
    ax.set_axisbelow(True)
    
    # Add legend
    ax.legend(fontsize=3, loc='upper left')
    
    # Adjust layout
    plt.subplots_adjust(left=0.15, right=0.85, top=0.9, bottom=0.15)

    # Save and show the plot
    plt.savefig(filename, dpi=600)
    plt.show()

def emissionsH2O_plot(title, time, H20, color, filename):
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot H20 emissions for each fuel type
    labels = ["Jet A1", "Ethane", "Methane", "Propane", "LNG"]
    for i in range(len(H20)):
        ax.plot(time, H20[i], label=labels[i], color=color[i], linewidth=1)

    # Labels and ticks
    ax.set_ylabel("H20 Emissions (CO2e kg)", fontsize=5.5, labelpad=0)
    ax.set_xlabel("Time (minutes)", fontsize=5.5, labelpad=0)
    
    # Set more ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(50))  # Major ticks every 50 seconds
    ax.tick_params(axis='x', which='major', labelsize=4.5)  # Font size for major x-axis ticks
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))   # Minor ticks every 10 seconds
    
    ax.yaxis.set_major_locator(plt.MultipleLocator(5e2))  # Major ticks every 500
    ax.tick_params(axis='y', which='major', labelsize=4.5)  # Font size for major y-axis ticks
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1e2))  # Minor ticks every 100
    
    # Major and minor gridlines
    ax.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.7)
    ax.grid(True, which='minor', linestyle='-', linewidth=0.3, color='gray', alpha=0.3)
    
    # Ensure gridlines are behind the lines
    ax.set_axisbelow(True)
    
    # Add legend
    ax.legend(fontsize=3, loc='upper left')
    
    # Adjust layout
    plt.subplots_adjust(left=0.15, right=0.85, top=0.9, bottom=0.15)

    # Save and show the plot
    plt.savefig(filename, dpi=600)
    plt.show()

def emissionIndexCO2_plot(title, time, EICO2, color, filename):
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot CO2 emissions for each fuel type
    labels = ["Jet A1", "Ethane", "Methane", "Propane", "LNG"]
    for i in range(len(EICO2)):
        ax.plot(time, EICO2[i], label=labels[i], color=color[i], linewidth=1)

    # Set labels 
    ax.set_ylabel(r"$\mathrm{CO}_2\ \text{Emission Index}\ \left(\frac{kg}{kg_{\mathrm{fuel}}}\right)$", fontsize=5.5, labelpad=0)
    ax.set_xlabel("Time (minutes)", fontsize=5.5, labelpad=0)
    
    # Set ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(50))  # Major ticks every 50 seconds
    ax.tick_params(axis='x', which='major', labelsize=4.5)  # Font size for major x-axis ticks
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))   # Minor ticks every 10 seconds
    
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.01))  # Major ticks every 0.01
    ax.tick_params(axis='y', which='major', labelsize=4.5)  # Font size for major y-axis ticks
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.002))  # Minor ticks every 0.002
    
    # Major and minor gridlines
    ax.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.7)
    ax.grid(True, which='minor', linestyle='-', linewidth=0.3, color='gray', alpha=0.3)
    
    # Ensure gridlines are behind the lines
    ax.set_axisbelow(True)
    
    # Add legend
    ax.legend(fontsize=3, loc='upper left')
    
    # Adjust layout
    plt.subplots_adjust(left=0.15, right=0.85, top=0.9, bottom=0.15)

    # Save and show the plot
    plt.savefig(filename, dpi=600)
    plt.show()

def emissionIndexH2O_plot(title, time, EIH2O, color, filename):
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot CO2 emissions for each fuel type
    labels = ["Jet A1", "Ethane", "Methane", "Propane", "LNG"]
    for i in range(len(EIH2O)):
        ax.plot(time, EIH2O[i], label=labels[i], color=color[i], linewidth=1)

    # Set labels 
    ax.set_ylabel(r"$\mathrm{H}_2O\ \text{Emission Index}\ \left(\frac{kg}{kg_{\mathrm{fuel}}}\right)$", fontsize=5.5, labelpad=0)
    ax.set_xlabel("Time (minutes)", fontsize=5.5, labelpad=0)
    
    # Set ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(50))  # Major ticks every 50 seconds
    ax.tick_params(axis='x', which='major', labelsize=4.5)  # Font size for major x-axis ticks
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))   # Minor ticks every 10 seconds
    
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.01))  # Major ticks every 0.01
    ax.tick_params(axis='y', which='major', labelsize=4.5)  # Font size for major y-axis ticks
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.002))  # Minor ticks every 0.002
    
    # Major and minor gridlines
    ax.grid(True, which='major', linestyle='-', linewidth=0.5, color='gray', alpha=0.7)
    ax.grid(True, which='minor', linestyle='-', linewidth=0.3, color='gray', alpha=0.3)
    
    # Ensure gridlines are behind the lines
    ax.set_axisbelow(True)
    
    # Add legend
    ax.legend(fontsize=3, loc='upper left')
    
    # Adjust layout
    plt.subplots_adjust(left=0.15, right=0.85, top=0.9, bottom=0.15)

    # Save and show the plot
    plt.savefig(filename, dpi=600)
    plt.show()


# Data
df1 = pd.read_csv('Emissions/Jet_A1_EmissionsB737.csv')

time = df1['Time'].to_numpy()
CO2_JetA1 = df1['emissions_CO2'].to_numpy()
H2O_JetA1 = df1['emissions_H2O'].to_numpy()
EICO2_JetA1 = df1['emission_index_CO2'].to_numpy()
EIH2O_JetA1 = df1['emission_index_H2O'].to_numpy()


df2 = pd.read_csv('Emissions/Ethane_EmissionsB737.csv')

CO2_Ethane = df2['emissions_CO2'].to_numpy()
H2O_Ethane = df2['emissions_H2O'].to_numpy()
EICO2_Ethane = df2['emission_index_CO2'].to_numpy()
EIH2O_Ethane = df2['emission_index_H2O'].to_numpy()

df3 = pd.read_csv('Emissions/Methane_EmissionsB737.csv')

CO2_Methane = df3['emissions_CO2'].to_numpy()
H2O_Methane = df3['emissions_H2O'].to_numpy()
EICO2_Methane = df3['emission_index_CO2'].to_numpy()
EIH2O_Methane = df3['emission_index_H2O'].to_numpy()

df4 = pd.read_csv('Emissions/Propane_EmissionsB737.csv')

CO2_Propane = df4['emissions_CO2'].to_numpy()
H2O_Propane = df4['emissions_H2O'].to_numpy()
EICO2_Propane = df4['emission_index_CO2'].to_numpy()
EIH2O_Propane = df4['emission_index_H2O'].to_numpy()

df5 = pd.read_csv('Emissions/Liquid_Natural_Gas_EmissionsB737.csv')

CO2_LNG = df5['emissions_CO2'].to_numpy()
H2O_LNG = df5['emissions_H2O'].to_numpy()
EICO2_LNG = df5['emission_index_CO2'].to_numpy()
EIH2O_LNG = df5['emission_index_H2O'].to_numpy()


CO2emissions = [CO2_JetA1, CO2_Ethane, CO2_Methane, CO2_Propane, CO2_LNG]
H20emissions = [H2O_JetA1, H2O_Ethane, H2O_Methane, H2O_Propane, H2O_LNG]
EICO2 = [EICO2_JetA1, EICO2_Ethane, EICO2_Methane, EICO2_Propane, EICO2_LNG]
EIH20 = [EIH2O_JetA1, EIH2O_Ethane, EIH2O_Methane, EIH2O_Propane, EIH2O_LNG]

EICO2_ = [EICO2_JetA1]

# Create individual plots
emissionsCO2_plot(" B737 CO2 Emissions", time, CO2emissions,  colors, "B737_CO2Emissions.png")
emissionsH2O_plot(" B737 H2O Emissions", time, H20emissions,  colors, "B737_H2OEmissions.png")
emissionIndexCO2_plot(" B737 CO2 Emission Index", time, EICO2_,  colors, "B737_CO2EmissionIndex.png")

# emissionIndexCO2_plot(" B737 CO2 Emission Index", time, EICO2,  colors, "B737_CO2EmissionIndex.png")
# emissionIndexH2O_plot(" B737 H2O Emission Index", time, EIH20,  colors, "B737_H2OEmissionIndex.png")
                      