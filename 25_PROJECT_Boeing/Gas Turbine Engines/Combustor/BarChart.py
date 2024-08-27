# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units , Data   
from RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor          import design_turbofan
from RCAIDE.Framework.Mission.Common      import  Conditions

# python imports 
import numpy as np  
import pickle
from copy import deepcopy
import matplotlib.pyplot as plt  
import os   
import matplotlib.cm as cm

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(): 
    # Define Engine  
    CO2_PSR = 2.5
    CO2_PFR = 2.6
    CO2_TOT = 3.16
    CO2_EI  = 3.16
    plot_results(CO2_PSR,CO2_PFR,CO2_TOT,CO2_EI)
    return

def plot_results(CO2_PSR,CO2_PFR,CO2_TOT,CO2_EI):
        # Get plot style
    ps = plot_style()

    # Data for plotting
    categories = ['PSR', 'PFR', 'PSR + PFR']
    values = [CO2_PSR, CO2_PFR, CO2_TOT]

    # Create the figure and axis
    fig, axis_1 = plt.subplots(figsize=(7, 6))

    # Create bars
    bars = axis_1.bar(categories, values, color=ps.color[0:3])

    # Add a horizontal dotted line for CO2 EI
    axis_1.axhline(y=CO2_EI, color='r', linestyle='--', label=f'EI CO2 = {CO2_EI:.2f}')

    # Set labels and title
    axis_1.set_ylabel('EI [kg/kg]')
    axis_1.set_title('Emission Index of CO2')

    # Add legend
    axis_1.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.05, 0.95),loc='upper left')     

    # Adjust layout
    fig.tight_layout()

    # Show plot
    plt.show()
    
    return

def plot_style(number_of_lines= 10): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  #figure.dpi': 1200
                  }


    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1.5  
    plot_parameters.line_style             = ['-','--']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '12'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size         = 16
    plot_parameters.title_font_size        = 16   
    plot_parameters.markers                =  ['o','x','o','v','P','p','^','D','*']
    plot_parameters.color                  = cm.inferno(np.linspace(0,0.9,number_of_lines)) 


    return plot_parameters


# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename): 
    pickle_file  =  filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return   

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results  



if __name__ == '__main__': 
    main()
    plt.show()