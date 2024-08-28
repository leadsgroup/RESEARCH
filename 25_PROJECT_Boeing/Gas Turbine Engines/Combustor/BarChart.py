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
    CO2_PSR = 1.58898801
    CO2_PFR = 2.62929774
    CO2_TOT = 3.161527
    CO2_EI  = 3.16
    CO2_err_PSR = ((np.abs(CO2_PSR - CO2_EI))/CO2_EI)*100
    CO2_err_PFR = ((np.abs(CO2_PFR - CO2_EI))/CO2_EI)*100
    CO2_err_TOT = ((np.abs(CO2_TOT - CO2_EI))/CO2_EI)*100
    
    CO_PSR = 0.00893172
    CO_PFR = 0.25463357
    CO_TOT = 0.000061
    CO_EI  = 0.00201
    CO_err_PSR = ((np.abs(CO_PSR - CO_EI))/CO_EI)*100
    CO_err_PFR = ((np.abs(CO_PFR - CO_EI))/CO_EI)*100
    CO_err_TOT = ((np.abs(CO_TOT - CO_EI))/CO_EI)*100
    
    H2O_PSR = 0.62397382
    H2O_PFR = 1.14962815
    H2O_TOT = 1.210699
    H2O_EI  = 1.34 
    H2O_err_PSR = ((np.abs(H2O_PSR - H2O_EI))/H2O_EI)*100
    H2O_err_PFR = ((np.abs(H2O_PFR - H2O_EI))/H2O_EI)*100
    H2O_err_TOT = ((np.abs(H2O_TOT - H2O_EI))/H2O_EI)*100   
    
    NOx_PSR = 0.0003973 + 0.08756719
    NOx_PFR = 6.94051579e-02 + 7.01306185e-05
    NOx_TOT = 0.155204 + 0.202018
    NOx_EI  = 0.01514 
    NOx_err_PSR = ((np.abs(NOx_PSR - NOx_EI))/NOx_EI)*100
    NOx_err_PFR = ((np.abs(NOx_PFR - NOx_EI))/NOx_EI)*100
    NOx_err_TOT = ((np.abs(NOx_TOT - NOx_EI))/NOx_EI)*100    
    
    plot_results(CO2_PSR,CO2_PFR,CO2_TOT,CO2_EI, CO2_err_PSR, CO2_err_PFR, CO2_err_TOT,CO_PSR,CO_PFR,CO_TOT,CO_EI, CO_err_PSR, CO_err_PFR, CO_err_TOT,H2O_PSR,H2O_PFR,H2O_TOT,H2O_EI, H2O_err_PSR, H2O_err_PFR, H2O_err_TOT,NOx_PSR,NOx_PFR,NOx_TOT,NOx_EI, NOx_err_PSR, NOx_err_PFR, NOx_err_TOT)
    return

def plot_results(CO2_PSR,CO2_PFR,CO2_TOT,CO2_EI, CO2_err_PSR, CO2_err_PFR, CO2_err_TOT,CO_PSR,CO_PFR,CO_TOT,CO_EI, CO_err_PSR, CO_err_PFR, CO_err_TOT,H2O_PSR,H2O_PFR,H2O_TOT,H2O_EI, H2O_err_PSR, H2O_err_PFR, H2O_err_TOT,NOx_PSR,NOx_PFR,NOx_TOT,NOx_EI, NOx_err_PSR, NOx_err_PFR, NOx_err_TOT):
        # Get plot style
    ps = plot_style()

    # Data for plotting
    categories = ['PSR', 'PFR', 'PSR + PFR']
    CO2values = [CO2_PSR, CO2_PFR, CO2_TOT]
    CO2errors = [CO2_err_PSR, CO2_err_PFR, CO2_err_TOT]
    COvalues = [CO_PSR, CO_PFR, CO_TOT]
    COerrors = [CO_err_PSR, CO_err_PFR, CO_err_TOT]
    H2Ovalues = [H2O_PSR, H2O_PFR, H2O_TOT]
    H2Oerrors = [H2O_err_PSR, H2O_err_PFR, H2O_err_TOT]
    NOxvalues = [NOx_PSR, NOx_PFR, NOx_TOT]
    NOxerrors = [NOx_err_PSR, NOx_err_PFR, NOx_err_TOT]    

    # Create the figure and axis
    fig, axis_1 = plt.subplots(figsize=(7, 6))
    fig2, axis_2 = plt.subplots(figsize=(7, 6))
    fig3, axis_3 = plt.subplots(figsize=(7, 6))
    fig4, axis_4 = plt.subplots(figsize=(7, 6))
    fig5, axis_5 = plt.subplots(figsize=(7, 6))
    fig6, axis_6 = plt.subplots(figsize=(7, 6))
    fig7, axis_7 = plt.subplots(figsize=(7, 6))
    fig8, axis_8 = plt.subplots(figsize=(7, 6))    
    
    # Generate gradient colors of blue
    cmap = plt.get_cmap("Blues")
    colors = [cmap((i / (len(CO2values) - 0.25)) + 0.5) for i in range(len(CO2values))]    

    # Create bars
    bars = axis_1.bar(categories, CO2values, color=colors)
    bars2 = axis_2.bar(categories, CO2errors, color=colors)
    bars3 = axis_3.bar(categories, COvalues, color=colors)
    bars4 = axis_4.bar(categories, COerrors, color=colors)
    bars5 = axis_5.bar(categories, H2Ovalues, color=colors)
    bars6 = axis_6.bar(categories, H2Oerrors, color=colors)
    bars7 = axis_7.bar(categories, NOxvalues, color=colors)
    bars8 = axis_8.bar(categories, NOxerrors, color=colors)

    # Add a horizontal dotted line for CO2 EI
    axis_1.axhline(y=CO2_EI, color='r', linestyle='--', label=f'EI CO2 = {CO2_EI:.2f}')
    axis_3.axhline(y=CO_EI, color='r', linestyle='--', label=f'EI CO = {CO_EI:.2f}')
    axis_5.axhline(y=H2O_EI, color='r', linestyle='--', label=f'EI H2O = {H2O_EI:.2f}')
    axis_7.axhline(y=NOx_EI, color='r', linestyle='--', label=f'EI NOx = {NOx_EI:.2f}')

    # Set labels and title
    axis_1.set_ylabel('EI [kg/kg]')
    axis_1.set_title('Emission Index of CO2')
    axis_2.set_ylabel('Percentage error [%]')
    axis_2.set_title('Percentage error for CO2 EI models')
    axis_3.set_ylabel('EI [kg/kg]')
    axis_3.set_title('Emission Index of CO')
    axis_4.set_ylabel('Percentage error [%]')
    axis_4.set_title('Percentage error for CO EI models')
    axis_5.set_ylabel('EI [kg/kg]')
    axis_5.set_title('Emission Index of H2O')
    axis_6.set_ylabel('Percentage error [%]')
    axis_6.set_title('Percentage error for H2O EI models')
    axis_7.set_ylabel('EI [kg/kg]')
    axis_7.set_title('Emission Index of NOx')
    axis_8.set_ylabel('Percentage error [%]')
    axis_8.set_title('Percentage error for NOx EI models')    

    # Add legend
    axis_1.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.05, 0.95),loc='upper left')
    axis_3.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.05, 0.95),loc='upper left') 
    axis_5.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.05, 0.95),loc='upper left') 
    axis_7.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.05, 0.95),loc='upper left') 

    # Adjust layout
    fig.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()
    fig6.tight_layout()
    fig7.tight_layout()
    fig8.tight_layout()
    

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