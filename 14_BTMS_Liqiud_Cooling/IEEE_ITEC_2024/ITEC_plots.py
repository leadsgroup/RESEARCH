
# Python Imports  
from RCAIDE.Core import Data
import numpy as np    
import matplotlib.pyplot as plt  
import matplotlib.cm as cm  
 
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  

    save_figure   = True
    file_type     = '.png'
    fig_name      = 'Figure_Name' 
    
    # define ploting parameters 
    PP           = define_plotting_parameters(data_range_length=5)  
        
    # Given data for the first set
    mass_flow_rate = [0.1, 0.5, 1, 1.5, 2]
    channel_width = [0.001068809, 0.001, 0.001, 0.001, 0.002990134]
    angle_of_contact_radians = [0.823273618, 0.829030395, 0.378364017, 0.27953316, 0.772574583]
    angle_of_contact_degrees = np.degrees(angle_of_contact_radians)
    
    # Plot Channel Width and Angle vs Mass Flow Rate
    fig_name = 'Figure_1'
    fig = plt.figure(fig_name)    
    fig.set_size_inches(PP.figure_width,PP.figure_height)  
    ax1 = fig.add_subplot(1,1,1)   
    
    # Plot Channel Width on the left Y-axis
    ax1.set_xlabel('Mass Flow Rate of Coolant (kg/s)')
    ax1.set_ylabel('Channel Width (m)', color='tab:blue')
    ax1.plot(mass_flow_rate, channel_width, color='tab:blue', marker='o', linewidth=PP.line_width, label='Channel Width')  # Added label
    ax1.tick_params(axis='y', labelcolor='tab:blue', direction='in')  # Added direction parameter
    
    # Create a secondary axis for Angle of Contact on the right Y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Angle of Contact (degrees)', color='black')
    ax2.plot(mass_flow_rate, angle_of_contact_degrees, color='black', linestyle='--', marker='s', linewidth=PP.line_width, label='Angle of Contact')  # Added label
    ax2.tick_params(axis='y', labelcolor='black', direction='in')  # Added direction parameter
    
    # Display the combined legend on the inside upper left
    #ax1.legend(loc='upper center')  # Adjusted location
    ax1.grid(True, linestyle='--', alpha=0.7)  # Added grid lines

    plt.tight_layout()
    
    # save results 
    if save_figure: 
        plt.savefig( fig_name+ file_type)    
        
    
    # Given data for the second set
    mass = [42.64409364, 42.75128108, 20.25481331, 15.32135184, 44.998556]
    power = [0.013889811, 0.320564344, 0.568661325, 0.951951201, 1.645415373]
    
    # Plot Power and Mass vs Mass Flow Rate 
    fig_name = 'Figure_2'
    fig = plt.figure(fig_name)    
    fig.set_size_inches(PP.figure_width,PP.figure_height)  
    ax1 = fig.add_subplot(1,1,1)    
    
    # Plot Mass on the left Y-axis
    ax1.set_xlabel('Mass Flow Rate of Coolant (kg/s)')
    ax1.set_ylabel('Mass (kg)', color='tab:green')
    ax1.plot(mass_flow_rate, mass, color='tab:green', marker='^', linewidth=PP.line_width, label='Mass')  # Added label
    ax1.tick_params(axis='y', labelcolor='tab:green', direction='in')  # Added direction parameter
    
    # Create a secondary axis for Power on the right Y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Power (W)', color='tab:orange')
    ax2.plot(mass_flow_rate, power, color='tab:orange', linestyle='--', marker='*', linewidth=PP.line_width, label='Power')  # Added label
    ax2.tick_params(axis='y', labelcolor='tab:orange', direction='in')  # Added direction parameter
    
    # Display the combined legend on the inside upper left
    #ax1.legend(loc='upper center')  # Adjusted location
    ax1.grid(True, linestyle='--', alpha=0.7)  # Added grid lines

    plt.tight_layout()

    # save results 
    if save_figure: 
        plt.savefig( fig_name+ file_type)       
    
    # Given data
    mass_flow_rate = [0.11, 0.5, 1, 1.5, 2, 2.5]
    height_of_stack = [0.010367271, 0.026737309, 0.040017113, 0.048559426, 0.053602384, 0.053677281]
    length_of_hot_side = [1.86420135, 1.647889512, 1.713028869, 1.924355652, 2.376950397, 3.452468908]
    length_of_cold_side = [1.378999049, 0.490059006, 0.335372408, 0.298018125, 0.30980726, 0.394728004]
    
    # Create a figure and axis 
    fig_name = 'Figure_3'
    fig = plt.figure(fig_name)    
    fig.set_size_inches(PP.figure_width,PP.figure_height)  
    ax1 = fig.add_subplot(1,1,1)     
    
    # Plot mass flow rate and height_of_stack on the left Y-axis
    ax1.set_xlabel('Mass Flow Rate (kg/s)')
    ax1.set_ylabel('Height of Stack (m)', color='tab:blue')
    ax1.plot(mass_flow_rate, height_of_stack, color='tab:blue', marker='o', linewidth=2.5, label='Height of Stack')  # Added label
    ax1.tick_params(axis='y', labelcolor='tab:blue', direction='in')  # Added direction parameter
    
    # Create a secondary axis for length_of_hot_side and length_of_cold_side on the right Y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Length (m)', color='black')  # Changed color to black
    ax2.plot(mass_flow_rate, length_of_hot_side, color='black', linestyle='--', marker='s', linewidth=2.5, label='Length of Hot Side')  # Added label, changed color
    ax2.plot(mass_flow_rate, length_of_cold_side, color='black', linestyle='-', marker='^', linewidth=2.5, label='Length of Cold Side')  # Added label, changed color
    ax2.tick_params(axis='y', labelcolor='black', direction='in')  # Added direction parameter
    
    # Display the combined legend on the inside upper left with decreased size
    ax2.legend(loc='upper left', fontsize='x-small')  # Adjusted location and decreased size
    ax1.grid(True, linestyle='--', alpha=0.7)  # Added grid lines

    plt.tight_layout()

    # save results 
    if save_figure: 
        plt.savefig( fig_name+ file_type)       
    
    
    mass_flow_rate_air = [0.11, 0.5, 1, 1.5, 2, 2.5]
    mass = [61.84426288, 50.10399404, 53.34773164, 64.6219048, 91.59547184, 169.7447745]
    power_watt = [47651.40283, 12827.62044, 10297.28773, 11069.56922, 12703.43629, 14682.51306]
    
    # Convert power to kilowatts
    power_kw = [p / 1000 for p in power_watt]
    
    # Create a figure and axis
    fig_name = 'Figure_4'
    fig = plt.figure(fig_name)    
    fig.set_size_inches(PP.figure_width,PP.figure_height)  
    ax1 = fig.add_subplot(1,1,1)    
    
    # Plot mass flow rate of air and mass on the left Y-axis
    ax1.set_xlabel('Mass Flow Rate of Air (kg/s)')
    ax1.set_ylabel('Mass (kg)', color='tab:green')
    ax1.plot(mass_flow_rate_air, mass, color='tab:green', marker='o', linewidth=2.5, label='Mass')  # Added label
    ax1.tick_params(axis='y', labelcolor='tab:green', direction='in')  # Added direction parameter
    
    # Create a secondary axis for power_kw on the right Y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Power (kW)', color='tab:orange')
    ax2.plot(mass_flow_rate_air, power_kw, color='tab:orange', linestyle='--', marker='s', linewidth=2.5, label='Power (kW)')  # Added label
    ax2.tick_params(axis='y', labelcolor='tab:orange', direction='in')  # Added direction parameter
    
    # Display the combined legend on the inside upper left
    #ax1.legend(loc='upper center')  # Adjusted location
    ax1.grid(True, linestyle='--', alpha=0.7)  # Added grid lines
    
    plt.tight_layout()

    # save results 
    if save_figure: 
        plt.savefig( fig_name+ file_type)       
    return  

    
def define_plotting_parameters(data_range_length):
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 12,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 2 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 3.5
    plot_parameters.figure_height    = 3
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 10 
    plot_parameters.plot_grid        = True   
    plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    plot_parameters.colors           = cm.inferno(np.linspace(0,1,data_range_length))                
    plot_parameters.legend_font      = 10                             # legend_font_size          
    
    return plot_parameters
  
if __name__ == '__main__': 
    main()    
    plt.show()   
 
 