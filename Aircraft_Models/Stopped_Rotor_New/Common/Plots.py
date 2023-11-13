'''
# Plots.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------  

from RCAIDE.Visualization import *  

# ----------------------------------------------------------------------
#   Plot Results
# ----------------------------------------------------------------------  
def plot_results(results,run_noise_model,save_figure_flag):  
    
    # Plots fligh conditions 
    plot_flight_conditions(results) 
    
    # Plot arcraft trajectory
    plot_flight_trajectory(results)
    
    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results)  
     
    # Plot Aircraft Stability
    plot_stability_coefficients(results) 
    
    # Plot Aircraft Electronics
    plot_battery_pack_conditions(results) 
    plot_battery_health_conditions(results)
    plot_battery_cell_conditions(results)
    plot_battery_C_rate(results)
    plot_battery_degradation(results) 
    
    # Plot Propeller Conditions 
    plot_rotor_conditions(results) 
    plot_disc_and_power_loading(results)
    
    # Plot Electric Motor and Propeller Efficiencies 
    plot_electric_motor_and_rotor_efficiencies(results) 
    
    # Plot Battery Degradation  
    plot_battery_degradation(results)   
    
    #if run_noise_model: 
        ## Plot noise level
        #plot_noise_level(results)
        
        ## Plot noise contour
        #plot_2D_noise_contour(results) 
                        
    return
