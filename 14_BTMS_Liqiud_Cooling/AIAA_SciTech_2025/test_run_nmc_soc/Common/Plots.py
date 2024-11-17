'''
# Plots.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------   
from RCAIDE.Library.Plots                                           import *  

# ----------------------------------------------------------------------
#   Plot Results
# ----------------------------------------------------------------------  
def plot_results(results,save_figure_flag):  
    
     #plot_propulsor_throttles(results)
    
    #plot_flight_conditions(results) 
    
    #plot_aerodynamic_forces(results)

    #plot_aerodynamic_coefficients(results)  
    
    #plot_aircraft_velocities(results)
    
    #plot_bus_conditions(results)
    
    plot_battery_cell_conditions(results)
    
    plot_thermal_management_component(results)

    plot_battery_degradation(results)

    #plot_rotor_conditions(results) 

    #plot_electric_propulsor_efficiencies(results) 
    
    #plot_bus_C_rates(results)
                        
    return
