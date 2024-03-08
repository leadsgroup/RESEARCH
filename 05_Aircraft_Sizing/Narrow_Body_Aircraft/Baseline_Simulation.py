 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from RCAIDE.Core import Units     
from RCAIDE.Visualization         import *     
from RCAIDE.Methods.Noise.Metrics import *  

# python imports 
import numpy as np   
import matplotlib.pyplot as plt   

import Vehicles
import Analyses
import Missions 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(): 
    # vehicle data
    vehicle  = Vehicles.vehicle_setup() 
    
    # Set up vehicle configs
    configs  = Vehicles.configs_setup(vehicle)

    # create analyses
    analyses = Analyses.setup(configs)
 
    # create mission instances (for multiple types of missions)
    missions = Missions.setup(analyses) 
     
    # mission analysis 
    baseline_results = missions.base_mission.evaluate()  
    TOL_results      = missions.takeoff_and_landing_noise_mission.evaluate()  
    
    # plot the results
    plot_baseline_mission(baseline_results) 
    plot_noise_mission(TOL_results) 
        
    return 

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------
 
def plot_baseline_mission(results):
    """This function plots the results of the mission analysis and saves those results to 
    png files."""

    # Plot Flight Conditions 
    plot_flight_conditions(results)
    
    # Plot Aerodynamic Forces 
    plot_aerodynamic_forces(results)
    
    # Plot Aerodynamic Coefficients 
    plot_aerodynamic_coefficients(results)
    
    # Plot Static Stability Coefficients 
    plot_stability_coefficients(results)    
    
    # Drag Components
    plot_drag_components(results)
    
    # Plot Altitude, sfc, vehicle weight 
    plot_altitude_sfc_weight(results)
    
    # Plot Velocities 
    plot_aircraft_velocities(results)  
    return

def plot_noise_mission(results):
    """This function plots the results of the mission analysis and saves those results to 
    png files.""" 
    flight_times = np.array(['06:00:00','07:00:00','08:00:00','09:00:00','10:00:00','11:00:00','12:00:00','13:00:00','14:00:00','15:00:00'])  
      
    noise_data      = post_process_noise_data(results)   
    noise_data      = DNL_noise_metric(noise_data, flight_times,time_period = 24*Units.hours)
    noise_data      = Equivalent_noise_metric(noise_data, flight_times,time_period = 15*Units.hours)
    noise_data      = SENEL_noise_metric(noise_data, flight_times,time_period = 24*Units.hours)  
    
    # Maximum Sound Pressure Level   
    plot_2D_noise_contour(noise_data,
                        noise_level      = np.max(noise_data.SPL_dBA,axis=0), 
                        min_noise_level  = 35,  
                        max_noise_level  = 90, 
                        noise_scale_label= 'SPL [dBA]',
                        save_filename    = "SPL_max_Noise_2D_Contour",
                        show_elevation   = True,
                        use_lat_long_coordinates= False,)   
                        

    # Day Night Average Noise Level 
    plot_2D_noise_contour(noise_data,
                        noise_level      = noise_data.DNL,
                        min_noise_level  = 35,  
                        max_noise_level  = 90, 
                        noise_scale_label= 'DNL',
                        save_filename    = "DNL_Noise_2D_Contour")  
        
    return
 
if __name__ == '__main__': 
    main()
    plt.show()