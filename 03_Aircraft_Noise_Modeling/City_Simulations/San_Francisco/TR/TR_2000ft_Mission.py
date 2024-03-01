 
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Core import Units, Data   
import pickle
from RCAIDE.Visualization.Performance.Aerodynamics.Vehicle                    import *  
from RCAIDE.Visualization.Performance.Mission                                 import *  
from RCAIDE.Visualization.Performance.Aerodynamics.Rotor                      import *      
from RCAIDE.Visualization.Performance.Energy.Battery                          import *   
from RCAIDE.Visualization.Performance.Noise                                   import *  
from RCAIDE.Visualization.Topography                                          import * 
from RCAIDE.Visualization.Geometry.Three_Dimensional.plot_3d_vehicle          import plot_3d_vehicle 
from RCAIDE.Methods.Noise.Fidelity_Zero.Noise_Tools.generate_microphone_points import generate_terrain_elevated_microphone_points
from RCAIDE.Methods.Missions.compute_point_to_point_geospacial_data           import compute_point_to_point_geospacial_data
from RCAIDE.Visualization.Geometry                                            import *

import time  
import numpy as np
import pylab as plt 
import sys 
sys.path.append('../../../../Aircraft_Models/Tiltrotor/Common')  

import Vehicle
import Analyses 
import Missions
import Plots  
try:
    #import vsp 
    from RCAIDE.Input_Output.OpenVSP.vsp_write import write 
except ImportError:
    # This allows RCAIDE to build without OpenVSP
    pass   

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():  
    # start simulation clock
    ti                         = time.time()
    RUN_NEW_MODEL_FLAG         = True
    
    # -------------------------------------------------------------------------------------------    
    # SET UP SIMULATION PARAMETERS   
    # -------------------------------------------------------------------------------------------  
    simulated_days             = 1               # number of days simulated 
    flights_per_day            = 1               # number of flights per day   
    recharge_battery           = False           # flag to simulate battery recharge  
    plot_mission               = True            # plot mission flag  
    resize_aircraft            = False 
    control_points             = 6               # number of control points per segment  
    
    # noise analysis parameters 
    run_noise_model            = False
    use_topology_flag          = True  

    # tag for city 
    city                       =  
    topography_file            =  
    
    # departure airport location  
    departure_location         =  
    
    # coordinates for departure airport 
    departure_coord            = 
    
    # destination airport location
    destination_location       =  
    
    # coordinates for destination airport 
    destination_coord          =  

    if RUN_NEW_MODEL_FLAG:     
        for leg in range(len(departure_location)):   
            # -------------------------------------------------------------------------------------------    
            # SET UP VEHICLE
            # -------------------------------------------------------------------------------------------  
            vehicle = Vehicle.vehicle_setup(resize_aircraft)  
        
            # -------------------------------------------------------------------------------------------    
            # SET UP CONFIGURATIONS 
            # -------------------------------------------------------------------------------------------
            configs           = Vehicle.configs_setup(vehicle)  
    
            microphone_terrain_data =  generate_terrain_elevated_microphone_points(topography_file   = topography_file,
                                                                   ground_microphone_x_resolution    = 200,  
                                                                   ground_microphone_y_resolution    = 500, 
                                                                   ground_microphone_x_stencil       = 20,   
                                                                   ground_microphone_y_stencil       = 20)    
            
           
            airport_geospacial_data            =  compute_point_to_point_geospacial_data(topography_file  = topography_file,
                                                        departure_tag                         = departure_location[leg],
                                                        destination_tag                       = destination_location[leg],
                                                        departure_coordinates                 = [departure_coord[leg,0], departure_coord[leg,1]],
                                                        destination_coordinates               = [destination_coord[leg,0], destination_coord[leg,1]])
            
            configs_analyses  = Analyses.analyses_setup(configs,run_noise_model,use_topology_flag,microphone_terrain_data,airport_geospacial_data)
        
            # -------------------------------------------------------------------------------------------    
            # SET UP MISSION PROFILE  
            # -------------------------------------------------------------------------------------------    
            base_mission      = Missions.direct_mission_setup_at_2000ft(configs_analyses,vehicle,simulated_days,flights_per_day,control_points,recharge_battery,microphone_terrain_data,airport_geospacial_data)
            missions_analyses = Missions.missions_setup(base_mission) 
        
            # -------------------------------------------------------------------------------------------    
            # DEFINE ANALYSES 
            # -------------------------------------------------------------------------------------------
            analyses          = RCAIDE.Analyses.Analysis.Container()
            analyses.configs  = configs_analyses
            analyses.missions = missions_analyses 
            
        
            # -------------------------------------------------------------------------------------------    
            # FINALIZE SIMULATION 
            # -------------------------------------------------------------------------------------------    
            configs.finalize()
            analyses.finalize()   
        
            # -------------------------------------------------------------------------------------------    
            # APPEND MISSION TO SIMULATION 
            # -------------------------------------------------------------------------------------------    
            mission           = analyses.missions.base
            
        
            # -------------------------------------------------------------------------------------------    
            # RUN SIMULATION !!
            # -------------------------------------------------------------------------------------------
            noise_results     = mission.evaluate() 
        
            # -------------------------------------------------------------------------------------------    
            # SAVE RESULTS
            # -------------------------------------------------------------------------------------------
            filename          = 'TR_2000ft_mission_' + city + '_' +  departure_location[leg] + '_to_' + destination_location[leg]
            
            if plot_mission: 
                plot_elevation_contours(topography_file, use_lat_long_coordinates = False,save_filename = filename) 
                        
            save_results(noise_results,filename)    
        
            elapsed_range = noise_results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]
            print('True Range     : ' + str(round(airport_geospacial_data.flight_range/Units.nmi,2))  + ' nmi')   
            print('Computed Range : ' + str(round(elapsed_range/Units.nmi,2)) + ' nmi')   
            print('Range Error    : ' + str(round(elapsed_range-airport_geospacial_data.flight_range,2)) + ' m or ' +  str(round((elapsed_range-airport_geospacial_data.flight_range)/Units.nmi,2))  + ' nmi')     
            
        else:
            filename      = 'TR_2000ft_mission_' + city + '_' +  departure_location[leg] + '_to_' + destination_location[leg]
            noise_results = load_results(filename) 
    
        if plot_mission: 
            Plots.plot_results(noise_results,run_noise_model,save_figure_flag = True)    

    tf = time.time() 
    print ('time taken: '+ str(round(((tf-ti)/60),3)) + ' mins')   
    
        
    return 


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
     
