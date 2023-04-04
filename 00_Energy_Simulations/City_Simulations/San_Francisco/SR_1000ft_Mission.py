 
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import MARC
from MARC.Core import Units, Data   
import pickle
from MARC.Visualization.Performance.Aerodynamics.Vehicle                    import *  
from MARC.Visualization.Performance.Mission                                 import *  
from MARC.Visualization.Performance.Aerodynamics.Rotor                      import *      
from MARC.Visualization.Performance.Energy.Battery                          import *   
from MARC.Visualization.Performance.Noise                                   import *  
from MARC.Visualization.Topography                                          import * 
from MARC.Visualization.Geometry.Three_Dimensional.plot_3d_vehicle          import plot_3d_vehicle 
from MARC.Methods.Noise.Fidelity_One.Noise_Tools.generate_microphone_points import generate_terrain_elevated_microphone_points
from MARC.Methods.Missions.compute_point_to_point_geospacial_data           import compute_point_to_point_geospacial_data
from MARC.Visualization.Geometry                                            import *

import time  
import numpy as np
import pylab as plt 
import sys 
sys.path.append('../../Aircraft_Models/Stopped_Rotor/Common')  

import Vehicle
import Analyses 
import Missions
import Plots  
try:
    #import vsp 
    from MARC.Input_Output.OpenVSP.vsp_write import write 
except ImportError:
    # This allows MARC to build without OpenVSP
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
    plot_mission               = False           # plot mission flag  
    resize_aircraft            = True 
    control_points             = 6               # number of control points per segment  
    
    # noise analysis parameters 
    run_noise_model            = True
    use_topology_flag          = True  

    # tag for city  
    city                       = # CHANGE TO ACRONYM OF CITY e.g. 'LA' for Los Angeles 
    topography_file            = # CHANGE TO FILE PATH OF MAP e.g '../../Maps_and_Scales/Los_Angeles/LA_Metropolitan_Area_5.txt'
    
    # departure airport location  
    departure_location         = # CREATE LIST OF CODES FOR DEPARTURE AIRPORTS. e.g. for Los Angeles  ['LAX','SNA','BUR','LAX','BUR','LAX']
    
    # coordinates for departure airport 
    departure_coord            = # CREATE LIST OF COORDIATES FOR DEPARTURE AIRPORTS. e.g. for Los Angeles  np.array([[ 33.94506045,-118.4106796 ],..., [34.19318238,-118.3541369],[ 33.94506045,-118.4106796 ],])
    
    # destination airport location
    destination_location       = # CREATE LIST OF CODES FOR DESTINATION AIRPORTS. e.g. for Los Angeles  ['SBD','ONT','ONT','BUR','DIS','DIS']
    
    # coordinates for destination airport 
    destination_coord          =  # CREATE LIST OF COORDIATES FOR DESTINATION AIRPORTS. e.g. for Los Angeles  np.array([[ 33.94506045,-118.4106796 ],..., [34.19318238,-118.3541369],[ 33.94506045,-118.4106796 ],])
    
    # To account climb and descent, we have to adjust the range of the cruise portion of flight 
    # Below are the adjusted distances (in nautical miles)to be subtracted from absolute distance 
    # between airport coordinates   
    adjusted_cruise_distance     = # RUN SCRIPT WITH ALL ZEROS FIRST TO GET CORRECT DISTANCES, THEN PLACE OUTPUT IN ARRAY np.array([0,0,...0])*Units.nmi

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
                                                                   ground_microphone_x_resolution    = 100,  
                                                                   ground_microphone_y_resolution    = 200, 
                                                                   ground_microphone_x_stencil       = 5,   
                                                                   ground_microphone_y_stencil       = 5)    
            
           
            airport_geospacial_data            =  compute_point_to_point_geospacial_data(topography_file  = topography_file,
                                                        departure_tag                         = departure_location[leg],
                                                        destination_tag                       = destination_location[leg],
                                                        departure_coordinates                 = [departure_coord[leg,0], departure_coord[leg,1]],
                                                        destination_coordinates               = [destination_coord[leg,0], destination_coord[leg,1]],
                                                        adjusted_cruise_distance              = adjusted_cruise_distance[leg])
            
            configs_analyses  = Analyses.analyses_setup(configs,run_noise_model,use_topology_flag,microphone_terrain_data,airport_geospacial_data)
        
            # -------------------------------------------------------------------------------------------    
            # SET UP MISSION PROFILE  
            # -------------------------------------------------------------------------------------------    
            base_mission      = Missions.full_mission_setup_at_1000ft(configs_analyses,vehicle,simulated_days,flights_per_day,control_points,recharge_battery,microphone_terrain_data,airport_geospacial_data)
            missions_analyses = Missions.missions_setup(base_mission) 
        
            # -------------------------------------------------------------------------------------------    
            # DEFINE ANALYSES 
            # -------------------------------------------------------------------------------------------
            analyses          = MARC.Analyses.Analysis.Container()
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
            filename          = 'LA_SR_1000ft_mission_' + city + '_' +  departure_location[leg] + '_to_' + destination_location[leg]
            
            if plot_mission: 
                plot_elevation_contours(topography_file, use_lat_long_coordinates = False, airport_geospacial_data  = airport_geospacial_data,save_filename = filename) 
                        
            save_results(noise_results,filename)    
        
            elapsed_range = noise_results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]
            print('True Range     : ' + str(round(airport_geospacial_data.flight_range/Units.nmi,2))  + ' nmi')   
            print('Computed Range : ' + str(round(elapsed_range/Units.nmi,2)) + ' nmi')   
            print('Range Error    : ' + str(round(elapsed_range-airport_geospacial_data.flight_range,2)) + ' m or ' +  str(round((elapsed_range-airport_geospacial_data.flight_range)/Units.nmi,2))  + ' nmi')     
            
        else:
            filename      = 'LA_SR_1000ft_mission_' + city + '_' +  departure_location[leg] + '_to_' + destination_location[leg]
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
     
