 
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
from RCAIDE.Methods.Noise.Fidelity_One.Noise_Tools.generate_microphone_points import generate_terrain_elevated_microphone_points
from RCAIDE.Methods.Missions.compute_point_to_point_geospacial_data           import compute_point_to_point_geospacial_data
from RCAIDE.Visualization.Geometry                                            import *

import time  
import numpy as np
import pylab as plt 
import sys 
sys.path.append('../../Aircraft_Models/Tiltrotor/Common')  

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
    resize_aircraft            = True
    control_points             = 6               # number of control points per segment  
    
    # noise analysis parameters 
    run_noise_model            = True
    use_topology_flag          = True  

       # tag for city 
    city                       = 'SF'
    topography_file            = # CHANGE TO FILE PATH OF MAP e.g '../../Maps_and_Scales/Los_Angeles/LA_Metropolitan_Area_5.txt'
    

    ''' Provided below is the information for 9 direct routes between two airports'''
    # departure airport location  
    departure_location         = # CREATE LIST OF CODES FOR DEPARTURE AIRPORTS. e.g. for Los Angeles  # ['LAX','SNA','BUR','LAX']  
    
    # coordinates for departure airport 
    departure_coord            =  # CREATE LIST OF COORDIATES FOR DEPARTURE AIRPORTS. e.g. for Los Angeles   # np.array([[33.9348,-118.39678], [33.81347,-118.15542], [34.1843,-118.36587], [33.9348,-118.39678]])
    
    # destination airport location
    destination_location       = # CREATE LIST OF CODES FOR DESTINATION AIRPORTS. e.g. for Los Angeles   ['SBD','ONT','ONT','BUR']
    
    # coordinates for destination airport 
    destination_coord          =  # CREATE LIST OF COORDIATES FOR DESTINATION AIRPORTS. e.g. for Los Angeles   # np.array([[34.0758,-117.2512], [34.04920,-117.59908],  [34.04920,-117.59908], [34.1843,-118.36587]])
    
    
    
    ''' Change the flight number i.e. "flight_no" below to simulate one of the 9 routes '''
    flight_no = 0    

    microphone_terrain_data =  generate_terrain_elevated_microphone_points(topography_file   = topography_file,
                                                           ground_microphone_x_resolution    = 225,  
                                                           ground_microphone_y_resolution    = 390, 
                                                           ground_microphone_x_stencil       = 5,   
                                                           ground_microphone_y_stencil       = 5)    
    
   
    airport_geospacial_data            =  compute_point_to_point_geospacial_data(topography_file  = topography_file,
                                                departure_tag                         = departure_location[flight_no],
                                                destination_tag                       = destination_location[flight_no],
                                                departure_coordinates                 = [departure_coord[flight_no,0], departure_coord[flight_no,1]],
                                                destination_coordinates               = [destination_coord[flight_no,0], destination_coord[flight_no,1]])
    if RUN_NEW_MODEL_FLAG:   
        # -------------------------------------------------------------------------------------------    
        # SET UP VEHICLE
        # -------------------------------------------------------------------------------------------  
        vehicle = Vehicle.vehicle_setup(resize_aircraft)  
    
        # -------------------------------------------------------------------------------------------    
        # SET UP CONFIGURATIONS 
        # -------------------------------------------------------------------------------------------
        configs           = Vehicle.configs_setup(vehicle)  

        
        configs_analyses  = Analyses.analyses_setup(configs,run_noise_model,use_topology_flag,microphone_terrain_data,airport_geospacial_data)
    
        # -------------------------------------------------------------------------------------------    
        # SET UP MISSION PROFILE  
        # -------------------------------------------------------------------------------------------    
        base_mission      = Missions.direct_mission_setup_at_1000ft(configs_analyses,vehicle,simulated_days,flights_per_day,control_points,recharge_battery,microphone_terrain_data,airport_geospacial_data)
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
        filename          = 'TR_1000ft_mission_' + city + '_' +  departure_location[flight_no] + '_to_' + destination_location[flight_no] 
        save_results(noise_results,filename)    
    
        elapsed_range = noise_results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]
        print('True Range     : ' + str(airport_geospacial_data.flight_range/Units.nmi) + ' nmi')   
        print('Computed Range : ' + str(elapsed_range/Units.nmi) + ' nmi')   
        print('Range Error    : ' + str(elapsed_range-airport_geospacial_data.flight_range) + ' m or ' +  str((elapsed_range-airport_geospacial_data.flight_range)/Units.nmi)  + ' nmi')     
        
    else: 
        filename      = 'TR_1000ft_mission_' + city + '_' +  departure_location[flight_no] + '_to_' + destination_location[flight_no]
        noise_results = load_results(filename) 

    if plot_mission: 
            plot_elevation_contours(topography_file, use_lat_long_coordinates = False,save_filename = filename) 
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
     
