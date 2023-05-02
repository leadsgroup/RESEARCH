 
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
sys.path.append('../../Aircraft_Models/Hexacopter/Common')  

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
    plot_mission               = True            # plot mission flag  
    resize_aircraft            = True 
    control_points             = 6               # number of control points per segment  
    
    # noise analysis parameters 
    run_noise_model            = True         
    use_topology_flag          = True  

    # tag for city 
    city                       = 'SF'
    topography_file            = '../../Maps_and_Scales/San_Francisco/SF_Metropolitan_Area_1.txt'
    

    ''' Provided below is the information for 9 direct routes between two airports'''
    # departure airport location  
    departure_location         = ['SQL']  
    
    # coordinates for departure airport 
    departure_coord            = np.array([[37.51406,-122.251304]])
    
    # destination airport location
    destination_location       = ['OAK']
    
    # coordinates for destination airport 
    destination_coord          = np.array([[37.7125,-122.2197]])
    
    
    ''' Change the flight number i.e. "flight_no" below to simulate one of the 9 routes '''
    flight_no = 0    
     
    microphone_terrain_data =  generate_terrain_elevated_microphone_points(topography_file   = topography_file,
                                                           ground_microphone_x_resolution    = 225,  
                                                           ground_microphone_y_resolution    = 390, 
                                                           ground_microphone_x_stencil       = 2,   
                                                           ground_microphone_y_stencil       = 2)    
    
   
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
        filename          = 'HC_1000ft_mission_' + city + '_' +  departure_location[flight_no] + '_to_' + destination_location[flight_no] 
        save_results(noise_results,filename)    
    
        elapsed_range = noise_results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]
        print('True Range     : ' + str(airport_geospacial_data.flight_range/Units.nmi) + ' nmi')   
        print('Computed Range : ' + str(elapsed_range/Units.nmi) + ' nmi')   
        print('Range Error    : ' + str(elapsed_range-airport_geospacial_data.flight_range) + ' m or ' +  str((elapsed_range-airport_geospacial_data.flight_range)/Units.nmi)  + ' nmi')     
        
    else: 
        filename      = 'HC_1000ft_mission_' + city + '_' +  departure_location[flight_no] + '_to_' + destination_location[flight_no]
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
     
