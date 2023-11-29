'''
# Simulation_Approach_and_Departure.py
#
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
 
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
    plot_mission               = True            # plot mission flag  
    resize_aircraft            = True
    control_points             = 10              # number of control points per segment  
    
    # noise analysis parameters 
    run_noise_model            = True 
    use_topology_flag          = False 
    
    flight_no                  = 0  # CHANGE ON SERVER 
        
    departure_angles           = [0,45,90,135,180]   
    microphone_terrain_data    = Data(ground_microphone_x_resolution          = 41  ,
                                      ground_microphone_y_resolution          = 41  ,   
                                      ground_microphone_min_x                 = 0 ,                
                                      ground_microphone_max_x                 = 2 * Units.nmi ,
                                      ground_microphone_min_y                 = -1 * Units.nmi ,   
                                      ground_microphone_max_y                 = 1 * Units.nmi + 1E-2, 
                                      ground_microphone_x_stencil             = 15,
                                      ground_microphone_y_stencil             = 15)  
    
    

    for flight_no in range(len(departure_angles)):    
        airport_geospacial_data    =  Data(true_course_angle                       = departure_angles[flight_no] * Units.degrees)
        
        
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
            base_mission      = Missions.approach_departure_mission_setup(configs_analyses,vehicle,simulated_days,flights_per_day,control_points,recharge_battery,microphone_terrain_data,airport_geospacial_data)
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
            filename          = 'SR_TG_Noise_Angle_' + str(departure_angles[flight_no]) 
            save_results(noise_results,filename)    
    
            elapsed_range = noise_results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]    
            print('Computed Range : ' + str(elapsed_range/Units.nmi) + ' nmi')    
    
        else: 
            filename          = 'SR_TG_Noise_Angle_' + str(departure_angles[flight_no]) 
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
     

