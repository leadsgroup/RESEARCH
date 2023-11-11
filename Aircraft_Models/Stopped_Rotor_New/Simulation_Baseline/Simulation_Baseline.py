'''
# Simulation_Baseline.py
#
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Core import Units, Data   
import pickle
from RCAIDE.Visualization import * 

import time  
import numpy as np
import pylab as plt 
import sys 
sys.path.append('../Common')  

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
    RUN_NEW_MODEL_FLAG         = False 
    
    # -------------------------------------------------------------------------------------------    
    # SET UP SIMULATION PARAMETERS   
    # -------------------------------------------------------------------------------------------  
    simulated_days             = 1               # number of days simulated 
    flights_per_day            = 1               # number of flights per day   
    recharge_battery           = False           # flag to simulate battery recharge  
    plot_mission               = True            # plot mission flag  
    resize_aircraft            = False 
    control_points             = 16              # number of control points per segment 
    aircraft_range             = 30 *Units.nmi   # total ground distance  
    
    # noise analysis parameters 
    run_noise_model            = False    
    use_topology_flag          = False 

    if RUN_NEW_MODEL_FLAG:    
         
        # -------------------------------------------------------------------------------------------    
        # SET UP VEHICLE
        # -------------------------------------------------------------------------------------------  
        vehicle = Vehicle.vehicle_setup(resize_aircraft)  
    
        # -------------------------------------------------------------------------------------------    
        # SET UP CONFIGURATIONS 
        # -------------------------------------------------------------------------------------------
        configs           = Vehicle.configs_setup(vehicle)   
    
        microphone_terrain_data    = Data(ground_microphone_x_resolution          = 61  ,
                                          ground_microphone_y_resolution          = 41  ,   
                                          ground_microphone_min_x                 = 0 ,                 
                                          ground_microphone_max_x                 = 6.2 * Units.nmi ,
                                          ground_microphone_min_y                 = -2 * Units.nmi ,   
                                          ground_microphone_max_y                 = 2 * Units.nmi + 1E-2, 
                                          ground_microphone_x_stencil             = 15,
                                          ground_microphone_y_stencil             = 20)  
        
        
        airport_geospacial_data    = Data(true_course_angle                       = 0.0 * Units.degrees,
                                          flight_range                            = aircraft_range)
         
        configs_analyses  = Analyses.analyses_setup(configs,run_noise_model,use_topology_flag,microphone_terrain_data,airport_geospacial_data)
    
        # -------------------------------------------------------------------------------------------    
        # SET UP MISSION PROFILE  
        # -------------------------------------------------------------------------------------------    
        base_mission      = Missions.baseline_mission_setup(configs_analyses,vehicle,simulated_days,flights_per_day,control_points,recharge_battery,airport_geospacial_data = airport_geospacial_data)
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
        filename          = 'SR_Baseline'
        save_results(noise_results,filename)   
        
        tf = time.time() 
        print ('time taken: '+ str(round(((tf-ti)/60),3)) + ' mins')    
        
        elapsed_range = noise_results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]
        print('True Range     : ' + str(round(airport_geospacial_data.flight_range/Units.nmi,2))  + ' nmi')   
        print('Computed Range : ' + str(round(elapsed_range/Units.nmi,2)) + ' nmi')   
            
    else:
        filename          = 'SR_Baseline'
        noise_results = load_results(filename) 
        
    if plot_mission: 
        Plots.plot_results(noise_results,run_noise_model,save_figure_flag = True)       
    
    
        
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
     

