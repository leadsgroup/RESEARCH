# TR_Study.py
# Created: Mar 2023 

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
from MARC.Visualization.Geometry.Three_Dimensional.plot_3d_vehicle          import plot_3d_vehicle 
from MARC.Methods.Noise.Fidelity_Zero.Noise_Tools.generate_microphone_points import generate_terrain_elevated_microphone_points
from MARC.Visualization.Geometry                                            import *

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
    resize_aircraft            = False
    control_points             = 30              # number of control points per segment  

    # noise analysis parameters 
    run_noise_model            = True   
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
        meta_data = Data(
            ground_microphone_x_resolution          = 51  ,
            ground_microphone_y_resolution          = 31  ,       
            ground_microphone_min_x                 = 0 ,             
            ground_microphone_max_x                 = 1000  ,
            ground_microphone_min_y                 = -500,       
            ground_microphone_max_y                 = 500,
            true_course_angle                       = 0.0, 
            ground_microphone_x_stencil             = 3,
            ground_microphone_y_stencil             = 15)

        configs_analyses  = Analyses.analyses_setup(configs,run_noise_model,use_topology_flag,meta_data)

        # -------------------------------------------------------------------------------------------    
        # SET UP MISSION PROFILE  
        # -------------------------------------------------------------------------------------------    
        base_mission      = Missions.flyover_at_200ft_mission_setup(configs_analyses,vehicle,simulated_days,flights_per_day,control_points,recharge_battery,meta_data)
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
        filename          = 'TR_Flyover_200ft'
        save_results(noise_results,filename)   

    else:
        filename          = 'TR_Flyover_200ft'
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
     
