'''
# Simulation_Repeated_Flight_Operations.py
#
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units, Data   
import pickle
from RCAIDE.Library.Plots                                           import *  

import time  
import numpy as np
import pylab as plt 
import sys 
sys.path.append('Common')  

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
    recharge_battery           = True            # flag to simulate battery recharge  
    plot_mission               = True            # plot mission flag  
    resize_aircraft            = True
    aircraft_range             = 50 *Units.nmi   # total ground distance  
    

    if RUN_NEW_MODEL_FLAG:    
         
        # -------------------------------------------------------------------------------------------    
        # SET UP VEHICLE
        # -------------------------------------------------------------------------------------------  
        vehicle = Vehicle.vehicle_setup(resize_aircraft, 'e-Twin-Otter')  
    
  
        # -------------------------------------------------------------------------------------------    
        # SET UP MISSION PROFILE  
        # -------------------------------------------------------------------------------------------    
        base_mission      = Missions.repeated_flight_operation_setup(vehicle,simulated_days,flights_per_day,recharge_battery)
        missions_analyses = Missions.missions_setup(base_mission)
       
    
        # -------------------------------------------------------------------------------------------    
        # RUN SIMULATION !!
        # -------------------------------------------------------------------------------------------
        results     = missions_analyses.base.evaluate() 
        
        # -------------------------------------------------------------------------------------------    
        # SAVE RESULTS
        # -------------------------------------------------------------------------------------------
        filename          = 'CTOL_Baseline'
        save_results(results,filename)   
    
    else:
        filename          = 'CTOL_Baseline'
        results = load_results(filename) 
        
    if plot_mission: 
        Plots.plot_results(results,save_figure_flag = True)       
    
    
    tf = time.time() 
    print ('time taken: '+ str(round(((tf-ti)/60),3)) + ' mins')   
    
    
    elapsed_range = results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]
    print('True Range     : ' + str(round(meta_data.flight_range/Units.nmi,2))  + ' nmi')   
    print('Computed Range : ' + str(round(elapsed_range/Units.nmi,2)) + ' nmi')   
        
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
     

