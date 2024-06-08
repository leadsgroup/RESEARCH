'''
# Simulation_Repeated_Flight_Operations.py
#
# Created: May 2019, M. Clarke
#          Sep 2020, M. Clarke
#          May 2024, S S. Shekar

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
import pandas as pd 
import sys 
sys.path.append('Common')  

import Vehicle
import Analyses 
import Missions
import Plots  
 

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
    recharge_battery           = True            # flag to simulate battery recharge  
    plot_mission               = True            # plot mission flag  
    resize_aircraft            = False
   
    airport = 'ORD'
    month = 'July'
    
    simulation_meta_data_path =  'Simulation_Plan.xlsx'
    flight_no, idle_time, taxi_time, cruise_distance, mean_temperature, tms_operation = get_flight_data(simulation_meta_data_path, airport, month)
  

    if RUN_NEW_MODEL_FLAG:    
         
        # -------------------------------------------------------------------------------------------    
        # SET UP VEHICLE
        # -------------------------------------------------------------------------------------------  
        vehicle = Vehicle.vehicle_setup(resize_aircraft,'e_twin_otter_vehicle')  
    
       
        # -------------------------------------------------------------------------------------------    
        # SET UP MISSION PROFILE  
        # -------------------------------------------------------------------------------------------    
        base_mission      = Missions.repeated_flight_operation_setup(vehicle,simulated_days,flight_no, idle_time, taxi_time, cruise_distance, mean_temperature,tms_operation, recharge_battery, airport, month)
        missions_analyses = Missions.missions_setup(base_mission)
       
    
        # -------------------------------------------------------------------------------------------    
        # RUN SIMULATION !!
        # -------------------------------------------------------------------------------------------
        results     = missions_analyses.base.evaluate() 
        
        # -------------------------------------------------------------------------------------------    
        # SAVE RESULTS
        # -------------------------------------------------------------------------------------------
        filename          = 'e_Twin_Otter_' + airport  + '_' + month
        save_results(results,filename)   
    
    else:
        filename          = 'e_Twin_Otter_' + airport + '_' + month
        results = load_results(filename) 
        
    if plot_mission: 
        Plots.plot_results(results,save_figure_flag = False)       
    
    
    tf = time.time() 
    print ('time taken: '+ str(round(((tf-ti)/60),3)) + ' mins')   
    
    
    #elapsed_range = results.segments[-1].conditions.frames.inertial.aircraft_range[-1,0]
    #print('True Range     : ' + str(round(meta_data.flight_range/Units.nmi,2))  + ' nmi')   
    #print('Computed Range : ' + str(round(elapsed_range/Units.nmi,2)) + ' nmi')   
        
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

# ------------------------------------------------------------------
#   Load Simulation Data
# ------------------------------------------------------------------   
def get_flight_data(file_path, airport, month):
    
    df = pd.read_excel(file_path, sheet_name=airport)
    
    # Fill forward missing values in 'Airport' and 'Month' columns
    df['Month'].fillna(method='ffill', inplace=True)
    
    # Filter the dataframe based on airport and month
    filtered_df = df[(df['Month'] == month)]
    
    # Store the data into variables
    flight_no = filtered_df['Flight_No'].tolist()
    idle_time = filtered_df['Idle_Time'].tolist()
    taxi_time = filtered_df['Taxi_Time'].tolist()
    cruise_distance = filtered_df['Cruise_Distance'].tolist()
    mean_temperature = filtered_df['Mean_Temperature'].tolist()
    
    tms_operation = {
           'no_hex_operation_has': filtered_df['no_hex_operation_has'].tolist(),
           'no_hex_operation_hex': filtered_df['no_hex_operation_hex'].tolist(),
           'max_hex_operation_has': filtered_df['max_hex_operation_has'].tolist(),
           'max_hex_operation_hex': filtered_df['max_hex_operation_hex'].tolist(),
           'hex_low_alt_climb_operation_has': filtered_df['hex_low_alt_climb_operation_has'].tolist(),
           'hex_low_alt_climb_operation_hex': filtered_df['hex_low_alt_climb_operation_hex'].tolist(),
           'hex_high_alt_climb_operation_has': filtered_df['hex_high_alt_climb_operation_has'].tolist(),
           'hex_high_alt_climb_operation_hex': filtered_df['hex_high_alt_climb_operation_hex'].tolist(),
           'hex_cruise_operation_has': filtered_df['hex_cruise_operation_has'].tolist(),
           'hex_cruise_operation_hex': filtered_df['hex_cruise_operation_hex'].tolist(),
           'hex_descent_operation_has': filtered_df['hex_descent_operation_has'].tolist(),
           'hex_descent_operation_hex': filtered_df['hex_descent_operation_hex'].tolist(),
           'recharge_has': filtered_df['recharge_has'].tolist(),
           'recharge_hex': filtered_df['recharge_hex'].tolist(),
       }    
    
    
    return flight_no, idle_time, taxi_time, cruise_distance, mean_temperature, tms_operation
if __name__ == '__main__': 
    main()    
    plt.show()
     

