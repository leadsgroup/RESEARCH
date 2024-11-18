#----------------------------------------------------------------------
#   Imports
# --------------------------------------------------------------------- 
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Plots import *        
from RCAIDE import  save
import  pickle

# python imports 
import os 
import pickle
import sys 
import pandas as pd
import numpy as  np 

local_path_1 =  os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0] 

sys.path.append( os.path.join(local_path_1, 'Post_Processing_Functions'))
from Aircraft_Noise_Emissions   import post_process_noise_data 

# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def main():           
           
    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator 
    routes_filepath = relative_path +  '..' + separator +  '..' + separator + 'UAM_City_Routes.xlsx'
    topography_file = relative_path +  '..' + separator +  'Topography' + separator + 'LA_Metropolitan_Area.txt'
    flight_data     = pd.read_excel(routes_filepath,sheet_name=['Los_Angeles'])
    LA_flight_data  =  flight_data['Los_Angeles']
    
    
    operation_flight_times = np.array(['06:00:00','06:15:00','06:30:00','06:45:00',
                                       '07:00:00','07:15:00','07:30:00','07:45:00',
                                       '08:00:00','08:15:00','08:30:00','08:45:00',
                                       '09:00:00','09:15:00','09:30:00','09:45:00',
                                       '10:00:00','10:15:00','10:30:00','10:45:00',
                                       '11:00:00','11:15:00','11:30:00','11:45:00', 
                                       '12:00:00','12:15:00','12:30:00','12:45:00',
                                       '13:00:00','13:15:00','13:30:00','13:45:00',
                                       '14:00:00','14:15:00','14:30:00','14:45:00',
                                       '15:00:00','15:15:00','15:30:00','15:45:00',
                                       '16:00:00','16:15:00','16:30:00','16:45:00',
                                       '17:00:00','17:15:00','17:30:00','17:45:00',
                                       '18:00:00','18:15:00','18:30:00','18:45:00',
                                       '19:00:00','19:15:00','19:30:00','19:45:00',
                                       '20:00:00','20:15:00','20:30:00','20:45:00',
                                       '21:00:00', ])
    operation_period  = ['06:00:00','22:00:00']
         

    mic_x_res                 = 1200
    mic_y_res                 = 1600 
    noise_timesteps           = 225   
    aircraft_code             = 'HC'
    city_code                 = 'LA' 
    cruise_altitude           = 1500*Units.feet
    
    for i in  range(len(LA_flight_data)):
        # Extract Data
        origin_code       = LA_flight_data['Origin Code'][i]   
        destination_code  = LA_flight_data['Destination Code'][i]
         
        try: 
            filename =  aircraft_code +'_mission' + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + str(cruise_altitude)   
            results = load_results(filename)
           
            # post process noise
            processed_noise_data =  post_process_noise_data(results, topography_file, operation_flight_times, operation_period , noise_timesteps,  mic_x_res, mic_y_res)
              
            # save data in json file 
            processed_filename =  'Processed_'  +  aircraft_code + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + str(cruise_altitude)
            
            # get total energy comsumed for each route and append it onto noise 
            Total_Energy = 0
            for network in results.segments[0].analyses.energy.vehicle.networks: 
                busses  = network.busses 
                for bus in busses:
                    for battery_module in enumerate(bus.battery_modules): 
                        E_start_module  = results.segments[0].conditions.energy[bus.tag].battery_modules[battery_module.tag]  
                        E_end_module    = results.segments[-1].conditions.energy[bus.tag].battery_modules[battery_module.tag]    
                        Total_Energy    += E_start_module - E_end_module  
            processed_filename.Route_Energy_Consumed = Total_Energy
            
            save(processed_noise_data, processed_filename + '.res')
        except:
            pass 
             
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