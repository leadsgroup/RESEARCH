#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE 
from RCAIDE.Framework.Core import Units, Data   
from  RCAIDE.Framework.Analyses.Geodesics.Geodesics import Calculate_Distance  
from RCAIDE.Library.Plots import *       
from RCAIDE import  load 
from RCAIDE import  save
import  pickle

# python imports 
import os 
import pickle
import sys 
import pandas as pd
import numpy as  np 

local_path_1 =  os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0]
local_path_2 =  os.path.split(os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0])[0]

sys.path.append( os.path.join(local_path_1, 'Flight_Path'))
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Hexacopter'))
from Hexacopter              import vehicle_setup, configs_setup
from compute_route_distances import compute_route_distances
from compute_terrain_points  import compute_terrain_points
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
    LA_flight_data         =  flight_data['Los_Angeles']
    
    
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
    mic_stencil               = 100
    aircraft_code             = 'HC'
    city_code                 = 'LA' 
    cruise_altitude           = 1500*Units.feet
    high_speed_climb_distance = 1000 # NEEDS BE UPDATED BASED ON AIRCRAFT 
    high_speed_descent_distance  = 1000 # NEEDS BE UPDATED BASED ON AIRCRAFT  
    radius_Vert1              = 3600*Units.feet # circular pattern radius around vertiport 1
    radius_Vert2              = 3600*Units.feet # circular pattern radius around vertiport 2
    dep_heading               = 200 *Units.degree # Heading [degrees] of the departure from vertiport 1
    app_heading               = 90  *Units.degree# Heading [degrees] of the approach to vertiport 2
    
    max_cruise_distance =  30 * Units.mile
     
    #with  open('results_data.pkl', 'rb') as  file:
        #results = pickle.load(file)            
    results = missions.base_mission.evaluate()
    
    filename =  aircraft_code +'_mission' + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + str(cruise_altitude)    # Aircraft_City_Frequency_Origin_Destination_Altitude            
    with  open((filename + '.pkl'), 'rb') as  file:
        pickle.load(filename)
   
    # post process noise 
    noise_data   = post_process_noise_data(results,
                                           flight_times = operation_flight_times,  
                                           time_period  = operation_period,
                                           evalaute_noise_metrics = False)  
  
    # save data
    filename_res =  aircraft_code + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + str(cruise_altitude)    # Aircraft_City_Frequency_Origin_Destination_Altitude
    save(noise_data, filename_res + '.res')
            
        
    return   
def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions 
 
if __name__ == '__main__': 
    main()     