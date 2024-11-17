#----------------------------------------------------------------------
#   Imports
# --------------------------------------------------------------------- 
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Plots import *       
from RCAIDE import  load 
from RCAIDE import  save
import  pickle

# python imports 
import os 
import pickle 
import pandas as pd
import numpy as np
 

# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def main():           
           
    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator 
    routes_filepath = relative_path +  '..' + separator +  '..' + separator + 'UAM_City_Routes.xlsx' 
    flight_data     = pd.read_excel(routes_filepath,sheet_name=['Los_Angeles'])
    LA_flight_data  =  flight_data['Los_Angeles']
     
    
    aircraft_code             = 'HC'
    city_code                 = 'LA' 
    cruise_altitude           = 1500*Units.feet 
    mic_x_res                 = 1200
    mic_y_res                 = 1600     
    
    # create data structures to store noise 
    L_eq_T     = np.empt((0,mic_x_res,mic_y_res))      
    L_eq_24hr  = np.empt((0,mic_x_res,mic_y_res))      
    L_dn       = np.empt((0,mic_x_res,mic_y_res))      
    L_max      = np.empt((0,mic_x_res,mic_y_res))
    
    # create data structure to store energy comsumed at each airport 
    unique_airports =  0 # TO CORRECT 
    E_origin        = np.zeros(len(unique_airports))
    
    for i in  range(len(LA_flight_data)):
        # Extract Data
        origin_code       = LA_flight_data['Origin Code'][i]   
        destination_code  = LA_flight_data['Destination Code'][i]
         
        try:  
            # load data
            processed_filename  =  'Processed_'  +  aircraft_code + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + str(cruise_altitude)    # Aircraft_City_Frequency_Origin_Destination_Altitude
            PND =  load(processed_filename + '.res') # this is a json file 
            
            # concatenate  L_eq,L_eq_24hr, L_dn
            L_eq_T    = np.concatenate((L_eq_T   ,PND.L_eq), axis = 0)
            L_eq_24hr = np.concatenate((L_eq_24hr,PND.L_eq_24hr), axis = 0)
            L_dn      = np.concatenate((L_dn     ,PND.L_dn), axis = 0)
            L_max     = np.concatenate((L_max    ,PND.L_max ), axis = 0)
            
            # add the energy consumed by each origin airport
            loc           =  np.where( unique_airports,origin_code)
            E_origin[loc] += PND.Route_Energy_Consumed
        except:
            pass 
    
    
    # Add all routes
    Total_L_eq      =  0 #NEEDS TO BE UPDATES 
    Total_L_eq_24hr =  0 #NEEDS TO BE UPDATES 
    Total_L_dn      =  0 #NEEDS TO BE UPDATES 
    Total_L_max     =  0 #NEEDS TO BE UPDATES 
    
    # create data structure
    cumulative_PND = Data(
        Total_L_eq = Total_L_eq,
        Total_L_eq_24hr = Total_L_eq_24hr, 
        Total_L_dn = Total_L_dn, 
        Total_L_max = Total_L_max,
        Total_Energy = E_origin, 
        Airports      = unique_airports, 
    )
    
    # save data 
    cumulative_processed_filename =  'Cumulative_'  +  aircraft_code + '_' + city_code + '_' + str(cruise_altitude)    # Aircraft_City_Frequency_Origin_Destination_Altitude
    save(cumulative_PND, cumulative_processed_filename + '.res')    
    return     
 
if __name__ == '__main__': 
    main()     