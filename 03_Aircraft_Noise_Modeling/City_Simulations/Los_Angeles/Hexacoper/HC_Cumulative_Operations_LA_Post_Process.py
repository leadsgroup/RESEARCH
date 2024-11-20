#----------------------------------------------------------------------
#   Imports
# --------------------------------------------------------------------- 
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Methods.Noise.Common.decibel_arithmetic   import SPL_arithmetic  
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
    Total_L_eq       = np.empt((mic_x_res,mic_y_res))      
    Total_L_eq_24hr  = np.empt((mic_x_res,mic_y_res))      
    Total_L_dn       = np.empt((mic_x_res,mic_y_res))      
    Total_L_max      = np.empt((mic_x_res,mic_y_res))
    
    # create data structure to store energy comsumed at each airport 
    unique_airports = LA_flight_data['Origin Code'].unique()
    E_origin        = np.zeros(len(unique_airports))
    
    processed_filename_list_name =  aircraft_code + '_' + city_code +  '_Single_Flights_Processed' 
    filename_list                = load(processed_filename_list_name)    
    for i in  range(len(LA_flight_data)):
        # Extract Data
        origin_code       = LA_flight_data['Origin Code'][i]    
         
        for filename in filename_list: 
            # load data 
            PND =  load(filename + '.res') # this is a json file 
            
            # concatenate  L_eq,L_eq_24hr,L_dn
            Total_L_eq      = SPL_arithmetic( np.concatenate((Total_L_eq[:,:,None]     ,PND.L_eq[:,:,None]), axis = 2), sum_axis=2)
            Total_L_eq_24hr = SPL_arithmetic(np.concatenate((Total_L_eq_24hr[:,:,None],PND.L_eq_24hr[:,:,None]), axis = 2), sum_axis=2)
            Total_L_dn      = SPL_arithmetic(np.concatenate((Total_L_dn[:,:,None]     ,PND.L_dn[:,:,None]), axis = 2), sum_axis=2)
            Total_L_max     = np.max(np.concatenate((Total_L_max[:,:,None]    ,PND.L_max[:,:,None] ), axis = 2), axis=2)
            
            # add the energy consumed by each origin airport
            loc           =  np.where( unique_airports,origin_code)
            E_origin[loc] += PND.Route_Energy_Consumed   
    
    # create data structure
    cumulative_PND = Data(
        Total_L_eq      = Total_L_eq,
        Total_L_eq_24hr = Total_L_eq_24hr, 
        Total_L_dn      = Total_L_dn, 
        Total_L_max     = Total_L_max,
        Total_Energy    = E_origin, 
        Airports        = unique_airports, 
    )
    
    # save data 
    cumulative_processed_filename =  'Cumulative_'  +  aircraft_code + '_' + city_code + '_' + str(cruise_altitude)    # Aircraft_City_Frequency_Origin_Destination_Altitude
    save(cumulative_PND, cumulative_processed_filename + '.res')    
    return     
 
if __name__ == '__main__': 
    main()     