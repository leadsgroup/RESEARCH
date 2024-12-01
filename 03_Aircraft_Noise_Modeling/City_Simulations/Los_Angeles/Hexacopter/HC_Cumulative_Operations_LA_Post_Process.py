#----------------------------------------------------------------------
#   Imports
# --------------------------------------------------------------------- 
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Methods.Noise.Common.decibel_arithmetic   import SPL_arithmetic  
from RCAIDE.Library.Plots import *
from RCAIDE.Library.Methods.Noise.Common.background_noise     import background_noise

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
    storage_dir = "/home/aidanrm2/storage/noise/hexacopter_single_mission_results/"   
    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator 
    routes_filepath = relative_path  +  '..' + separator + 'UAM_City_Routes.xlsx' 
    flight_data     = pd.read_excel(routes_filepath,sheet_name=['Los_Angeles'])
    LA_flight_data  =  flight_data['Los_Angeles']
     
    
    aircraft_code             = 'HC'
    city_code                 = 'LA' 
    cruise_altitude           = 1000*Units.feet 
    mic_x_res                 = 1200
    mic_y_res                 = 2700
    
    # create data structures to store noise 
    Total_L_eq       = np.empty((mic_x_res,mic_y_res))      
    Total_L_eq_24hr  = np.empty((mic_x_res,mic_y_res))      
    Total_L_dn       = np.empty((mic_x_res,mic_y_res))      
    Total_L_max      = np.empty((mic_x_res,mic_y_res))
    
    # create data structure to store energy comsumed at each airport 
    unique_airports = LA_flight_data['Origin Code'].unique()
    E_origin        = np.zeros(len(unique_airports))
    
    processed_filename_list_name =  aircraft_code + '_' + city_code +  '_Single_Flights_Processed' 
    file_name_dict               =  load(processed_filename_list_name + '.res')
    
    for filename in file_name_dict.filename_list:     
        # load data

        origin_code = filename.split('_')[3]
        
        PND =  load_results(filename,storage_dir)  
        
        L_eq = PND.L_eq[:,:,None]
        L_max = PND.L_max[:,:,None]
        L_dn = PND.L_dn[:,:,None]
        L_eq_24hr = PND.L_eq_24hr[:,:,None]
        mask = L_dn < (background_noise() + 0.1)
        L_dn[mask] = 0
        L_eq[mask] = 0
        L_max[mask] = 0
        L_eq_24hr[mask] = 0
        
        # concatenate  L_eq,L_eq_24hr,L_dn
        Total_L_eq      = SPL_arithmetic( np.concatenate((Total_L_eq[:,:,None]     ,L_eq), axis = 2), sum_axis=2)
        Total_L_eq_24hr = SPL_arithmetic(np.concatenate((Total_L_eq_24hr[:,:,None],L_eq_24hr), axis = 2), sum_axis=2)
        Total_L_dn      = SPL_arithmetic(np.concatenate((Total_L_dn[:,:,None]     ,L_dn), axis = 2), sum_axis=2)
        Total_L_max     = np.max(np.concatenate((Total_L_max[:,:,None]    ,L_max ), axis = 2), axis=2)
        
        # add the energy consumed by each origin airport
        loc           =  list(unique_airports).index(origin_code)
        E_origin[loc] += PND.energy_consumed   
    
    # create data structure
    cumulative_PND = Data(
        Total_L_eq      = Total_L_eq,
        Total_L_eq_24hr = Total_L_eq_24hr, 
        Total_L_dn      = Total_L_dn, 
        Total_L_max     = Total_L_max,
        Total_Energy    = E_origin, 
        Airports        = unique_airports, 
    )
    
    mask = Total_L_dn < (background_noise() + 0.1)
    Total_L_dn[mask] = background_noise()
    Total_L_eq[mask] = background_noise()
    Total_L_max[mask] = background_noise()
    Total_L_eq_24hr[mask] = background_noise()    
    
    # save data 
    cumulative_processed_filename =  'Cumulative_'  +  aircraft_code + '_' + city_code + '_' + str(int(np.ceil(round(cruise_altitude/Units.feet)))) +'ft'    # Aircraft_City_Frequency_Origin_Destination_Altitude
    save_results(cumulative_PND, cumulative_processed_filename, storage_dir)    
    return     
 
def save_results(results, filename, storage_dir):
    save_dir = storage_dir
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  # Create the directory if it doesn't exist
    res_file = os.path.join(save_dir, f"{filename}.res")
    save(results, res_file)
    return

def load_results(filename, storage_dir=os.getcwd()):
    load_dir = storage_dir
    load_file = os.path.join(load_dir, f"{filename}.res")
    if os.path.exists(load_file):
        results = load(load_file)
        return results
    else:
        raise FileNotFoundError(f"File {load_file} not found.")
 
if __name__ == '__main__': 
    main()
    
