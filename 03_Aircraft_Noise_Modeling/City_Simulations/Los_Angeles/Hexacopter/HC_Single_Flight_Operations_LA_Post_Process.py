#----------------------------------------------------------------------
#   Imports
# --------------------------------------------------------------------- 
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Plots import *        
from RCAIDE import  save     
from RCAIDE import  load   

# python imports
import time as t
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
    routes_filepath = relative_path  +  '..' + separator + 'UAM_City_Routes.xlsx'
    topography_file = relative_path +  '..' + separator +  'Topography' + separator + 'LA_Metropolitan_Area.txt'
    flight_data     = pd.read_excel(routes_filepath,sheet_name=['Los_Angeles'])
    LA_flight_data  =  flight_data['Los_Angeles']
    route_count     = {}
    
    ti = t.time()
    
    operation_flight_times = np.array(['06:00:00',
                                       '07:00:00',
                                       '08:00:00',
                                       #'09:00:00',
                                       #'10:00:00',
                                       #'11:00:00', 
                                       #'12:00:00',
                                       #'13:00:00',
                                       #'14:00:00',
                                       #'15:00:00',
                                       #'16:00:00',
                                       #'17:00:00',
                                       #'18:00:00',
                                       #'19:00:00',
                                       #'20:00:00',
                                       #'21:00:00', 
                                       ])
    operation_period  = ['06:00:00','22:00:00']
         

    mic_x_res                       = 600
    mic_y_res                       = 1350
    aircraft_code                   = 'HC'
    city_code                       = 'LA' 
    cruise_altitude                 = 1000 * Units.feet
    noise_evaluation_pitch          = 300 * Units.feet
    number_of_microphone_in_stencil = 25
     

    filename_list_name =  aircraft_code + '_' + city_code +  '_Single_Flights_Raw'
    file_name_dict     =  load(filename_list_name + '.res' )
    processed_filename_list  = []
             
    for filename in file_name_dict.filename_list:  
        results = load(filename + '.res')
        
        origin_code = filename.split('_')[3]
        destination_code = filename.split('_')[4]
        
        if  origin_code in route_count:
            route_count[origin_code] += 1
        else:
            route_count[origin_code] = 1
            
        # Update flight departure time
        for i in range(len(operation_flight_times)):
            time = operation_flight_times[i]
            time = time.split(':')
            time[1] = int(time[1]) + 5 * (route_count[origin_code] - 1)
            if time[1] < 10:
                time[1] = '0' + str(time[1])
            else:
                time[1] = str(time[1])
            operation_flight_times[i] = time[0] +':' + time[1] + ':'+ time[2]
        
        # post process noise
        processed_noise_data =  post_process_noise_data(results, topography_file, operation_flight_times, operation_period ,number_of_microphone_in_stencil, noise_evaluation_pitch,  mic_x_res, mic_y_res)
          
        # save data in json file 
        processed_filename =  'Processed_'  +  aircraft_code + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + str(int(round(cruise_altitude/Units.feet,0)))+ 'ft'
        
        print("saving results")
        save(processed_noise_data, processed_filename + '.res')
        processed_filename_list.append(processed_filename)
        
    processed_filename_list_name =  aircraft_code + '_' + city_code +  '_Single_Flights_Processed'
    F =  Data(filename_list_name=processed_filename_list_name)
    save(F, processed_filename_list_name + '.res')
    tf = t.time() 
    print("Total time: "+str(tf-ti))
    return     

 
if __name__ == '__main__': 
    main()     