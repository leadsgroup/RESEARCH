# Script to convert cumulative .res files to .csv tables for post processing

import  pandas as  pd
import  numpy as  np
import os
from RCAIDE import  load 
from RCAIDE import  save
import  sys

local_path_1 =  os.path.split(sys.path[0])[0]

sys.path.append( os.path.join(local_path_1, 'Post_Processing_Functions'))
from Aircraft_Noise_Emissions   import generate_terrain_microphone_locations


def main():
    store_dir = '/Users/aidanmolloy/Documents/LEADS/RESEARCH/03_Aircraft_Noise_Modeling/City_Simulations/Los_Angeles/Tilt_Stopped_Rotor' # Change this as needed
    filename = 'Cumulative_TRS_LA_1000ft' # Change this as needed
    
    microphone_x_resolution                 = 1200 
    microphone_y_resolution                 = 2700
    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator 
    topography_file = relative_path +  '..' + separator +  'City_Simulations' + separator + 'Los_Angeles' + separator + 'Topography' + separator + 'LA_Metropolitan_Area.txt'    
    
    microphone_locations , microphone_coordinates= generate_terrain_microphone_locations(topography_file, microphone_x_resolution, microphone_y_resolution)
    mic_loc =    np.reshape(microphone_locations, (microphone_x_resolution,microphone_y_resolution,3))
    mic_coord =     np.reshape(microphone_coordinates,   (microphone_x_resolution,microphone_y_resolution,3))
    
    flight_data = load_results(filename, store_dir)
    
    Total_L_dn = flight_data['Total_L_dn'].flatten(order='C')
    Total_L_max = flight_data['Total_L_max'].flatten(order='C')
    Total_L_eq = flight_data['Total_L_eq'].flatten(order='C')
    Total_L_eq_24hr = flight_data['Total_L_eq_24hr'].flatten(order='C')
    elevation = mic_coord[:,:, 2].flatten(order='C')
    Latitude = mic_coord[:, :, 0].flatten(order='C')
    Longitude = mic_coord[:, :, 1].flatten(order='C')
    SEL = np.zeros(np.shape(Latitude))
    L_Aeq_jetliner = np.zeros(np.shape(Latitude))
    
    numpy_coord_num_data = np.array([elevation, Latitude, Longitude, Total_L_max, Total_L_eq, Total_L_eq_24hr, SEL, Total_L_dn, L_Aeq_jetliner])
    numpy_coord_num_data = numpy_coord_num_data.transpose()
    
    coord_num_data = pd.DataFrame(numpy_coord_num_data, columns=['elevation', 'Latitude', 'Longitude', 'L_Amax', 'L_AeqT', 'L_AeqT_24hr', 'SEL', 'L_dn', 'L_Aeq_jetliner'])
    coord_num_data.to_csv((filename+'.csv'))
    
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