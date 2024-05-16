import numpy as np
from RCAIDE.Framework.Core import Units , Data 
import matplotlib.pyplot as plt
from matplotlib.image import imread
from scipy.interpolate import griddata  
import matplotlib.colors
import matplotlib.colors as colors  
from geopy.distance import geodesic as GD 
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main(): 
    # tag for city  
    city                        = 'LA' 
    aircraft_models             = ['HC','SR','TR']
    
    # altitude 
    altitudes                   = ['1000','2000']
    
    # departure airport location  
    dep_loc         = ['LAX','LGB','BUR','LAX','BUR','LAX','DIS','SNA','SNA','SBD','ONT','ONT','BUR','DIS','DIS','SBD','BUR','SBD'] 
    
    # destination airport location
    des_loc         = ['SBD','ONT','ONT','BUR','DIS','DIS','SBD','BUR','SBD','LAX','LGB','BUR','LAX','BUR','LAX','DIS','SNA','SNA'] 
     
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            altitude      = altitudes[alt_i]
            aircraft      = aircraft_models[m_i]   
            

            for flight_no in range(len(dep_loc)):   
                # create file name path
                flight_path         =  dep_loc[flight_no] + '_to_' + des_loc[flight_no] 
                filename            = 'Raw_Data_' + aircraft + '/' + aircraft + '_' + altitude + 'ft_mission_' + city + '_' +  flight_path
                results             = load_results(filename)         
                  
                energy_per_flight   = results.segments[-1].conditions.propulsion.battery.pack.energy[0,0] - results.segments[-1].conditions.propulsion.battery.pack.energy[-1,0]
                
                print(aircraft)
                print(flight_path)
                print('Flight Time: ',  str(round(results.segments[-1].conditions.frames.inertial.time[-1,0]/60,2)) + ' mins')
                print('Consumed Energy : ' + str(round(energy_per_flight*0.000277778,4)) + ' Wh \n')
    return  

 
# ------------------------------------------------------------------
# SUPPORTING FUNCTIONS 
# ------------------------------------------------------------------

def get_topography_data(topography_file,N_lat,N_long ):  
    
    data = np.loadtxt(topography_file)
    Long = data[:,0]
    Lat  = data[:,1]
    Elev = data[:,2]  
    
    x_min_coord = np.min(Lat)
    x_max_coord = np.max(Lat)
    y_min_coord = np.min(Long)
    y_max_coord = np.max(Long)
    if np.min(Long)>180: 
        y_min_coord = np.min(Long)-360
    if np.max(Long)>180:
        y_max_coord = np.max(Long)-360  
    
    top_left_map_coords      = np.array([x_max_coord,y_min_coord])
    bottom_left_map_coords   = np.array([x_min_coord,y_min_coord])  
    bottom_right_map_coords  = np.array([x_min_coord,y_max_coord]) 
    
    x_dist_max = GD(top_left_map_coords,bottom_left_map_coords).m 
    y_dist_max = GD(bottom_right_map_coords,bottom_left_map_coords).m  
     
    [long_dist,lat_dist]  = np.meshgrid(np.linspace(0,y_dist_max,N_long),np.linspace(0,x_dist_max,N_lat))
    [long_deg,lat_deg]    = np.meshgrid(np.linspace(np.min(Long),np.max(Long),N_long),np.linspace(np.min(Lat),np.max(Lat),N_lat)) 
    elevation             = griddata((Lat,Long), Elev, (lat_deg, long_deg), method='linear')     
     
    return long_dist,lat_dist,long_deg,lat_deg,elevation 
  
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
