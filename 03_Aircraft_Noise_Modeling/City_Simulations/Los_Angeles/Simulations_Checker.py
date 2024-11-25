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

local_path_1 =  os.path.split(os.path.split(sys.path[0])[0])[0] 
local_path_2 =  os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0] 

sys.path.append( os.path.join(local_path_1, 'Flight_Path_Functions'))
sys.path.append( os.path.join(local_path_1, 'Post_Processing_Functions'))

from compute_terrain_points   import compute_terrain_points
from compute_route_distances  import compute_route_distances
# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def main():           

    # ----------------------------------------------------------------------------------------------------------------------
    # FILE IMPORTS 
    # ----------------------------------------------------------------------------------------------------------------------            
    ospath                = os.path.abspath(__file__)
    separator             = os.path.sep
    relative_path         = os.path.dirname(ospath) + separator 
    routes_filepath       = relative_path + 'UAM_City_Routes.xlsx'
    topography_file       = relative_path +  'Topography' + separator + 'LA_Metropolitan_Area.txt'
    flight_data           = pd.read_excel(routes_filepath,sheet_name=['Los_Angeles'])
    LA_flight_data        = flight_data['Los_Angeles']
     
    # ----------------------------------------------------------------------------------------------------------------------
    #  SIMULATION SETTINGS 
    # ----------------------------------------------------------------------------------------------------------------------  
    high_speed_climb_distance    = 1736 # DONE as of 11/18. CHANGE FOR EACH AIRCRAFT  
    radius_Vert1                 = 3600* Units.feet # circular pattern radius around vertiport 1
    radius_Vert2                 = 3600* Units.feet # circular pattern radius around vertiport 2
    dep_heading                  = 200 * Units.degree # Heading [degrees] of the departure from vertiport 1
    app_heading                  = 90  * Units.degree# Heading [degrees] of the approach to vertiport 2   
    
    for i in range(len(LA_flight_data)):
        # Extract Data
        origin_code       = LA_flight_data['Origin Code'][i]   
        destination_code  = LA_flight_data['Destination Code'][i] 
        origin_coord      = [LA_flight_data['Origin Latitude'][i], LA_flight_data['Origin Longitude'][i]]
        destination_coord = [LA_flight_data['Destination Latitude'][i], LA_flight_data['Destination Longitude'][i]]
        

        terrain_data =  compute_terrain_points(topography_file, 
                               number_of_latitudinal_points  = 100,
                               number_of_longitudinal_points = 100) 
    
        y0_coord               = np.array([origin_coord[0], terrain_data['bottom_left_map_coordinates'][1]]) # verify this
        bottom_left_map_coords = terrain_data['bottom_left_map_coordinates']
        x0_coord               = np.array([terrain_data['bottom_left_map_coordinates'][0], origin_coord[1]])
        
        # -------------------------------------
        #   Lat-lon to X-Y. Assume that vertiport 1 is at 0,0 and then calcualte vertiport two lcoation. We'll calcualte everything in this frame, find the distnaces and then the program when it converts the mission profile back will handle it on that side. 
        # -------------------------------------
        y1 = Calculate_Distance(x0_coord,bottom_left_map_coords) * Units.kilometers # Correct
        x1 = Calculate_Distance(y0_coord,bottom_left_map_coords) * Units.kilometers # Correct
        y2 = Calculate_Distance([bottom_left_map_coords[0], destination_coord[1]],bottom_left_map_coords) * Units.kilometers # Double check
        x2 = Calculate_Distance([destination_coord[0], bottom_left_map_coords[1]],bottom_left_map_coords) * Units.kilometers # Double check
        
        # -------------------------------------
        #   Calculate Distance
        # -------------------------------------
        _, _, _, _ = compute_route_distances(x1, y1, x2, y2, radius_Vert1, radius_Vert2, dep_heading, app_heading,high_speed_climb_distance)
         
      
    return
 
if __name__ == '__main__': 
    main()     