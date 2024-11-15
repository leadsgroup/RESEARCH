#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE 
from RCAIDE.Framework.Core import Units, Data   
from  RCAIDE.Framework.Analyses.Geodesics.Geodesics import Calculate_Distance  
from RCAIDE.Library.Plots import *       
from RCAIDE import  load 
from RCAIDE import  save

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
    
    for i in  range(len(LA_flight_data)):
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
        x1 = Calculate_Distance(x0_coord,bottom_left_map_coords) * Units.kilometers # Correct
        y1 = Calculate_Distance(y0_coord,bottom_left_map_coords) * Units.kilometers # Correct
        x2 = Calculate_Distance([bottom_left_map_coords[0], destination_coord[1]],bottom_left_map_coords) * Units.kilometers # Double check
        y2 = Calculate_Distance([destination_coord[0], bottom_left_map_coords[1]],bottom_left_map_coords) * Units.kilometers # Double check
        
        # -------------------------------------
        #   Calculate Distance
        # -------------------------------------
        total_cruise_distance, path_heading, dep_sector, app_sector = compute_route_distances(x1, y1, x2, y2, radius_Vert1, radius_Vert2, dep_heading, app_heading,high_speed_climb_distance,high_speed_descent_distance)
        
        vehicle  = vehicle_setup(redesign_rotors= False) 
        
        # Set up configs
        configs  = configs_setup(vehicle)
        
        # vehicle analyses
        analyses = analyses_setup(configs, origin_coord,destination_coord ,mic_x_res, mic_y_res ,noise_timesteps ,mic_stencil)
        
        # mission analyses 
        mission = mission_setup(analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, total_cruise_distance,cruise_altitude)        
        missions = missions_setup(mission) 
         
        if max_cruise_distance > total_cruise_distance:
            results = missions.base_mission.evaluate() 
            
            # post process noise 
            noise_data   = post_process_noise_data(results,
                                                   flight_times = operation_flight_times,  
                                                   time_period  = operation_period,
                                                   evalaute_noise_metrics = False)  
          
            # save data
            filename =  aircraft_code + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + cruise_altitude    # Aircraft_City_Frequency_Origin_Destination_Altitude
            save(noise_data, filename + '.res') 
        
    return  

def analyses_setup(configs, origin_coord,destination_coord ,mic_x_res, mic_y_res ,noise_timesteps ,mic_stencil):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config, origin_coord,destination_coord ,mic_x_res, mic_y_res ,noise_timesteps ,mic_stencil)
        analyses[tag] = analysis

    return analyses

# ------------------------------------------------------------------
# Base Analysis
# ------------------------------------------------------------------
def base_analysis(vehicle, origin_coord,destination_coord ,mic_x_res, mic_y_res ,noise_timesteps ,mic_stencil):
    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator 
    topography_file = relative_path +  '..' + separator +  'Topography' + separator + 'LA_Metropolitan_Area.txt'    
    
    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle() 
    
    # ------------------------------------------------------------------
    #  Weights
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics         = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle = vehicle 
    analyses.append(aerodynamics)
     
    # ------------------------------------------------------------------
    #  Stability Analysis
    stability         = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method() 
    stability.vehicle = vehicle 
    analyses.append(stability)    

    #  Noise Analysis   
    noise = RCAIDE.Framework.Analyses.Noise.Frequency_Domain_Buildup()   
    noise.vehicle = vehicle
    noise.settings.mean_sea_level_altitude          = False         
    noise.settings.aircraft_origin_coordinates      = origin_coord  
    noise.settings.aircraft_destination_coordinates = destination_coord  
    noise.settings.microphone_x_resolution          = mic_x_res       
    noise.settings.microphone_y_resolution          = mic_y_res        
    noise.settings.noise_times_steps                = noise_timesteps 
    noise.settings.number_of_microphone_in_stencil  = mic_stencil     
    noise.settings.topography_file                  = topography_file
    analyses.append(noise)
 
    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses    



# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def mission_setup(analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, total_cruise_distance,cruise_altitude): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment()
    
           
            
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 

    segment                                                = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag                                            = "Vertical_Climb" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                                 = 0.0  * Units.ft 
    segment.altitude_end                                   = 200.  * Units.ft  
    segment.climb_rate                                     = 300. * Units['ft/min']   
    segment.initial_battery_state_of_charge                = 1.0
    
    # define flight dynamics to model  
    segment.flight_dynamics.force_z                        = True 

    # define flight controls  
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #  First Transition Segment
    # ------------------------------------------------------------------  

    segment                                  = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag                              = "Vertical_Transition"  
    segment.analyses.extend( analyses.vertical_transition) 
    segment.altitude                         = 200.  * Units.ft       
    segment.air_speed_start                  = 300. * Units['ft/min'] 
    segment.air_speed_end                    = 35 * Units['mph']    
    segment.acceleration                     = 0.5
    segment.pitch_initial                    = 0. * Units.degrees
    segment.pitch_final                      = 0. * Units.degrees
 
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True 
    
    mission.append_segment(segment)

    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Climb_1"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 600. * Units['ft/min']
    segment.air_speed_start                  = 35.   * Units['mph']
    segment.air_speed_end                    = 55.  * Units['mph']       
    segment.altitude_start                   = 200.0 * Units.ft  
    segment.altitude_end                     = 500.0 * Units.ft       
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True 
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Climb_2"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 500. * Units['ft/min']
    segment.air_speed_start                  = 55.   * Units['mph']
    segment.air_speed_end                    = 75.  * Units['mph']       
    segment.altitude_start                   = 500.0 * Units.ft     
    segment.altitude_end                     = 2500.0 * Units.ft
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True     
    mission.append_segment(segment)                

    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                              = "Cruise"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude                         = 2500.0 * Units.ft      
    segment.air_speed                        = 75. * Units['mph']      
    segment.distance                         = 10*Units.nmi

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True     
    mission.append_segment(segment)      
    
                
    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Descent_1"  
    segment.analyses.extend(analyses.forward_flight)
    segment.climb_rate                       = -200. * Units['ft/min']
    segment.air_speed_start                  = 75. * Units['mph']      
    segment.air_speed_end                    = 35. * Units['mph']      
    segment.altitude_start                   = 2500.0 * Units.ft 
    segment.altitude_end                     = 200.0 * Units.ft

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True     
    mission.append_segment(segment)     
        
     
    # ------------------------------------------------------------------
    #  Third Transition Segment
    # ------------------------------------------------------------------

    segment                           = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag                       = "Decent_Transition" 
    segment.analyses.extend( analyses.descent_transition) 
    segment.altitude                  = 200.  * Units.ft 
    segment.air_speed_start           = 35.  * Units['mph'] 
    segment.air_speed_end             = 300. * Units['ft/min']
    segment.acceleration              = -0.5307 
    segment.pitch_initial             = 1. * Units.degrees
    segment.pitch_final               = 2. * Units.degrees

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True     
    mission.append_segment(segment)

    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                           = Segments.Vertical_Flight.Descent(base_segment)
    segment.tag                       = "Vertical_Descent"  
    segment.analyses.extend( analyses.vertical_flight) 
    segment.altitude_start            = 200.0  * Units.ft  
    segment.altitude_end              = 0.  * Units.ft  
    segment.descent_rate              = 300. * Units['ft/min']

    # define flight dynamics to model  
    segment.flight_dynamics.force_z                        = True 

    # define flight controls  
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]      
    mission.append_segment(segment)
             
  
    return mission 
 

def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions 
 
if __name__ == '__main__': 
    main()     