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

sys.path.append( os.path.join(local_path_1, 'Flight_Path_Functions'))
sys.path.append( os.path.join(local_path_1, 'Post_Processing_Functions'))
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Hexacopter'))
from Hexacopter               import vehicle_setup, configs_setup
from compute_route_distances  import compute_route_distances
from compute_terrain_points   import compute_terrain_points
from Aircraft_Noise_Emissions import read_flight_simulation_results 
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
    routes_filepath       = relative_path  +  '..' + separator + 'UAM_City_Routes.xlsx'
    topography_file       = relative_path +  '..' + separator +  'Topography' + separator + 'LA_Metropolitan_Area.txt'
    flight_data           = pd.read_excel(routes_filepath,sheet_name=['Los_Angeles'])
    LA_flight_data_total  = flight_data['Los_Angeles']
    

    # ----------------------------------------------------------------------------------------------------------------------
    #  BATCH SETTINGS 
    # ----------------------------------------------------------------------------------------------------------------------      
    number_of_batches = 1
    batch_number      = 0 # THIS HAS TO BE CHANGED ON THE SERVER BEFORE YOU RUN IT 
    
    n_sims_total   = len(LA_flight_data_total)
    n_sims_group   = int(np.ceil(n_sims_total / number_of_batches))
    start          = batch_number * n_sims_group
    end            = (batch_number + 1) * n_sims_group
    LA_flight_data = LA_flight_data_total[start:end]
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  SIMULATION SETTINGS 
    # ----------------------------------------------------------------------------------------------------------------------  
    aircraft_code                = 'HC' # CHANGE FOR EACH AIRCRAFT 
    city_code                    = 'LA' 
    cruise_altitude              = 1000*Units.feet
    high_speed_climb_distance    = 1736 # DONE as of 11/18. CHANGE FOR EACH AIRCRAFT  
    radius_Vert1                 = 3600* Units.feet # circular pattern radius around vertiport 1
    radius_Vert2                 = 3600* Units.feet # circular pattern radius around vertiport 2
    dep_heading                  = 200 * Units.degree # Heading [degrees] of the departure from vertiport 1
    app_heading                  = 90  * Units.degree# Heading [degrees] of the approach to vertiport 2 
    max_cruise_distance          = 58*Units.nmi #CHANGE FOR EACH AIRCRAFT
    number_of_cpts               = 10  
    

    # ----------------------------------------------------------------------------------------------------------------------
    #  RUN BASELINE MISSION TO GET NOISE 
    # ----------------------------------------------------------------------------------------------------------------------     
    noise_results = run_noise_mission(number_of_cpts ) # Run noise simulation
    

    # ----------------------------------------------------------------------------------------------------------------------
    #  USE BASE MISSION TO GET NOISE OF ALL OPERATIONS 
    # ----------------------------------------------------------------------------------------------------------------------     
    filename_list = [] 
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
        total_cruise_distance, path_heading, dep_sector, app_sector = compute_route_distances(x1, y1, x2, y2, radius_Vert1, radius_Vert2, dep_heading, app_heading,high_speed_climb_distance)
        
        vehicle  = vehicle_setup(redesign_rotors= False) 
        
        # Set up configs
        configs  = configs_setup(vehicle)
        
        # vehicle analyses
        analyses = unconverged_analyses_setup(configs, origin_coord,destination_coord)
        
        # mission analyses 
        mission  = unconverged_mission_setup(number_of_cpts, analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, total_cruise_distance,cruise_altitude)        
        missions = missions_setup(mission) 
         
        if (max_cruise_distance > total_cruise_distance):
            
            # evaluate mission, not that it purposely does not converge
            results  = missions.base_mission.evaluate()
             
            N_segs = len(results.segments) 
            N_cpts = results.segments[0].state.numerics.number_of_control_points  
            for seg in range(N_segs):
                for i in range(N_cpts):  
                    results.segments[seg].state.conditions.noise  = noise_results.segments[seg].state.conditions.noise          
             
            filename =  aircraft_code +'_mission' + '_' + city_code + '_' + origin_code + '_' +  destination_code  + '_' + str(int(np.ceil(round(cruise_altitude/Units.feet,0)))) + 'ft'   
            res      =  read_flight_simulation_results(results, noise_results,  origin_coord, destination_coord)
            
            # save results 
            save(res, filename + '.res')
            
            filename_list.append(filename)
                
            filename_list_name =  aircraft_code + '_' + city_code +  '_Single_Flights_Raw'
            F =  Data(filename_list=filename_list)
            save(F, filename_list_name + '.res')
              
    return


def run_noise_mission(number_of_cpts):           
    vehicle  = vehicle_setup(redesign_rotors = False)     
    # Set up configs
    configs  = configs_setup(vehicle)

    # vehicle analyses
    analyses = noise_analyses_setup(configs)

    # mission analyses
    mission  = noise_mission_setup(number_of_cpts, analyses)
    missions = missions_setup(mission) 
     
    results = missions.base_mission.evaluate() 
    return results

def noise_analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = noise_base_analysis(config)
        analyses[tag] = analysis

    return analyses

def unconverged_analyses_setup(configs, origin_coord,destination_coord):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = unconverged_base_analysis(config)
        analyses[tag] = analysis

    return analyses

# ------------------------------------------------------------------
# Base Analysis
# ------------------------------------------------------------------
def noise_base_analysis(vehicle, origin_coord=[[0, 0]],destination_coord=[[0, 0]] ,mic_x_res=1000, mic_y_res=1000):
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
    #  Stability Analysis'HC_mission_LA_ONT_BUR_1000ft'
    stability         = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method() 
    stability.vehicle = vehicle 
    analyses.append(stability)    

    #  Noise Analysis   
    noise = RCAIDE.Framework.Analyses.Noise.Frequency_Domain_Buildup()   
    noise.vehicle = vehicle
    noise.settings.mean_sea_level_altitude          = False         
    #noise.settings.aircraft_origin_coordinates      = origin_coord  
    #noise.settings.aircraft_destination_coordinates = destination_coord  
    #noise.settings.microphone_x_resolution          = mic_x_res       
    #noise.settings.microphone_y_resolution          = mic_y_res         
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

def unconverged_base_analysis(vehicle):   
    
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
def noise_mission_setup(number_of_cpts, analyses, radius_Vert1=3600*Units.ft, radius_Vert2=3600*Units.ft, dep_heading=200*Units.degrees, app_heading=90*Units.degrees, dep_sector=90*Units.degrees, app_sector=90*Units.degrees, path_heading = 100, level_cruise_distance=10*Units.miles,cruise_altitude=1000*Units.ft): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment()
    base_segment.state.numerics.number_of_control_points    = number_of_cpts
    
    
    # ------------------------------------------------------------------
    #   Mission Constants
    # ------------------------------------------------------------------
    
    pattern_speed    = 55  * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
    cruise_speed     = 75. * Units['mph'] #CHANGE FOR EACH AIRCRAFT
    transition_speed = 35. * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
            
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 

    segment                                                = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag                                            = "Vertical_Climb" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                                 = 0.0  * Units.ft 
    segment.altitude_end                                   = 50.  * Units.ft  
    segment.climb_rate                                     = 500. * Units['ft/min']   
    segment.initial_battery_state_of_charge                = 1.0
    segment.true_course                                    = dep_heading   # This is in radians)     
    
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
    segment.altitude                         = 50.  * Units.ft       
    segment.air_speed_start                  = 300. * Units['ft/min'] 
    segment.air_speed_end                    = transition_speed    
    segment.acceleration                     = 1.0
    segment.pitch_initial                    = 0. * Units.degrees
    segment.pitch_final                      = 0. * Units.degrees
    segment.true_course                      = dep_heading     
 
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True 
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------'HC_mission_LA_ONT_BUR_1000ft'
    #   First Climb Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Departure_Climb"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 551.181  * Units['ft/min']
    segment.air_speed_start                  = transition_speed 
    segment.air_speed_end                    = pattern_speed       
    segment.altitude_start                   = 50.0 * Units.ft  
    segment.altitude_end                     = 500.0 * Units.ft
    segment.true_course                      = dep_heading    
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True 
    mission.append_segment(segment)
    
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Circular departure pattern 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Cruise.Curved_Constant_Radius_Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Departure_Pattern_Curve"     
    segment.analyses.extend(analyses.climb) 
    segment.altitude    = 500.0 * Units.ft  
    segment.air_speed   = pattern_speed       
    segment.turn_radius = radius_Vert1  
    segment.true_course = dep_heading + (90 * Units.degree)
    segment.turn_angle  = dep_sector
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                                             = True    
    segment.flight_dynamics.force_z                                             = True    
    segment.flight_dynamics.force_y                                             = True     
    segment.flight_dynamics.moment_y                                            = True 
    segment.flight_dynamics.moment_x                                            = True
    segment.flight_dynamics.moment_z                                            = True 

    # define flight controls              
    segment.assigned_control_variables.throttle.active                          = True
     
    segment.assigned_control_variables.throttle.assigned_propulsors             = [['rotor_propulsor_1','rotor_propulsor_2'],['rotor_propulsor_3'],
                                                                                    ['rotor_propulsor_4','rotor_propulsor_5'],['rotor_propulsor_6']]  
    segment.assigned_control_variables.body_angle.active                        = True     
    segment.assigned_control_variables.bank_angle.active                        = True    
    segment.assigned_control_variables.bank_angle.initial_guess_values          = [[20.0 * Units.degree]]
    
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Cruise_Climb"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 500. * Units['ft/min']
    segment.air_speed_start                  = pattern_speed
    segment.air_speed_end                    = cruise_speed      
    segment.altitude_start                   = 500.0 * Units.ft     
    segment.altitude_end                     = cruise_altitude
    segment.true_course                      = path_heading    
    
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
    segment.altitude                         = cruise_altitude      
    segment.air_speed                        = cruise_speed    
    segment.distance                         = level_cruise_distance
    segment.true_course                      = path_heading

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
    #   Cruise Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Cruise_Descent"  
    segment.analyses.extend(analyses.forward_flight)
    segment.climb_rate                       = -500. * Units['ft/min']
    segment.air_speed_start                  = cruise_speed      
    segment.air_speed_end                    = pattern_speed     
    segment.altitude_start                   = cruise_altitude
    segment.altitude_end                     = 500.0 * Units.ft
    segment.true_course                      = path_heading    

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True     
    mission.append_segment(segment)     
        
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Circular approach pattern 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Cruise.Curved_Constant_Radius_Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Approach_Pattern_Curve"     
    segment.analyses.extend(analyses.climb) 
    segment.altitude       = 500.0 * Units.ft  
    segment.air_speed      = pattern_speed     
    segment.turn_radius    = radius_Vert2  
    segment.true_course    = path_heading - (90 *Units.degrees)  
    segment.turn_angle     = app_sector  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                                             = True    
    segment.flight_dynamics.force_z                                             = True    
    segment.flight_dynamics.force_y                                             = True     
    segment.flight_dynamics.moment_y                                            = True 
    segment.flight_dynamics.moment_x                                            = True
    segment.flight_dynamics.moment_z                                            = True 

    # define flight controls              
    segment.assigned_control_variables.throttle.active                          = True           
    segment.assigned_control_variables.throttle.assigned_propulsors             = [['rotor_propulsor_1','rotor_propulsor_2'],['rotor_propulsor_3'],
                                                                                    ['rotor_propulsor_4','rotor_propulsor_5'],['rotor_propulsor_6']]  
    segment.assigned_control_variables.body_angle.active                        = True     
    segment.assigned_control_variables.bank_angle.active                        = True    
    segment.assigned_control_variables.bank_angle.initial_guess_values          = [[20.0 * Units.degree]]
    
    mission.append_segment(segment) 
 
     
       # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Approach_Descent"  
    segment.analyses.extend(analyses.forward_flight)
    segment.climb_rate                       = -551.181  * Units['ft/min']
    segment.air_speed_start                  = pattern_speed      
    segment.air_speed_end                    = transition_speed      
    segment.altitude_start                   = 500.0 * Units.ft  
    segment.altitude_end                     = 50.0 * Units.ft
    segment.true_course                      = app_heading

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
    segment.tag                       = "Final_Approach_Transition" 
    segment.analyses.extend( analyses.descent_transition) 
    segment.altitude                  = 50.  * Units.ft 
    segment.air_speed_start           = transition_speed  
    segment.air_speed_end             = 300. * Units['ft/min']
    segment.acceleration              = -1.0 
    segment.pitch_initial             = 1. * Units.degrees
    segment.pitch_final               = 2. * Units.degrees
    segment.true_course               = app_heading 

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
    segment.altitude_start            = 50.0  * Units.ft  
    segment.altitude_end              = 0.  * Units.ft  
    segment.descent_rate              = 300. * Units['ft/min']
    segment.true_course               = app_heading 
    

    # define flight dynamics to model  
    segment.flight_dynamics.force_z                        = True 

    # define flight controls  
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]      
    mission.append_segment(segment)      
    
    return mission 
# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def unconverged_mission_setup(number_of_cpts,analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, level_cruise_distance,cruise_altitude): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment()
    base_segment.state.numerics.number_of_control_points    = number_of_cpts
    
    
    # ------------------------------------------------------------------
    #   Mission Constants
    # ------------------------------------------------------------------
    
    pattern_speed    = 55  * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
    cruise_speed     = 75. * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
    transition_speed = 35. * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
            
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 

    segment                                                = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag                                            = "Vertical_Climb" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                                 = 0.0  * Units.ft
    del segment.process.converge
    segment.altitude_end                                   = 50.  * Units.ft  
    segment.climb_rate                                     = 500. * Units['ft/min']   
    segment.initial_battery_state_of_charge                = 1.0
    segment.true_course                                    = dep_heading   # This is in radians)     
    
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
    segment.altitude                         = 50.  * Units.ft
    del segment.process.converge
    segment.air_speed_start                  = 300. * Units['ft/min'] 
    segment.air_speed_end                    = transition_speed    
    segment.acceleration                     = 1.0
    segment.pitch_initial                    = 0. * Units.degrees
    segment.pitch_final                      = 0. * Units.degrees
    segment.true_course                      = dep_heading     
 
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
    #   First Climb Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Departure_Climb"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 551.181  * Units['ft/min']
    segment.air_speed_start                  = transition_speed
    del segment.process.converge
    segment.air_speed_end                    = pattern_speed       
    segment.altitude_start                   = 50.0 * Units.ft  
    segment.altitude_end                     = 500.0 * Units.ft
    segment.true_course                      = dep_heading      
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True 
    mission.append_segment(segment)
    
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Circular departure pattern 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Cruise.Curved_Constant_Radius_Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Departure_Pattern_Curve"     
    segment.analyses.extend(analyses.climb) 
    segment.altitude    = 500.0 * Units.ft
    del segment.process.converge
    segment.air_speed   = pattern_speed       
    segment.turn_radius = radius_Vert1  
    segment.true_course = dep_heading + (90 * Units.degree)
    segment.turn_angle  = dep_sector
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                                             = True    
    segment.flight_dynamics.force_z                                             = True    
    segment.flight_dynamics.force_y                                             = True     
    segment.flight_dynamics.moment_y                                            = True 
    segment.flight_dynamics.moment_x                                            = True
    segment.flight_dynamics.moment_z                                            = True 

    # define flight controls              
    segment.assigned_control_variables.throttle.active                          = True
     
    segment.assigned_control_variables.throttle.assigned_propulsors             = [['rotor_propulsor_1','rotor_propulsor_2'],['rotor_propulsor_3'],
                                                                                    ['rotor_propulsor_4','rotor_propulsor_5'],['rotor_propulsor_6']]  
    segment.assigned_control_variables.body_angle.active                        = True     
    segment.assigned_control_variables.bank_angle.active                        = True    
    segment.assigned_control_variables.bank_angle.initial_guess_values          = [[20.0 * Units.degree]]
    
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Cruise_Climb"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 500. * Units['ft/min']
    segment.air_speed_start                  = pattern_speed
    del segment.process.converge
    segment.air_speed_end                    = cruise_speed      
    segment.altitude_start                   = 500.0 * Units.ft     
    segment.altitude_end                     = cruise_altitude
    segment.true_course                      = path_heading    
    
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
    del segment.process.converge
    segment.altitude                         = cruise_altitude      
    segment.air_speed                        = cruise_speed    
    segment.distance                         = level_cruise_distance
    segment.true_course                      = path_heading

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
    #   Cruise Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Cruise_Descent"  
    segment.analyses.extend(analyses.forward_flight)
    del segment.process.converge
    segment.climb_rate                       = -500. * Units['ft/min']
    segment.air_speed_start                  = cruise_speed      
    segment.air_speed_end                    = pattern_speed     
    segment.altitude_start                   = cruise_altitude
    segment.altitude_end                     = 500.0 * Units.ft
    segment.true_course                      = path_heading    

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True     
    mission.append_segment(segment)     
        
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Circular approach pattern 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Cruise.Curved_Constant_Radius_Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Approach_Pattern_Curve"     
    segment.analyses.extend(analyses.climb)
    del segment.process.converge
    segment.altitude       = 500.0 * Units.ft  
    segment.air_speed      = pattern_speed     
    segment.turn_radius    = radius_Vert2  
    segment.true_course    = path_heading - (90 *Units.degrees)     # CHECK UNITS
    segment.turn_angle     = app_sector  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                                             = True    
    segment.flight_dynamics.force_z                                             = True    
    segment.flight_dynamics.force_y                                             = True     
    segment.flight_dynamics.moment_y                                            = True 
    segment.flight_dynamics.moment_x                                            = True
    segment.flight_dynamics.moment_z                                            = True 

    # define flight controls              
    segment.assigned_control_variables.throttle.active                          = True           
    segment.assigned_control_variables.throttle.assigned_propulsors             = [['rotor_propulsor_1','rotor_propulsor_2'],['rotor_propulsor_3'],
                                                                                    ['rotor_propulsor_4','rotor_propulsor_5'],['rotor_propulsor_6']]  
    segment.assigned_control_variables.body_angle.active                        = True     
    segment.assigned_control_variables.bank_angle.active                        = True    
    segment.assigned_control_variables.bank_angle.initial_guess_values          = [[20.0 * Units.degree]]
    
    mission.append_segment(segment) 
 
     
       # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Approach_Descent"  
    segment.analyses.extend(analyses.forward_flight)
    del segment.process.converge
    segment.climb_rate                       = -551.181  * Units['ft/min']
    segment.air_speed_start                  = pattern_speed      
    segment.air_speed_end                    = transition_speed      
    segment.altitude_start                   = 500.0 * Units.ft 
    segment.altitude_end                     = 50.0 * Units.ft
    segment.true_course                      = app_heading

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
    segment.tag                       = "Final_Approach_Transition" 
    segment.analyses.extend( analyses.descent_transition)
    del segment.process.converge
    segment.altitude                  = 50.  * Units.ft 
    segment.air_speed_start           = transition_speed  
    segment.air_speed_end             = 300. * Units['ft/min']
    segment.acceleration              = -1.0 
    segment.pitch_initial             = 1. * Units.degrees
    segment.pitch_final               = 2. * Units.degrees
    segment.true_course               = app_heading 

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
    del segment.process.converge
    segment.altitude_start            = 50.0  * Units.ft  
    segment.altitude_end              = 0.  * Units.ft  
    segment.descent_rate              = 300. * Units['ft/min']
    segment.true_course               = app_heading 
    

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