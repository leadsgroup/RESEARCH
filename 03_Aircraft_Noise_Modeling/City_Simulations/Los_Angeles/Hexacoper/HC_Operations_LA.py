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
import  pandas as pd 

sys.path.append(os.path.join( os.path.split(os.path.split(sys.path[0])[0])[0], 'Aircraft' + os.sep + 'Hexacopter'))
from Hexacopter    import vehicle_setup, configs_setup
from compute_route_distances import  compute_route_distances
# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def main():           
         

    separator            = os.path.sep  
    technology_filename  = 'Data' + separator + 'Technology' + separator + 'Technology_Data.xlsx'
    flight_data          = pd.read_excel(technology_filename,sheet_name=['Los_Angeles'])  
    
    
    operation_flight_times = np.array(['06:00:00','06:05:00','06:10:00','06:15:00',
                             '06:20:00','06:25:00','06:30:00','06:35:00',
                             '06:40:00','06:45:00','06:50:00','06:55:00',
                             '07:00:00','07:05:00','07:10:00','07:15:00',
                             '07:20:00','07:25:00','07:30:00','07:35:00',
                             '07:40:00','07:45:00','07:50:00','07:55:00', 
                             '08:00:00'])
    operation_period  = ['06:00:00','09:00:00']
        
    # Read in Excel Sheet
    flight_data =
    
    Hexacopter_Max_Cruise_Distance =  30 * Units.mile
    
    for i in  range(flight_data):
        # Extract Data
        departure_vertiport_code =
        origin_vertiport_code    =
        departure_coords         =  
        origin_coords            =    
        distance                 = 
         
       
        # -------------------------------------
        #   Inputs
        # -------------------------------------
        vert_Loc1 = np.array([33.932868, -118.380858]) # Lat, Long coordinates of vertiport 1
        vert_Loc2 = np.array([34.057759, -118.221603]) # Lat, Long coordinates of vertiport 2
        radius_Vert1 = 1 * Units.km # circular pattern radius around vertiport 1
        radius_Vert2 = 1 * Units.km # circular pattern radius around vertiport 2
        dep_heading = 200 # Heading [degrees] of the departure from vertiport 1
        app_heading = 90 # Heading [degrees] of the approach to vertiport 2
        
        # -------------------------------------
        #   Lat-lon to X-Y. Assume that vertiport 1 is at 0,0 and then calcualte vertiport two lcoation. We'll calcualte everything in this frame, find the distnaces and then the program when it converts the mission profile back will handle it on that side. 
        # -------------------------------------
        x1 = 0 #Calculate_Distance(x0_coord,bottom_left_map_coords) * Units.kilometers
        y1 = 0 #Calculate_Distance(y0_coord,bottom_left_map_coords) * Units.kilometers
        x2 = Calculate_Distance([vert_Loc1[0], vert_Loc2[1]],vert_Loc1) * Units.kilometers # Double check
        y2 = Calculate_Distance([vert_Loc1[0], vert_Loc2[1]],vert_Loc2) * Units.kilometers # Double check
        # -------------------------------------
        #   Calculate Distance
        # -------------------------------------
        total_cruise_distance, path_heading, dep_sector, app_sector = compute_route_distances(x1, y1, x2, y2, radius_Vert1, radius_Vert2, dep_heading, app_heading)
       
        
        vehicle  = vehicle_setup() 
        
        # Set up configs
        configs  = configs_setup(vehicle)
        
        # vehicle analyses
        analyses = analyses_setup(configs)
        
        # mission analyses 
        mission = mission_setup(analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, total_cruise_distance)        
        missions = missions_setup(mission) 
         
        if Hexacopter_Max_Cruise_Distance < distance:
            results = missions.base_mission.evaluate() 
            
            # post process noise 
            noise_data   = post_process_noise_data(results,
                                                   flight_times = operation_flight_times,  
                                                   time_period  = operation_period,
                                                   evalaute_noise_metrics = False)  
          
            # save data
            filename = ''   # Aircraft_City_Frequency_Origin_Destination_Altitude
            save(noise_data, filename + '.res')
        
        
    return  

def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

# ------------------------------------------------------------------
# Base Analysis
# ------------------------------------------------------------------
def base_analysis(vehicle):

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
def mission_setup(analyses): 
    
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

    segment                                            = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag                                        = "Vertical_Climb" 
    segment.analyses.extend( analyses.vertical_flight) 
    segment.altitude_start                             = 0.0  * Units.ft 
    segment.altitude_end                               = 200.  * Units.ft  
    segment.climb_rate                                 = 300. * Units['ft/min']   
    segment.initial_battery_state_of_charge            = 1.0
    
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