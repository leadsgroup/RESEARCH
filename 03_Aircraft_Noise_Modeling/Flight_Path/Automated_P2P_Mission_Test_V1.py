# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
from  RCAIDE.Framework.Analyses.Geodesics.Geodesics import Calculate_Distance
from RCAIDE.Framework.Core import Units
import RCAIDE
from RCAIDE.Library.Methods.Performance.estimate_stall_speed                   import estimate_stall_speed

import  numpy as  np
import  math
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def  Create_P2P_Mission(analyses, vehicle):
    
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
    total_cruise_distance, path_heading, dep_sector, app_sector = Route_Distances(x1, y1, x2, y2, radius_Vert1, radius_Vert2, dep_heading, app_heading)
  
    # -------------------------------------
    #   Create Mission
    # -------------------------------------
    mission = mission_setup(analyses,vehicle, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, total_cruise_distance)
    
    return(mission)

    
def mission_setup(analyses,vehicle, radius_Vert1, radius_Vert2, departure_heading, approach_heading, dep_sector, app_sector, path_heading, total_cruise_distance): 
    
    # ------------------------------------------------------------------
    #   Mission Constants
    # ------------------------------------------------------------------    
    transition_altitude = 20 *  Units.ft
    dep_pattern_altitude =  500 *  Units.ft
    app_pattern_altitude = 500 * Units.ft
    cruise_altitude = 1000 * Units.ft
    cruise_speed = 130 * Units.kts
    pattern_speed = 90 * Units.kts
    transition_speed = 10 * Units.kts

    departure_climb_angle = np.arctan((dep_pattern_altitude-transition_altitude)/radius_Vert1)
    approach_descent_angle  = np.arctan((app_pattern_altitude-transition_altitude)/radius_Vert2)
    climb_rate            = 720.  * Units['ft/min'] # https://www.icas.org/ICAS_ARCHIVE/ICAS2022/data/papers/ICAS2022_0408_paper.pdf
    descent_rate          = -720. * Units['ft/min'] # https://www.icas.org/ICAS_ARCHIVE/ICAS2022/data/papers/ICAS2022_0408_paper.pdf

    descent_deceleration = np.arctan((app_pattern_altitude-transition_altitude)/radius_Vert2) *(pattern_speed - transition_speed) * abs(descent_rate) / (app_pattern_altitude-transition_altitude)
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'baseline_mission' 
    
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    # base segment           
    base_segment  = Segments.Segment()   
    ''' NOT USED?
    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)      
    '''
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Vertical Climb 
    #------------------------------------------------------------------------------------------------------------------------------------  
    segment     = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )  
    segment.altitude_start                                = 0.0  * Units.ft  
    segment.altitude_end                                  = transition_altitude 
    segment.initial_battery_state_of_charge               = 1.0 
    segment.climb_rate                                    = 500. * Units['ft/min']   
            
    # define flight dynamics to model  
    segment.flight_dynamics.force_z                       = True     
    
    '''Change the following absed on the aircraft type'''
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                              'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']] 
       
    mission.append_segment(segment)
      
    #------------------------------------------------------------------------------------------------------------------------------------  
    # High-Speed Climbing Transition
    #------------------------------------------------------------------------------------------------------------------------------------  
    segment                                               = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                           = "High_Speed_Climbing_Transition" 
    segment.analyses.extend( analyses.transition_flight)    
    segment.altitude_start                                = transition_altitude  
    segment.altitude_end                                  = dep_pattern_altitude
    segment.climb_angle                                   = departure_climb_angle
    segment.air_speed                                     = 20    * Units.mph
    segment.acceleration                                  = 0.25  * Units['m/s/s'] 
    segment.pitch_initial                                 = 0.    * Units.degrees 
    segment.pitch_final                                   = 7.    * Units.degrees
    segment.true_course                                   = departure_heading

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2'],
                                                             ['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                            'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']]
    mission.append_segment(segment) 

    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Circular departure pattern 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Cruise.Curved_Constant_Radius_Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Departure_pattern_curve"       
    
    segment.altitude    = dep_pattern_altitude
    segment.air_speed   = pattern_speed   
    segment.turn_radius = radius_Vert1
    segment.true_course = departure_heading + 90 * Units.deg
    segment.turn_angle  = dep_sector
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2'],
                                                             ['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                            'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']]
    mission.append_segment(segment) 
 
    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Cruise climb Cstant_Speed_Linear_Altitude. Add heading = 
    #------------------------------------------------------------------------------------------------------------------------------------  
    segment                                               = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                           = "Cruise_Climb"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                                = dep_pattern_altitude  
    segment.altitude_end                                  = cruise_altitude 
    segment.climb_rate                                    = climb_rate
    segment.air_speed_end                                 = cruise_speed
    segment.true_course                                   = path_heading
            
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    '''Change based on the aircraft'''
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2']] 
    segment.assigned_control_variables.body_angle.active             = True                
                
    mission.append_segment(segment)  

    # Calculate level cruise distance.
    time_to_climb = (cruise_altitude - dep_pattern_altitude) / climb_rate
    climb_distance = (math.sqrt(pattern_speed**2-climb_rate**2) + math.sqrt(cruise_speed**2-climb_rate**2)) /2 * time_to_climb
    
    time_to_descent = abs((cruise_altitude - app_pattern_altitude) / descent_rate)
    descent_distance = (math.sqrt(pattern_speed**2-descent_rate**2) + math.sqrt(cruise_speed**2-descent_rate**2)) /2 * time_to_descent
    level_cruise_distance = total_cruise_distance - climb_distance - descent_distance
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Cruise 
    #------------------------------------------------------------------------------------------------------------------------------------  
    segment                                               = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Cruise"  
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                      = cruise_altitude  
    segment.air_speed                                     = cruise_speed
    segment.distance                                      = level_cruise_distance
    segment.true_course                                   = path_heading
            
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2']] 
    segment.assigned_control_variables.body_angle.active             = True                
         
    mission.append_segment(segment)   
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Cruise descent
    #------------------------------------------------------------------------------------------------------------------------------------   
    segment                                               = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                           = "Descent"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                                = cruise_altitude
    segment.altitude_end                                  = app_pattern_altitude
    segment.climb_rate                                    = descent_rate
    segment.air_speed_start                               = cruise_speed 
    segment.air_speed_end                                 = pattern_speed
    segment.true_course                                   = path_heading
            
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2']] 
    segment.assigned_control_variables.body_angle.active             = True                
       
    mission.append_segment(segment)  
    
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Circular approach pattern 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Cruise.Curved_Constant_Radius_Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Approach_pattern_curve"       
    
    segment.altitude    = app_pattern_altitude
    segment.air_speed   = pattern_speed   
    segment.turn_radius = radius_Vert2
    segment.true_course = path_heading - 90 * Units.deg
    segment.turn_angle  = app_sector
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2'],
                                                             ['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                            'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']]
    mission.append_segment(segment) 
    
            
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # High-Speed Descending Transition
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                           = "High_Speed_Descending_Transition"  
    segment.analyses.extend( analyses.transition_flight )
    segment.air_speed_start                               = pattern_speed
    segment.altitude_start                                = app_pattern_altitude   
    segment.altitude_end                                  = transition_altitude  
    segment.climb_angle                                   = approach_descent_angle
    segment.acceleration                                  = descent_deceleration    
    segment.pitch_initial                                 = 4.3  * Units.degrees  
    segment.pitch_final                                   = 7. * Units.degrees
    segment.true_course                                   = approach_heading

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2'],
                                                             ['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                            'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8'] ]
        
    mission.append_segment(segment)  

    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Vertical Descent 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    segment                                               = Segments.Vertical_Flight.Descent(base_segment)
    segment.tag                                           = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                                = transition_altitude   
    segment.altitude_end                                  = 0.   * Units.ft  
    segment.descent_rate                                  = 300. * Units['ft/min']  
    
    # define flight dynamics to model  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                              'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']] 
            
    mission.append_segment(segment)  
  
    # DELETED CHARGE SEGMENT
    return mission 

def Route_Distances(x1, y1, x2, y2, radius_Vert1, radius_Vert2, dep_heading, app_heading):
    
    dep_distance = radius_Vert1
    app_distance = radius_Vert2
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Path heading:
    # ---------------------------------------------------------------------------------------------------------------------- 
    path_heading = (np.arctan2((x2 - x1), (y2 - y1)) *180 /np.pi) % 360
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Departure sector: follows a right-hand pattern
    # ---------------------------------------------------------------------------------------------------------------------- 
    if (path_heading > dep_heading):
        dep_sector = path_heading - dep_heading # negative means counter clockwise, positive means clockwise
    else:
        dep_sector = 360 + path_heading - dep_heading
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Departure pattern distance:
    # ---------------------------------------------------------------------------------------------------------------------- 
    dep_pattern_distance = abs(dep_sector*np.pi/180) *radius_Vert1
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Path distance:
    # ---------------------------------------------------------------------------------------------------------------------- 
    path_distance = np.sqrt((y1-y2)**2 + (x1-x2)**2) - radius_Vert1 - radius_Vert2
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Path approach heading: The angle measures on the circle that is the traffic apttern are 180 degrees offset from the heading 
    # ---------------------------------------------------------------------------------------------------------------------- 
    path_heading_angle = (path_heading + 180) % 360 # in
    app_angle = (app_heading + 180) % 360 #out 
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Approach sector: assumes a right hand pattern 
    # ---------------------------------------------------------------------------------------------------------------------- 
    if path_heading_angle > app_angle:
        app_sector = 360 - (path_heading_angle - app_angle)
    else :
        app_sector = app_angle - path_heading_angle
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Approach pattern distance:
    # ---------------------------------------------------------------------------------------------------------------------- 
    app_pattern_distance = abs(app_sector*np.pi/180) * radius_Vert2
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Total distance:
    # ---------------------------------------------------------------------------------------------------------------------- 
    total_distance = dep_distance + dep_pattern_distance + path_distance + app_pattern_distance + app_distance
    
    # ---------------------------------------------------------------------------------------------------------------------- 
    # Point-to-point distance:
    # ---------------------------------------------------------------------------------------------------------------------- 
    p2p_distance = np.sqrt((y1-y2)**2 + (x1-x2)**2)
    
    return path_distance, path_heading, dep_sector, app_sector

if __name__ == '__main__':       
    Create_P2P_Mission() 