# Missions.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from RCAIDE.Framework.Core import Units  

#   Define the Mission
def stick_fixed_stability_setup(analyses,vehicle,cruise_velocity,cruise_altitude,angle_of_attack = 0 ,bank_angle= 0): 
    missions                     =RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier         = 1.0 # this multiplier is used to compute V_max from V_nominal
    missions.stick_fixed_cruise  = base_mission_setup(analyses.stick_fixed_cruise,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle) 
 
    return missions   

def elevator_sizing_setup(analyses,vehicle,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle): 
    missions =RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier                = 1.4 # this multiplier is used to compute V_max from V_nominal
    missions.elevator_sizing_pull_up    = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle)   

    max_speed_multiplier                = 1.0
    missions.elevator_sizing_push_over  = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle)   
    return missions   

def aileron_rudder_sizing_setup(analyses,vehicle,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle): 
    missions = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier      = 1.0     
    missions.aileron_sizing   = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle)  
    max_speed_multiplier      = 1.4   # this multiplier is used to compute V_max from V_nominal   
    missions.turn_criteria    = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle) 
 
    return missions   
    
def flap_sizing_setup(analyses,vehicle,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle): 
    missions = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier     = 1.0      
    missions.flap_sizing     = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle)   
    return missions        
    

# ------------------------------------------------------------------
#   Initialize the Mission
# ------------------------------------------------------------------    
    
def base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,bank_angle):   
    '''
    This sets up the nominal cruise of the aircraft
    '''
     
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'
  
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    #   Cruise Segment: constant Speed, constant altitude 
    segment                           = Segments.Untrimmed.Untrimmed()
    segment.analyses.extend( analyses )   
    segment.tag                       = "cruise"
    segment.angle_of_attack           = angle_of_attack
    segment.bank_angle                = bank_angle
    segment.altitude                  = cruise_altitude
    segment.air_speed                 = cruise_velocity * max_speed_multiplier

    segment.flight_dynamics.force_x   = True    
    segment.flight_dynamics.force_z   = True    
    segment.flight_dynamics.force_y   = True     
    segment.flight_dynamics.moment_y  = True 
    segment.flight_dynamics.moment_x  = True
    segment.flight_dynamics.moment_z  = True
    
    mission.append_segment(segment)     
    
    return mission
