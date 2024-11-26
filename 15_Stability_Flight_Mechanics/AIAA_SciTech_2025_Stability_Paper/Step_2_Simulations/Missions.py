# Missions.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE

def stick_fixed_stability_setup(analyses,cruise_velocity,cruise_altitude,angle_of_attack = 0 ,sideslip_angle=0,bank_angle= 0, roll_rate = 0,pitch_rate = 0,yaw_rate = 0): 
    missions                     = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier         = 1.0 
    missions.stick_fixed_cruise  = base_mission_setup(analyses.stick_fixed_cruise,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate) 
 
    return missions    
      
def base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate):   
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
    segment.sideslip_angle            = sideslip_angle
    segment.bank_angle                = bank_angle
    segment.altitude                  = cruise_altitude
    segment.air_speed                 = cruise_velocity * max_speed_multiplier
    segment.roll_rate                 = roll_rate 
    segment.pitch_rate                = pitch_rate  
    segment.yaw_rate                  = yaw_rate   

    segment.flight_dynamics.force_x   = True    
    segment.flight_dynamics.force_z   = True    
    segment.flight_dynamics.force_y   = True     
    segment.flight_dynamics.moment_y  = True 
    segment.flight_dynamics.moment_x  = True
    segment.flight_dynamics.moment_z  = True
    
    mission.append_segment(segment)     
    
    return mission
 
 