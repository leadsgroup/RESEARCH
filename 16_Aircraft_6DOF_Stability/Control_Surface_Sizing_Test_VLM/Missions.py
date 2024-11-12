# Missions.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from RCAIDE.Framework.Core import Units  

#   Define the Mission
def stick_fixed_stability_setup(analyses,vehicle,cruise_velocity,cruise_altitude): 
    missions                     =RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier         = 1.0 # this multiplier is used to compute V_max from V_nominal
    missions.stick_fixed_cruise  = base_mission_setup(analyses.stick_fixed_cruise,max_speed_multiplier,cruise_velocity,cruise_altitude) 
 
    return missions   

def elevator_sizing_setup(analyses,vehicle,cruise_velocity,cruise_altitude): 
    missions =RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier      = 1.4 # this multiplier is used to compute V_max from V_nominal
    missions.elevator_sizing  = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude)   
 
    return missions   

def aileron_rudder_sizing_setup(analyses,vehicle,cruise_velocity,cruise_altitude): 
    missions = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier      = 1.0     
    missions.aileron_sizing   = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude)  
    max_speed_multiplier      = 1.4   # this multiplier is used to compute V_max from V_nominal   
    missions.turn_criteria    = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude) 
 
    return missions   
    
def flap_sizing_setup(analyses,vehicle,cruise_velocity,cruise_altitude): 
    missions = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier     = 1.0      
    missions.flap_sizing     = base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude)   
    return missions        
    

# ------------------------------------------------------------------
#   Initialize the Mission
# ------------------------------------------------------------------    
    
def base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude):   
    '''
    This sets up the nominal cruise of the aircraft
    '''
     
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'
  
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    # base segment
    base_segment = Segments.Segment() 
    base_segment.state.numerics.number_control_points    = 3
 
    #   Cruise Segment: constant Speed, constant altitude 
    segment                           = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.analyses.extend( analyses )   
    segment.tag                       = "cruise"   
    segment.altitude                  = cruise_altitude
    segment.air_speed                 = cruise_velocity * max_speed_multiplier
    segment.distance                  =  20.   * Units.nautical_mile   
    mission.append_segment(segment)     
    
    return mission
