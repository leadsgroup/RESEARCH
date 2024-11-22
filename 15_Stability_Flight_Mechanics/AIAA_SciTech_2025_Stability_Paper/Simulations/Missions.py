# Missions.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from RCAIDE.Framework.Core import  Units 
from RCAIDE.Library.Methods.Performance.estimate_stall_speed       import estimate_stall_speed 


# ############################################################################################################################################################     
#  MISSION CASES 
# ############################################################################################################################################################     
def stick_fixed_stability_setup(analyses,cruise_velocity,cruise_altitude,angle_of_attack = 0 ,sideslip_angle=0,bank_angle= 0, roll_rate = 0,pitch_rate = 0,yaw_rate = 0): 
    missions                     = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier         = 1.0 
    missions.stick_fixed_cruise  = base_mission_setup(analyses.stick_fixed_cruise,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate) 
 
    return missions   

def elevator_sizing_setup(analyses,cruise_velocity,cruise_altitude,angle_of_attack = 0 ,sideslip_angle=0,bank_angle= 0, roll_rate = 0,pitch_rate = 0,yaw_rate = 0): 
    missions                            = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier                = 1.87 
    vehicle_mass                        = analyses.elevator_sizing_pull_up.aerodynamics.vehicle.mass_properties.max_takeoff
    reference_area                      = analyses.elevator_sizing_pull_up.aerodynamics.vehicle.reference_area 
    stall_velocity                      = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)   
    missions.elevator_sizing_pull_up    = base_mission_setup(analyses.elevator_sizing_pull_up,max_speed_multiplier,stall_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate)   

    max_speed_multiplier                = 1.87
    missions.elevator_sizing_push_over  = base_mission_setup(analyses.elevator_sizing_push_over,max_speed_multiplier,stall_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate)   
    return missions   

def aileron_rudder_sizing_setup(analyses,cruise_velocity,cruise_altitude,angle_of_attack = 0,sideslip_angle=0 ,bank_angle= 0, roll_rate = 0,pitch_rate = 0,yaw_rate = 0): 
    missions                                    = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier                        = 1.0 
    missions.roll_maneuver                      = base_mission_setup(analyses.aileron_rudder_roll_sizing,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate)    
    
    max_speed_multiplier                        = 1.0
    missions.crosswind_maneuver                 = base_mission_setup(analyses.aileron_rudder_crosswind_sizing,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate)
    return missions   
     

def aileron_rudder_oei_sizing_setup(analyses,cruise_velocity,cruise_altitude,angle_of_attack = 0,sideslip_angle=0 ,bank_angle= 0, roll_rate = 0,pitch_rate = 0,yaw_rate = 0): 
    missions                                    = RCAIDE.Framework.Mission.Missions() 
    max_speed_multiplier                        = 1.0
    missions.cruise_oei                         = cruise_oei_mission_setup(analyses.aileron_rudder_oei_sizing,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate)     
    return missions   


def flap_sizing_setup(analyses,cruise_velocity,cruise_altitude,angle_of_attack = 0,sideslip_angle=0 ,bank_angle= 0, roll_rate = 0,pitch_rate = 0,yaw_rate = 0): 
    missions                            = RCAIDE.Framework.Mission.Missions()
    max_speed_multiplier                = 1.4
    angle_of_attack                     = 12.*Units.degrees  
    vehicle_mass                        = analyses.flap_sizing_flaps_up.aerodynamics.vehicle.mass_properties.max_takeoff
    reference_area                      = analyses.flap_sizing_flaps_up.aerodynamics.vehicle.reference_area 
    stall_velocity                      = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)
    
    missions.flap_sizing_flaps_up       = base_mission_setup(analyses.flap_sizing_flaps_up,max_speed_multiplier,stall_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate) 
    
    missions.flap_sizing_flaps_down     = base_mission_setup(analyses.flap_sizing_flaps_down,max_speed_multiplier,stall_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate)  
    
    return missions        
    
    
def hover_oei_setup(analyses,cruise_velocity,cruise_altitude,angle_of_attack = 0,sideslip_angle=0 ,bank_angle= 0, roll_rate = 0,pitch_rate = 0,yaw_rate = 0): 
    missions                            = RCAIDE.Framework.Mission.Missions() 
    missions.hover_prop_rotor_oei       = hover_oei_mission_setup(analyses.hover_prop_rotor_oei)  
    missions.hover_lift_rotor_oei       = hover_oei_mission_setup(analyses.hover_lift_rotor_oei)   
    return missions

# ############################################################################################################################################################     
#  MISSION SETUP  
# ############################################################################################################################################################
 
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

def cruise_oei_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate):   
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

    segment.assigned_control_variables.throttle.active               = True
    segment.assigned_control_variables.throttle.initial_guess_values = [[0.5]]
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['prop_rotor_propulsor_1','prop_rotor_propulsor_2','prop_rotor_propulsor_3','prop_rotor_propulsor_4','prop_rotor_propulsor_5']]
    
        
    mission.append_segment(segment)     
    
    return mission



def hover_oei_mission_setup(analyses):   
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
    segment.tag                       = "hover_oei"  
    segment.altitude                  = 20 * Units.feet 
 
    segment.flight_dynamics.force_x   = True    
    segment.flight_dynamics.force_z   = True    
    segment.flight_dynamics.force_y   = True         
    segment.flight_dynamics.moment_y  = True 
    segment.flight_dynamics.moment_x  = True
    segment.flight_dynamics.moment_z  = True

    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['prop_rotor_propulsor_1','prop_rotor_propulsor_2','prop_rotor_propulsor_3',
                                                                         'prop_rotor_propulsor_4','prop_rotor_propulsor_5','prop_rotor_propulsor_6'],
                                                                        ['lift_rotor_propulsor_1','lift_rotor_propulsor_2','lift_rotor_propulsor_3',
                                                                         'lift_rotor_propulsor_4', 'lift_rotor_propulsor_5', 'lift_rotor_propulsor_6']]
        
    mission.append_segment(segment)     
      
    return mission

 