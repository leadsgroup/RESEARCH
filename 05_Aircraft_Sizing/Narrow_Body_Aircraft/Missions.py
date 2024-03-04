# Missions.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from   RCAIDE.Core                                        import Units
from   RCAIDE.Methods.Performance.estimate_stall_speed    import estimate_stall_speed
import numpy                                              as np

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------
    
def setup(analyses):
    """This allows multiple missions to be incorporated if desired, but only one is used here."""

    missions     = RCAIDE.Analyses.Mission.Missions()  

    base_mission = base_mission_setup(analyses)
    missions.base = base_mission 
    
    Initial_OEI_mission_1 = Initial_OEI_mission_setup_1(analyses)
    missions.Initial_OEI_1 = Initial_OEI_mission_1        
    
    Initial_OEI_mission_2 = Initial_OEI_mission_setup_2(analyses)
    missions.Initial_OEI_2 = Initial_OEI_mission_2
    
    Transitional_OEI_mission_1 = Transitional_OEI_mission_setup_1(analyses)
    missions.Transitional_OEI_1 = Transitional_OEI_mission_1        
    
    Transitional_OEI_mission_2 = Transitional_OEI_mission_setup_2(analyses)
    missions.Transitional_OEI_2 = Transitional_OEI_mission_2
    
    Established_OEI_mission_1 = Established_OEI_mission_setup_1(analyses)
    missions.Established_OEI_1 = Established_OEI_mission_1        
    
    Established_OEI_mission_2 = Established_OEI_mission_setup_2(analyses)
    missions.Established_OEI_2 = Established_OEI_mission_2
    
    Cruise_Climb_OEI_mission_1 = Cruise_Climb_OEI_mission_setup_1(analyses)
    missions.Cruise_Climb_OEI_1 = Cruise_Climb_OEI_mission_1        
    
    Cruise_Climb_OEI_mission_2 = Cruise_Climb_OEI_mission_setup_2(analyses)
    missions.Cruise_Climb_OEI_2 = Cruise_Climb_OEI_mission_2    
    
    Aborted_landing_AEO_mission_1 = Aborted_landing_AEO_mission_setup_1(analyses)
    missions.Aborted_landing_AEO_1 = Aborted_landing_AEO_mission_1
    
    Aborted_landing_AEO_mission_2 = Aborted_landing_AEO_mission_setup_2(analyses)
    missions.Aborted_landing_AEO_2 = Aborted_landing_AEO_mission_2   
    
    Aborted_landing_OEI_mission_1 = Aborted_landing_OEI_mission_setup_1(analyses)
    missions.Aborted_landing_OEI_1 = Aborted_landing_OEI_mission_1
    
    Aborted_landing_OEI_mission_2 = Aborted_landing_OEI_mission_setup_2(analyses)
    missions.Aborted_landing_OEI_2 = Aborted_landing_OEI_mission_2 
     
    missions.append(base_mission)
    missions.append(Initial_OEI_mission_1)
    missions.append(Initial_OEI_mission_2)
    missions.append(Transitional_OEI_mission_1)
    missions.append(Transitional_OEI_mission_2)
    missions.append(Established_OEI_mission_1)
    missions.append(Established_OEI_mission_2)
    missions.append(Cruise_Climb_OEI_mission_1)
    missions.append(Cruise_Climb_OEI_mission_2)
    missions.append(Aborted_landing_AEO_mission_1)
    missions.append(Aborted_landing_AEO_mission_2)
    missions.append(Aborted_landing_OEI_mission_1)
    missions.append(Aborted_landing_OEI_mission_2)

    return missions 
    
def base_mission_setup(analyses):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 

    # -------------------------------------------------------------------
    #   Takeoff Roll
    # -------------------------------------------------------------------

    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff" 
    segment.analyses.extend( analyses.takeoff )
    segment.velocity_start           = 100.* Units.knots
    segment.velocity_end             = 150 * Units.knots
    segment.friction_coefficient     = 0.04
    segment.altitude                 = 0.0   
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 3.048 * Units.km
    segment.air_speed                                    = 78.1956 * Units['m/s']
    segment.climb_rate                                   = 3000. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 3.657 * Units.km
    segment.air_speed                                    = 168.0 * Units['m/s']
    segment.climb_rate                                   = 2500. * Units['ft/min']     
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                  
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 25000. * Units.ft
    segment.air_speed                                    = 200.0  * Units['m/s']
    segment.climb_rate                                   = 1800. * Units['ft/min']   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Fourth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_4"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    segment.altitude_end                                  = 32000. * Units.ft
    segment.air_speed                                     = 230.0* Units['m/s']
    segment.climb_rate                                    = 900. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    # add to mission
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------
    #   Fifth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_5"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    segment.altitude_end                                 = 37000. * Units.ft
    segment.air_speed                                    = 230.0  * Units['m/s']
    segment.climb_rate                                   = 300.   * Units['ft/min']

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active              = True           
    segment.flight_controls.throttle.assigned_propulsors = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active            = True                
    
    # add to mission
    mission.append_segment(segment)    

    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.cruise )  
    segment.air_speed                                     = 450. * Units.knots
    segment.distance                                      = 2050. * Units.nmi
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1" 
    segment.analyses.extend( analyses.cruise )  
    segment.altitude_end                                  = 9.31  * Units.km
    segment.air_speed                                     = 440.0 * Units.knots
    segment.descent_rate                                  = 2600. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag  = "descent_2" 
    segment.analyses.extend( analyses.cruise_spoilers ) 

    segment.altitude_end                                  = 3.657 * Units.km
    segment.air_speed                                     = 365.0 * Units.knots
    segment.descent_rate                                  = 2300. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"  
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end                                  = 0.0   * Units.km
    segment.air_speed                                     = 250.0 * Units.knots
    segment.descent_rate                                  = 1500. * Units['ft/min'] 
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ----------------------------------------------------------------- 
    #   Landing Roll 
    # ----------------------------------------------------------------- 

    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"

    segment.analyses.extend( analyses.landing )
    segment.velocity_start                                      = 150 * Units.knots
    segment.velocity_end                                        = 100 * Units.knots 
    segment.friction_coefficient                                = 0.4
    segment.altitude                                            = 0.0   
    segment.flight_controls.elapsed_time.active                 = True  
    segment.flight_controls.elapsed_time.initial_guess_values   = [[30.]]  
    mission.append_segment(segment)     

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------
    
    #------------------------------------------------------------------
    ###         Reserve mission
    #------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Throttle
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_climb" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start = 0.0    * Units.km
    segment.altitude_end   = 15000. * Units.ft
    segment.air_speed      = 138.0  * Units['m/s']
    segment.climb_rate     = 3000.  * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)
      
    # ------------------------------------------------------------------
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude(base_segment)
    segment.tag = "reserve_cruise"
    
    segment.analyses.extend( analyses.cruise )
    
    segment.mach_number                                   = 0.5
    segment.distance                                      = 140.0 * Units.nautical_mile    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Loiter Segment: constant mach, constant time
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude_Loiter(base_segment)
    segment.tag = "reserve_loiter"

    segment.analyses.extend(analyses.cruise) 
    segment.mach_number                                   = 0.5
    segment.time                                          = 30.0 * Units.minutes
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #  Final Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------
    
    segment = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "reserve_descent_1"
    
    segment.analyses.extend( analyses.landing )
    
    segment.altitude_end                                  = 0.0  * Units.km
    segment.descent_rate                                  = 3.0  * Units['m/s']
    segment.mach_number_end                               = 0.24
    segment.mach_number_start                             = 0.3
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    #------------------------------------------------------------------
    ###         Reserve mission completed
    #------------------------------------------------------------------
    
    return mission 

def Initial_OEI_mission_setup_1(analyses):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'Initial_OEI_mission_setup_1'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 

    # -------------------------------------------------------------------
    #   Takeoff Roll
    # -------------------------------------------------------------------

    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff" 
    segment.analyses.extend( analyses.takeoff )
    segment.velocity_start           = 100.* Units.knots
    segment.velocity_end             = 150 * Units.knots
    segment.friction_coefficient     = 0.04
    segment.altitude                 = 0.0   
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Past A First Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1a" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 0.4572 * Units.km
    segment.air_speed                                    = 138.0 * Units['m/s']
    segment.climb_angle                                  = 0.6878 * Units.deg
    segment.true_course_angle                            = 0.0 * Units.degrees     
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Part B First Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1b" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start                               = 0.4572   * Units.km
    segment.altitude_end                                 = 3.048 * Units.km
    segment.air_speed                                    = 138.0 * Units['m/s']
    segment.climb_rate                                   = 3000. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)    

    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 3.657 * Units.km
    segment.air_speed                                    = 168.0 * Units['m/s']
    segment.climb_rate                                   = 2500. * Units['ft/min']     
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                  
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 25000. * Units.ft
    segment.air_speed                                    = 200.0  * Units['m/s']
    segment.climb_rate                                   = 1800. * Units['ft/min']   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Fourth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_4"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    segment.altitude_end                                  = 32000. * Units.ft
    segment.air_speed                                     = 230.0* Units['m/s']
    segment.climb_rate                                    = 900. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    # add to mission
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------
    #   Fifth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_5"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    segment.altitude_end                                 = 37000. * Units.ft
    segment.air_speed                                    = 230.0  * Units['m/s']
    segment.climb_rate                                   = 300.   * Units['ft/min']

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active              = True           
    segment.flight_controls.throttle.assigned_propulsors = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active            = True                
    
    # add to mission
    mission.append_segment(segment)    

    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.cruise )  
    segment.air_speed                                     = 450. * Units.knots
    segment.distance                                      = 2050. * Units.nmi
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1" 
    segment.analyses.extend( analyses.cruise )  
    segment.altitude_end                                  = 9.31  * Units.km
    segment.air_speed                                     = 440.0 * Units.knots
    segment.descent_rate                                  = 2600. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag  = "descent_2" 
    segment.analyses.extend( analyses.cruise_spoilers ) 

    segment.altitude_end                                  = 3.657 * Units.km
    segment.air_speed                                     = 365.0 * Units.knots
    segment.descent_rate                                  = 2300. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"  
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end                                  = 0.0   * Units.km
    segment.air_speed                                     = 250.0 * Units.knots
    segment.descent_rate                                  = 1500. * Units['ft/min'] 
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)

    # ----------------------------------------------------------------- 
    #   Landing Roll 
    # ----------------------------------------------------------------- 

    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"

    segment.analyses.extend( analyses.landing )
    segment.velocity_start                                      = 150 * Units.knots
    segment.velocity_end                                        = 100 * Units.knots 
    segment.friction_coefficient                                = 0.4
    segment.altitude                                            = 0.0   
    segment.flight_controls.elapsed_time.active                 = True  
    segment.flight_controls.elapsed_time.initial_guess_values   = [[30.]]  
    mission.append_segment(segment)     

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------
    
    #------------------------------------------------------------------
    ###         Reserve mission
    #------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Throttle
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_climb" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start = 0.0    * Units.km
    segment.altitude_end   = 15000. * Units.ft
    segment.air_speed      = 138.0  * Units['m/s']
    segment.climb_rate     = 3000.  * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)
      
    # ------------------------------------------------------------------
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude(base_segment)
    segment.tag = "reserve_cruise"
    
    segment.analyses.extend( analyses.cruise )
    
    segment.mach_number                                   = 0.5
    segment.distance                                      = 140.0 * Units.nautical_mile    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Loiter Segment: constant mach, constant time
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude_Loiter(base_segment)
    segment.tag = "reserve_loiter"

    segment.analyses.extend(analyses.cruise) 
    segment.mach_number                                   = 0.5
    segment.time                                          = 30.0 * Units.minutes
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #  Final Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------
    
    segment = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "reserve_descent_1"
    
    segment.analyses.extend( analyses.landing )
    
    segment.altitude_end                                  = 0.0  * Units.km
    segment.descent_rate                                  = 3.0  * Units['m/s']
    segment.mach_number_end                               = 0.24
    segment.mach_number_start                             = 0.3
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    #------------------------------------------------------------------
    ###         Reserve mission completed
    #------------------------------------------------------------------
    
    return mission 

def AEO_mission_setup(analyses, vehicle):
    
   
    
    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)    
    

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 
    
    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 10.)
    base_segment.temperature_deviation = atmo_data.temperature[0][0]     

 
    # -------------------------------------------------------------------
    #   Takeoff Roll
    # -------------------------------------------------------------------

    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff" 
    segment.analyses.extend( analyses.takeoff )
    segment.velocity_start           = 100.* Units.knots
    segment.velocity_end             = 150 * Units.knots
    segment.friction_coefficient     = 0.04
    segment.altitude                 = 0.0   
    mission.append_segment(segment)
    
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 3.048 * Units.km
    segment.air_speed                                    = 138.0 * Units['m/s']
    segment.climb_rate                                   = 3000. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 3.657 * Units.km
    segment.air_speed                                    = 168.0 * Units['m/s']
    segment.climb_rate                                   = 2500. * Units['ft/min']     
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                  
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 25000. * Units.ft
    segment.air_speed                                    = 200.0  * Units['m/s']
    segment.climb_rate                                   = 1800. * Units['ft/min']   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    # ------------------------------------------------------------------
    #   Fourth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_4"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    
    segment.altitude_end                                  = 32000. * Units.ft
    segment.air_speed                                     = 230.0* Units['m/s']
    segment.climb_rate                                    = 900. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    
    # add to mission
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------
    #   Fifth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_5"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    segment.altitude_end                                 = 37000. * Units.ft
    segment.air_speed                                    = 230.0  * Units['m/s']
    segment.climb_rate                                   = 300.   * Units['ft/min']

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active              = True           
    segment.flight_controls.throttle.assigned_propulsors = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active            = True                
    
    
    # add to mission
    mission.append_segment(segment)    



    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.cruise )  
    segment.air_speed                                     = 450. * Units.knots
    segment.distance                                      = 2050. * Units.nmi
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1" 
    segment.analyses.extend( analyses.cruise )  
    segment.altitude_end                                  = 9.31  * Units.km
    segment.air_speed                                     = 440.0 * Units.knots
    segment.descent_rate                                  = 2600. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag  = "descent_2" 
    segment.analyses.extend( analyses.cruise_spoilers ) 

    segment.altitude_end                                  = 3.657 * Units.km
    segment.air_speed                                     = 365.0 * Units.knots
    segment.descent_rate                                  = 2300. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"  
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end                                  = 0.0   * Units.km
    segment.air_speed                                     = 250.0 * Units.knots
    segment.descent_rate                                  = 1500. * Units['ft/min'] 
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
 

    # ----------------------------------------------------------------- 
    #   Landing Roll 
    # ----------------------------------------------------------------- 

    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"

    segment.analyses.extend( analyses.landing )
    segment.velocity_start                                      = 150 * Units.knots
    segment.velocity_end                                        = 100 * Units.knots 
    segment.friction_coefficient                                = 0.4
    segment.altitude                                            = 0.0   
    segment.flight_controls.elapsed_time.active                 = True  
    segment.flight_controls.elapsed_time.initial_guess_values   = [[30.]]  
    mission.append_segment(segment)     

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------
    
    
    #------------------------------------------------------------------
    ###         Reserve mission
    #------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Throttle + FAR Part 25.119 requirement
    # ------------------------------------------------------------------


    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_climb" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start = 0.0    * Units.km
    segment.altitude_end   = 15000. * Units.ft
    segment.air_speed      = 138.0  * Units['m/s']    #(1.23*V_SR0)
    segment.climb_rate     = 868.85  * Units['ft/min'] #(Vv = V*sin(gamma))
    # gradient of climb minimum 3.2 percent = 0.03199 rad
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)
      
    # ------------------------------------------------------------------
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude(base_segment)
    segment.tag = "reserve_cruise"
    
    segment.analyses.extend( analyses.cruise )
    
    segment.mach_number                                   = 0.5
    segment.distance                                      = 140.0 * Units.nautical_mile    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Loiter Segment: constant mach, constant time
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude_Loiter(base_segment)
    segment.tag = "reserve_loiter"

    segment.analyses.extend(analyses.cruise) 
    segment.mach_number                                   = 0.5
    segment.time                                          = 30.0 * Units.minutes
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #  Final Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------
    
    segment = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "reserve_descent_1"
    
    segment.analyses.extend( analyses.landing )
    
    segment.altitude_end                                  = 0.0  * Units.km
    segment.descent_rate                                  = 3.0  * Units['m/s']
    segment.mach_number_end                               = 0.24
    segment.mach_number_start                             = 0.3
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    #------------------------------------------------------------------
    ###         Reserve mission completed
    #------------------------------------------------------------------
    
    return mission 

def OEI_mission_setup(analyses):
    

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 

 
    # -------------------------------------------------------------------
    #   Takeoff Roll
    # -------------------------------------------------------------------

    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff" 
    segment.analyses.extend( analyses.takeoff )
    segment.velocity_start           = 100.* Units.knots
    segment.velocity_end             = 150 * Units.knots
    segment.friction_coefficient     = 0.04
    segment.altitude                 = 0.0   
    mission.append_segment(segment)
    
    
    # ------------------------------------------------------------------
    #   First Climb Segment A: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1A" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 0.4572 * Units.km
    segment.air_speed                                    = 138.0 * Units['m/s']
    segment.climb_rate                                   = 325.96 * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Climb Segment B: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1B" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start                               = 0.4572* Units.km
    segment.altitude_end                                 = 3.048 * Units.km
    segment.air_speed                                    = 138.0 * Units['m/s']
    segment.climb_rate                                   = 3000. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)    


    # ------------------------------------------------------------------
    #   Second Climb Segment B: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 3.657 * Units.km
    segment.air_speed                                    = 168.0 * Units['m/s']
    segment.climb_rate                                   = 2500. * Units['ft/min']     
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                  
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3" 
    segment.analyses.extend( analyses.cruise ) 

    segment.altitude_end                                 = 25000. * Units.ft
    segment.air_speed                                    = 200.0  * Units['m/s']
    segment.climb_rate                                   = 1800. * Units['ft/min']   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    # ------------------------------------------------------------------
    #   Fourth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_4"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    
    segment.altitude_end                                  = 32000. * Units.ft
    segment.air_speed                                     = 230.0* Units['m/s']
    segment.climb_rate                                    = 900. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    
    # add to mission
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------
    #   Fifth Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Constant_Speed_Constant_Rate()
    segment.tag = "climb_5"
    
    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )
     
    segment.altitude_end                                 = 37000. * Units.ft
    segment.air_speed                                    = 230.0  * Units['m/s']
    segment.climb_rate                                   = 300.   * Units['ft/min']

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active              = True           
    segment.flight_controls.throttle.assigned_propulsors = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active            = True                
    
    
    # add to mission
    mission.append_segment(segment)    



    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.cruise )  
    segment.air_speed                                     = 450. * Units.knots
    segment.distance                                      = 2050. * Units.nmi
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1" 
    segment.analyses.extend( analyses.cruise )  
    segment.altitude_end                                  = 9.31  * Units.km
    segment.air_speed                                     = 440.0 * Units.knots
    segment.descent_rate                                  = 2600. * Units['ft/min']
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag  = "descent_2" 
    segment.analyses.extend( analyses.cruise_spoilers ) 

    segment.altitude_end                                  = 3.657 * Units.km
    segment.air_speed                                     = 365.0 * Units.knots
    segment.descent_rate                                  = 2300. * Units['ft/min']
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed Constant Rate  
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"  
    segment.analyses.extend( analyses.landing ) 
    segment.altitude_end                                  = 0.0   * Units.km
    segment.air_speed                                     = 250.0 * Units.knots
    segment.descent_rate                                  = 1500. * Units['ft/min'] 
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
 

    # ----------------------------------------------------------------- 
    #   Landing Roll 
    # ----------------------------------------------------------------- 

    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"

    segment.analyses.extend( analyses.landing )
    segment.velocity_start                                      = 150 * Units.knots
    segment.velocity_end                                        = 100 * Units.knots 
    segment.friction_coefficient                                = 0.4
    segment.altitude                                            = 0.0   
    segment.flight_controls.elapsed_time.active                 = True  
    segment.flight_controls.elapsed_time.initial_guess_values   = [[30.]]  
    mission.append_segment(segment)     

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------
    
    
    #------------------------------------------------------------------
    ###         Reserve mission
    #------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Throttle + FAR Part 25.119 requirement
    # ------------------------------------------------------------------


    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_climb" 
    segment.analyses.extend( analyses.takeoff ) 
    segment.altitude_start = 0.0    * Units.km
    segment.altitude_end   = 15000. * Units.ft
    segment.air_speed      = 138.0  * Units['m/s']    #(1.23*V_SR0)
    segment.climb_rate     = 868.85  * Units['ft/min'] #(Vv = V*sin(gamma))
    # gradient of climb minimum 3.2 percent = 0.03199 rad
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)
      
    # ------------------------------------------------------------------
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude(base_segment)
    segment.tag = "reserve_cruise"
    
    segment.analyses.extend( analyses.cruise )
    
    segment.mach_number                                   = 0.5
    segment.distance                                      = 140.0 * Units.nautical_mile    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Loiter Segment: constant mach, constant time
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Mach_Constant_Altitude_Loiter(base_segment)
    segment.tag = "reserve_loiter"

    segment.analyses.extend(analyses.cruise) 
    segment.mach_number                                   = 0.5
    segment.time                                          = 30.0 * Units.minutes
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #  Final Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------
    
    segment = Segments.Descent.Linear_Mach_Constant_Rate(base_segment)
    segment.tag = "reserve_descent_1"
    
    segment.analyses.extend( analyses.landing )
    
    segment.altitude_end                                  = 0.0  * Units.km
    segment.descent_rate                                  = 3.0  * Units['m/s']
    segment.mach_number_end                               = 0.24
    segment.mach_number_start                             = 0.3
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                
    
    mission.append_segment(segment)
    
    #------------------------------------------------------------------
    ###         Reserve mission completed
    #------------------------------------------------------------------
    
    return mission 