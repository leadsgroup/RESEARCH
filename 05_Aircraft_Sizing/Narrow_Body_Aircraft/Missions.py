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
    missions.append(base_mission)
    
    # FAA Rule FAR 25.111c_3 
    OEI_mission_111c_1 = OEI_mission_setup_111c_1(analyses) # standard atm, 15 degree C 
    missions.append(OEI_mission_111c_1)    
    OEI_mission_111c_2 = OEI_mission_setup_111c_2(analyses) # 25 degree C 
    missions.append(OEI_mission_111c_2)
     
    # FAA Rule FAR 25.121a_b_c     
    OEI_mission_121_1 = OEI_mission_setup_121_1(analyses) # standard atm, 15 degree C   
    missions.append(OEI_mission_121_1)  
    OEI_mission_121_2 = OEI_mission_setup_121_2(analyses)# 25 degree C
    missions.append(OEI_mission_121_2)   

    # FAA Rule FAR 25.119    
    Aborted_landing_AEO_mission_119_1 = Aborted_landing_AEO_mission_setup_119_1(analyses) # standard atm, 15 degree C 
    missions.append(Aborted_landing_AEO_mission_119_1) 
    Aborted_landing_AEO_mission_119_2 = Aborted_landing_AEO_mission_setup_119_2(analyses) # 25 degree C 
    missions.append(Aborted_landing_AEO_mission_119_2)
    
    # FAA Rule FAR 25.121d     
    Aborted_landing_OEI_mission_121d_1 = Aborted_landing_OEI_mission_setup_121d_1(analyses) # standard atm, 15 degree C
    missions.append(Aborted_landing_OEI_mission_121d_1)
    Aborted_landing_OEI_mission_121d_2 = Aborted_landing_OEI_mission_setup_121d_2(analyses) # 25 degree C 
    missions.append(Aborted_landing_OEI_mission_121d_2)
    
    # Takeoff distance requirement
    TO_mission_setup_flap15_h0_mission = TO_mission_setup_flap15_h0(analyses)
    missions.append(TO_mission_setup_flap15_h0_mission) 

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

def OEI_mission_setup_111c_1(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'OEI_mission_setup_111c_1'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 
    
    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 0.0)
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
    #   Past A First Climb Segment: Constant Speed Constant Angle 
    # ------------------------------------------------------------------
    
    # VSTALL Calculation
    vehicle_mass   = analyses.takeoff_gear_up.aerodynamics.geometry.mass_properties.max_takeoff
    reference_area = analyses.takeoff_gear_up.aerodynamics.geometry.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)          

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1a" 
    segment.analyses.extend( analyses.takeoff_gear_up ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 1500 * Units.feet
    segment.air_speed                                    = 1.2*Vstall * Units['m/s']
    segment.climb_angle                                  = 0.6878 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)    
  
    
    return mission 

def OEI_mission_setup_111c_2(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'OEI_mission_setup_111c_2'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 
    
    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 10.0)
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
    #   Past A First Climb Segment: Constant Speed Constant Angle 
    # ------------------------------------------------------------------
    
    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)          

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1a" 
    segment.analyses.extend( analyses.takeoff_gear_up ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 1500 * Units.feet
    segment.air_speed                                    = 1.2*Vstall * Units['m/s']
    segment.climb_angle                                  = 0.6878 * Units.deg 
     
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
    segment.altitude_end                                 = 3.048 * Units.km
    segment.air_speed                                    = 78.1956 * Units['m/s']
    segment.climb_rate                                   = 3000. * Units['ft/min']    
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)    
    return mission

def OEI_mission_setup_121_1(analyses,vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'OEI_mission_setup_121_1'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 
    
    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 0.0)
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
    #   Past A First Climb Segment: Constant Speed Constant Angle 
    # ------------------------------------------------------------------
    
    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)          

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1a" 
    segment.analyses.extend( analyses.takeoff_gear_down ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 100 
    segment.air_speed                                    = 1.1*Vstall * Units['m/s']
    segment.climb_angle                                  = 0.1 * Units.deg    
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Part B First Climb Segment: Constant Speed Constant Angle 
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1b" 
    segment.analyses.extend( analyses.takeoff_gear_up ) 
    segment.altitude_end                                 = 0.4572 * Units.km
    segment.air_speed                                    = 1.2*Vstall * Units['m/s']
    segment.climb_angle                                  = 1.3755321172 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment) 
    
    # ------------------------------------------------------------------
    #   Part C First Climb Segment: Constant Speed Constant Angle  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1c" 
    segment.analyses.extend( analyses.takeoff_no_flaps )  
    segment.altitude_end                                 = 3.048 * Units.km
    segment.air_speed                                    = 1.2*Vstall * Units['m/s']
    segment.climb_angle                                  = 0.68786507006 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)       

    
    return mission

def OEI_mission_setup_121_2(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'OEI_mission_setup_121_2'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 
    
    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 10.0)
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
    #   Past A First Climb Segment: Constant Speed Constant Angle 
    # ------------------------------------------------------------------
    
    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)          

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1a" 
    segment.analyses.extend( analyses.takeoff_gear_down ) 
    segment.altitude_start                               = 0.0   * Units.km
    segment.altitude_end                                 = 100 
    segment.air_speed                                    = 1.1*Vstall * Units['m/s']
    segment.climb_angle                                  = 0.1 * Units.deg    
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Part B First Climb Segment: Constant Speed Constant Angle 
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1b" 
    segment.analyses.extend( analyses.takeoff_gear_up ) 
    segment.altitude_end                                 = 0.4572 * Units.km
    segment.air_speed                                    = 1.2*Vstall * Units['m/s']
    segment.climb_angle                                  = 1.3755321172 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment) 
    
    # ------------------------------------------------------------------
    #   Part C First Climb Segment: Constant Speed Constant Angle  
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "climb_1c" 
    segment.analyses.extend( analyses.takeoff_no_flaps )  
    segment.altitude_end                                 = 3.048 * Units.km
    segment.air_speed                                    = 1.2*Vstall * Units['m/s']
    segment.climb_angle                                  = 0.68786507006 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)       
      
    
    return mission

def Aborted_landing_AEO_mission_setup_119_1(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'Aborted_landing_AEO_mission_setup_119_1'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 
    
    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 0.0)
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
    
    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff*0.85
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)        
 
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Angle
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "reserve_climb" 
    segment.analyses.extend( analyses.landing_AEO ) 
    segment.altitude_start                               = 0.0    * Units.km
    segment.altitude_end                                 = 15000. * Units.ft
    segment.air_speed                                    = 1.3*Vstall  * Units['m/s']
    segment.climb_angle                                  = 1.8337691465 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment) 
   
    return mission 

def Aborted_landing_AEO_mission_setup_119_2(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'Aborted_landing_AEO_mission_setup_119_2'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 
    
    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 10.0)
    base_segment.temperature_deviation = atmo_data.temperature[0][0]      

   
    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff*0.85
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)        
 
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Angle
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "reserve_climb" 
    segment.analyses.extend( analyses.landing_AEO ) 
    segment.altitude_start                               = 0.0    * Units.km
    segment.altitude_end                                 = 15000. * Units.ft
    segment.air_speed                                    = 1.3*Vstall  * Units['m/s']
    segment.climb_angle                                  = 1.8337691465 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment) 
    return mission 

def Aborted_landing_OEI_mission_setup_121d_1(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'Aborted_landing_OEI_mission_setup_121d_1'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 

    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 0.0)
    base_segment.temperature_deviation = atmo_data.temperature[0][0]  

    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff*0.85
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)        
 
    # ------------------------------------------------------------------
    #  A First Climb Segment: Constant Speed, Constant Angle
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "reserve_climb_a" 
    segment.analyses.extend( analyses.landing_OEI ) 
    segment.altitude_start                               = 0.0    * Units.km
    segment.altitude_end                                 = 1500. * Units.ft
    segment.air_speed                                    = 1.4*Vstall  * Units['m/s']
    segment.climb_angle                                  = 1.20364474013 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)  
   
    return mission 

def Aborted_landing_OEI_mission_setup_121d_2(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'Aborted_landing_OEI_mission_setup_121d_2'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 

    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data          = atmosphere.compute_values(altitude = 0,temperature_deviation= 10.0)
    base_segment.temperature_deviation = atmo_data.temperature[0][0]      

    # VSTALL Calculation
    vehicle_mass   = vehicle.mass_properties.max_takeoff*0.85
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)    
  
    # ------------------------------------------------------------------
    #  A First Climb Segment: Constant Speed, Constant Angle
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Angle(base_segment)
    segment.tag = "reserve_climb_a" 
    segment.analyses.extend( analyses.landing_OEI ) 
    segment.altitude_start                               = 0.0    * Units.km
    segment.altitude_end                                 = 1500. * Units.ft
    segment.air_speed                                    = 1.4*Vstall  * Units['m/s']
    segment.climb_angle                                  = 1.20364474013 * Units.deg 
     
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor']] 
    segment.flight_controls.body_angle.active             = True                 
    
    mission.append_segment(segment)  

    return mission 

def TO_mission_setup_flap15_h0(analyses, vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'TO_mission_flap15_h0'
  
    Segments = RCAIDE.Analyses.Mission.Segments 
    base_segment = Segments.Segment() 

    # -------------------------------------------------------------------
    #   Takeoff Roll
    # -------------------------------------------------------------------

    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "TO_mission_flap15_h0" 
    segment.analyses.extend( analyses.TO_mission_flap15_h0 )
    segment.velocity_start           = 0.* Units.knots
    S_TO_lim                         = 5699.38 * Units.feet
    vehicle.mass_properties.max_takeoff = 114224.9491 * Units.pounds
    W                                = vehicle.mass_properties.max_takeoff
    T                                = 2*vehicle.turbofan.design_thrust
    segment.friction_coefficient     = 0.04
    C_D_TO                           = 0.02
    C_L_TO                           = 2.2
    rho                              = 0.002
    g                                = 32.17405
    segment.velocity_end             = np.sqrt(((2*g*S_TO_lim*T/W) - 
                                                (segment.friction_coefficient*2*g*S_TO_lim))/
                                               (1 + (S_TO_lim*g*rho*vehicle.wing.areas.reference*C_D_TO/(2*W) - 
                                                     g*S_TO_lim*segment.friction_coefficient*rho*
                                                     vehicle.wing.areas.reference*C_L_TO/(2*W)))) * Units.knots
    segment.altitude                 = 0.0   
    mission.append_segment(segment)
    
    return mission