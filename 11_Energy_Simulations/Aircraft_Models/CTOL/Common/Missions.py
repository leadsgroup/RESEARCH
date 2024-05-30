'''

# Missions.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke  

'''

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units    
import numpy as np 
from RCAIDE.Library.Methods.Performance.estimate_stall_speed   import estimate_stall_speed 
from RCAIDE.Library.Methods.Utilities.Chebyshev  import chebyshev_data
from RCAIDE.Library.Methods.Utilities.Chebyshev  import linear_data

# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def baseline_mission_setup(analyses,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'baseline_mission'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport      

    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976() 
    atmo_data          = atmosphere.compute_values(altitude = airport.altitude,temperature_deviation= 1.)     

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                             = Segments.Segment() 
    ones_row                                                 = base_segment.state.ones_row
    base_segment.battery_discharge                           = True  
    base_segment.state.numerics.number_control_points        = control_points 
    base_segment.state.numerics.discretization_method        = linear_data
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability     = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle        = analyses.base.aerodynamics.geometry
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)         

    # ------------------------------------------------------------------
    #   Takeoff
    # ------------------------------------------------------------------      
    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff"  
    segment.analyses.extend( analyses.base )
    segment.velocity_start                                   = Vstall*0.1  
    segment.velocity_end                                     = Vstall  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 10.            
    segment.altitude                                         = 0.0  
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle   
    segment.battery_energy                                   = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)         
    mission.append_segment(segment) 
     
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Departure_End_of_Runway'       
    segment.analyses.extend( analyses.base )           
    segment.altitude_start                                   = 0.0 * Units.feet
    segment.altitude_end                                     = 50.0 * Units.feet
    segment.air_speed_start                                  = Vstall  
    segment.air_speed_end                                    = Vstall*1.1        
    segment.climb_rate                                       = 600 * Units['ft/min'] 
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3])    
    mission.append_segment(segment)  
                
    #------------------------------------------------------------------
    #  Initial Climb Area Segment Flight 1  
    #------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Initial_CLimb_Area'  
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                   = 50.0 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1     
    segment.air_speed_end                                    = Vstall*1.2  
    segment.climb_rate                                       = 600 * Units['ft/min']   
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle        
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])  
    mission.append_segment(segment)  
     
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Climb'       
    segment.analyses.extend( analyses.base )         
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 1500 * Units.feet 
    segment.air_speed_start                                  = Vstall*1.2  
    segment.air_speed_end                                    = 150.* Units['mph']    
    segment.climb_rate                                       = 500* Units['ft/min']
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle       
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])   
    mission.append_segment(segment)  
    
    # ------------------------------------------------------------------
    #   Cruise Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag = 'Cruise'  
    segment.analyses.extend(analyses.base) 
    segment.air_speed                                        = 150.* Units['mph']        
    cruise_distance                                          = airport_geospacial_data.flight_range - 15.5415*Units.nmi  
    segment.altitude                                         = 1500 * Units.feet 
    segment.distance                                         = cruise_distance          
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.8]) 
    mission.append_segment(segment)   
 
    
    # ------------------------------------------------------------------
    #   Descent Segment Flight 1   
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Descent'   
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                   = 1500 * Units.feet  
    segment.altitude_end                                     = 1000.0 * Units.feet
    segment.air_speed_start                                  = 150.* Units['mph']    
    segment.air_speed_end                                    = Vstall*1.3     
    segment.climb_rate                                       = -300 * Units['ft/min']    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.8])  
    mission.append_segment(segment)  
     

    # ------------------------------------------------------------------
    #  Downleg_Altitude Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment) 
    segment.tag = 'Downleg' 
    segment.analyses.extend(analyses.base) 
    segment.air_speed_start                                  = Vstall*1.3 
    segment.air_speed_end                                    = Vstall*1.2       
    segment.distance                                         =  6000 * Units.feet
    segment.acceleration                                     = -0.05 * Units['m/s/s']     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.7])  
    mission.append_segment(segment)        
    
    # ------------------------------------------------------------------
    #  Baseleg Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = 'Baseleg' 
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                                   = 1000 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.2
    segment.air_speed_end                                    = Vstall*1.1  
    segment.climb_rate                                       = -300 * Units['ft/min']  
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment ,  initial_throttles = [0.7])
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #  Final Approach Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment_name = 'Final_Approach' 
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)            
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 00.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1  
    segment.air_speed_end                                    = Vstall 
    segment.climb_rate                                       = -300 * Units['ft/min'] 
    segment.state.unknowns.throttle                          =  0.8 * ones_row(1)    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3],  initial_throttles = [0.8] )  
    mission.append_segment(segment)   
    

    # ------------------------------------------------------------------
    #   Landing  
    # ------------------------------------------------------------------  
    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"   
    segment.analyses.extend( analyses.base) 
    segment.velocity_start                                   = Vstall  
    segment.velocity_end                                     = Vstall*0.1  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 30.            
    segment.altitude                                         = 0.0  
    segment.state.unknowns.velocity_x                        = 0.1* Vstall * ones_row(1)    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.3])    
    mission.append_segment(segment)   

            
    return mission 


# ------------------------------------------------------------------
#   Full Mission 1500
# ------------------------------------------------------------------
def direct_mission_setup_at_1500ft(analyses,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'direct_mission_setup_at_1500ft'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport      

    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976() 
    atmo_data          = atmosphere.compute_values(altitude = airport.altitude,temperature_deviation= 1.)     

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                             = Segments.Segment() 
    ones_row                                                 = base_segment.state.ones_row
    base_segment.battery_discharge                           = True  
    base_segment.state.numerics.number_control_points        = control_points 
    base_segment.state.numerics.discretization_method        = linear_data
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability     = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle        = analyses.base.aerodynamics.geometry
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)         

    # ------------------------------------------------------------------
    #   Takeoff
    # ------------------------------------------------------------------      
    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff"  
    segment.analyses.extend( analyses.base )
    segment.velocity_start                                   = Vstall*0.1  
    segment.velocity_end                                     = Vstall  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 10.            
    segment.altitude                                         = 0.0  
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle   
    segment.battery_energy                                   = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)         
    mission.append_segment(segment) 
     
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Departure_End_of_Runway'       
    segment.analyses.extend( analyses.base )           
    segment.altitude_start                                   = 0.0 * Units.feet
    segment.altitude_end                                     = 50.0 * Units.feet
    segment.air_speed_start                                  = Vstall  
    segment.air_speed_end                                    = Vstall*1.1        
    segment.climb_rate                                       = 600 * Units['ft/min'] 
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3])    
    mission.append_segment(segment)  
                
    #------------------------------------------------------------------
    #  Initial Climb Area Segment Flight 1  
    #------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Initial_CLimb_Area'  
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                   = 50.0 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1     
    segment.air_speed_end                                    = Vstall*1.2  
    segment.climb_rate                                       = 600 * Units['ft/min']   
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle        
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])  
    mission.append_segment(segment)  
     
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Climb'       
    segment.analyses.extend( analyses.base )         
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 1500 * Units.feet 
    segment.air_speed_start                                  = Vstall*1.2  
    segment.air_speed_end                                    = 150.* Units['mph']    
    segment.climb_rate                                       = 500* Units['ft/min']
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle       
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])   
    mission.append_segment(segment)  
    
    # ------------------------------------------------------------------
    #   Cruise Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag = 'Cruise'  
    segment.analyses.extend(analyses.base) 
    segment.air_speed                                        = 150.* Units['mph']        
    cruise_distance                                          = airport_geospacial_data.flight_range - 15.5415*Units.nmi  
    segment.altitude                                         = 1500 * Units.feet 
    segment.distance                                         = cruise_distance          
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.8]) 
    mission.append_segment(segment)   
 
    
    # ------------------------------------------------------------------
    #   Descent Segment Flight 1   
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Descent'   
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                   = 1500 * Units.feet  
    segment.altitude_end                                     = 1000.0 * Units.feet
    segment.air_speed_start                                  = 150.* Units['mph']    
    segment.air_speed_end                                    = Vstall*1.3     
    segment.climb_rate                                       = -300 * Units['ft/min']    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.8])  
    mission.append_segment(segment)  
     

    # ------------------------------------------------------------------
    #  Downleg_Altitude Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment) 
    segment.tag = 'Downleg' 
    segment.analyses.extend(analyses.base) 
    segment.air_speed_start                                  = Vstall*1.3 
    segment.air_speed_end                                    = Vstall*1.2       
    segment.distance                                         =  6000 * Units.feet
    segment.acceleration                                     = -0.05 * Units['m/s/s']     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.7])  
    mission.append_segment(segment)        
    
    # ------------------------------------------------------------------
    #  Baseleg Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = 'Baseleg' 
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                                   = 1000 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.2
    segment.air_speed_end                                    = Vstall*1.1  
    segment.climb_rate                                       = -300 * Units['ft/min']  
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment ,  initial_throttles = [0.7])
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #  Final Approach Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment_name = 'Final_Approach' 
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)            
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 00.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1  
    segment.air_speed_end                                    = Vstall 
    segment.climb_rate                                       = -300 * Units['ft/min'] 
    segment.state.unknowns.throttle                          =  0.8 * ones_row(1)    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3],  initial_throttles = [0.8] )  
    mission.append_segment(segment)   
    

    # ------------------------------------------------------------------
    #   Landing  
    # ------------------------------------------------------------------  
    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"   
    segment.analyses.extend( analyses.base) 
    segment.velocity_start                                   = Vstall  
    segment.velocity_end                                     = Vstall*0.1  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 30.            
    segment.altitude                                         = 0.0  
    segment.state.unknowns.velocity_x                        = 0.1* Vstall * ones_row(1)    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.3])    
    mission.append_segment(segment)   

            
    return mission 



# ------------------------------------------------------------------
#   Full Mission 1500
# ------------------------------------------------------------------
def direct_mission_setup_at_2000ft(analyses,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'direct_mission_setup_at_2000ft'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport      

    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976() 
    atmo_data          = atmosphere.compute_values(altitude = airport.altitude,temperature_deviation= 1.)     

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                             = Segments.Segment() 
    ones_row                                                 = base_segment.state.ones_row
    base_segment.battery_discharge                           = True  
    base_segment.state.numerics.number_control_points        = control_points 
    base_segment.state.numerics.discretization_method        = linear_data
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability     = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle        = analyses.base.aerodynamics.geometry
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)         

    # ------------------------------------------------------------------
    #   Takeoff
    # ------------------------------------------------------------------      
    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff"  
    segment.analyses.extend( analyses.base )
    segment.velocity_start                                   = Vstall*0.1  
    segment.velocity_end                                     = Vstall  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 10.            
    segment.altitude                                         = 0.0  
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle   
    segment.battery_energy                                   = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)         
    mission.append_segment(segment) 
     
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Departure_End_of_Runway'       
    segment.analyses.extend( analyses.base )           
    segment.altitude_start                                   = 0.0 * Units.feet
    segment.altitude_end                                     = 50.0 * Units.feet
    segment.air_speed_start                                  = Vstall  
    segment.air_speed_end                                    = Vstall*1.1        
    segment.climb_rate                                       = 600 * Units['ft/min'] 
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3])    
    mission.append_segment(segment)  
                
    #------------------------------------------------------------------
    #  Initial Climb Area Segment Flight 1  
    #------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Initial_CLimb_Area'  
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                   = 50.0 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1     
    segment.air_speed_end                                    = Vstall*1.2  
    segment.climb_rate                                       = 600 * Units['ft/min']   
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle        
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])  
    mission.append_segment(segment)  
     
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Climb'       
    segment.analyses.extend( analyses.base )         
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 2000 * Units.feet 
    segment.air_speed_start                                  = Vstall*1.2  
    segment.air_speed_end                                    = 150.* Units['mph']    
    segment.climb_rate                                       = 500* Units['ft/min']
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle       
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])   
    mission.append_segment(segment)  
    
    # ------------------------------------------------------------------
    #   Cruise Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag = 'Cruise'  
    segment.analyses.extend(analyses.base) 
    segment.air_speed                                        = 150.* Units['mph']        
    cruise_distance                                          = airport_geospacial_data.flight_range - 15.5415*Units.nmi  
    segment.altitude                                         = 2000 * Units.feet 
    segment.distance                                         = cruise_distance          
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.8]) 
    mission.append_segment(segment)   
 
    
    # ------------------------------------------------------------------
    #   Descent Segment Flight 1   
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Descent'   
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                   = 2000 * Units.feet  
    segment.altitude_end                                     = 1000.0 * Units.feet
    segment.air_speed_start                                  = 150.* Units['mph']    
    segment.air_speed_end                                    = Vstall*1.3     
    segment.climb_rate                                       = -300 * Units['ft/min']    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.8])  
    mission.append_segment(segment)  
     

    # ------------------------------------------------------------------
    #  Downleg_Altitude Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment) 
    segment.tag = 'Downleg' 
    segment.analyses.extend(analyses.base) 
    segment.air_speed_start                                  = Vstall*1.3 
    segment.air_speed_end                                    = Vstall*1.2       
    segment.distance                                         =  6000 * Units.feet
    segment.acceleration                                     = -0.05 * Units['m/s/s']     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.7])  
    mission.append_segment(segment)        
    
    # ------------------------------------------------------------------
    #  Baseleg Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = 'Baseleg' 
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                                   = 1000 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.2
    segment.air_speed_end                                    = Vstall*1.1  
    segment.climb_rate                                       = -300 * Units['ft/min']  
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment ,  initial_throttles = [0.7])
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #  Final Approach Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment_name = 'Final_Approach' 
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)            
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 00.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1  
    segment.air_speed_end                                    = Vstall 
    segment.climb_rate                                       = -300 * Units['ft/min'] 
    segment.state.unknowns.throttle                          =  0.8 * ones_row(1)    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3],  initial_throttles = [0.8] )  
    mission.append_segment(segment)   
    

    # ------------------------------------------------------------------
    #   Landing  
    # ------------------------------------------------------------------  
    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"   
    segment.analyses.extend( analyses.base) 
    segment.velocity_start                                   = Vstall  
    segment.velocity_end                                     = Vstall*0.1  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 30.            
    segment.altitude                                         = 0.0  
    segment.state.unknowns.velocity_x                        = 0.1* Vstall * ones_row(1)    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.3])    
    mission.append_segment(segment)   

            
    return mission 

# ------------------------------------------------------------------
#   Flyover at 200 ft 
# ------------------------------------------------------------------
def flyover_at_200ft_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'flyover_at_200ft'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0 
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.base )                  
    segment.altitude                                 = 200. * Units.ft 
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energy                           = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
    mission.append_segment(segment)     
             
    return mission 

# ------------------------------------------------------------------
#   Flyover at 500 ft 
# ------------------------------------------------------------------
def flyover_at_500ft_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'flyover_at_500ft'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0 
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.base )                  
    segment.altitude                                 = 500. * Units.ft 
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energy                           = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
    mission.append_segment(segment)     
             
    return mission 

# ------------------------------------------------------------------
#   Flyover at 1000 ft 
# ------------------------------------------------------------------
def flyover_at_1000ft_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'flyover_at_1000ft'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0 
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.base )                  
    segment.altitude                                 = 1000. * Units.ft 
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energy                           = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
    mission.append_segment(segment)     
             
    return mission 

# ------------------------------------------------------------------
#   Flyover at 1500 ft 
# ------------------------------------------------------------------
def flyover_at_1500ft_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'flyover_at_1500ft'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0 
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.base )                  
    segment.altitude                                 = 1000. * Units.ft  
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000   
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energy                           = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
    mission.append_segment(segment)     
             
    return mission 

# ------------------------------------------------------------------
#   Flyover at 2000 ft 
# ------------------------------------------------------------------
def flyover_at_2000ft_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'flyover_at_2000ft'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0 
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.base )                  
    segment.altitude                                 = 2000. * Units.ft  
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energy                           = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
    mission.append_segment(segment)     
             
    return mission 


# ------------------------------------------------------------------
#   Repeated Flight Operation Setup
# ------------------------------------------------------------------
def repeated_flight_operation_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission' 

    atmosphere         = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976() 
    atmo_data          = atmosphere.compute_values(altitude = 0 ,temperature_deviation= 1.)
        
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment() 
  

    # VSTALL Calculation  
    vehicle        = analyses.base.aerodynamics.geometry
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2) 
        
    bat                   = vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc
    Charging_C_Rate       = 1
    pack_charging_current = bat.cell.nominal_capacity * Charging_C_Rate *  bat.pack.electrical_configuration.series
    pack_capacity_Ah      =  bat.cell.nominal_capacity * bat.pack.electrical_configuration.series
    charging_time         = 0.9 * (pack_capacity_Ah)/pack_charging_current * Units.hrs  
    
    for day in range(simulated_days):
        
        # compute daily temperature in san francisco: link: https://www.usclimatedata.com/climate/san-francisco/california/united-states/usca0987/2019/1
        daily_temp = (13.5 + (day)*(-0.00882) + (day**2)*(0.00221) + (day**3)*(-0.0000314) + (day**4)*(0.000000185)  + \
                      (day**5)*(-0.000000000483)  + (day**6)*(4.57E-13)) + 273.2
        
        base_segment.temperature_deviation = daily_temp - atmo_data.temperature[0][0]
        
        print(' ***********  Day ' + str(day) + ' ***********  ')
        for f_idx in range(flights_per_day): 
            flight_no = f_idx + 1     
            ## ------------------------------------------------------------------
            ##   Takeoff Roll
            ## ------------------------------------------------------------------
        
            #segment = Segments.Ground.Takeoff(base_segment)
            #segment.tag = "Takeoff"  + "_F_" + str(flight_no) + "_D_" + str (day)  
            #segment.analyses.extend( analyses.base )
        
            #segment.velocity_start           = 100.* Units.knots #  Vstall*0.1  
            #segment.velocity_end             = 150 * Units.knots #  Vstall   
            #segment.friction_coefficient                             = 0.04 
            #segment.battery_pack_temperature                         = atmo_data.temperature[0,0]
            
            #segment.initial_battery_state_of_charge                  = 1.0
            #if (day == 0) and (f_idx == 0):        
                #segment.initial_battery_resistance_growth_factor     = 1
                #segment.initial_battery_capacity_fade_factor         = 1 
                   
            ## add to misison
            #mission.append_segment(segment) 
             
            # ------------------------------------------------------------------
            #   Departure End of Runway Segment Flight 1 : 
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = "DER"  + "_F_" + str(flight_no) + "_D_" + str (day)  
            segment.analyses.extend( analyses.max_hex_operation )  
            segment.altitude_start                                = 0.0 * Units.feet
            segment.altitude_end                                  = 50.0 * Units.feet
            segment.air_speed_start                               = 45  * Units['m/s'] 
            segment.air_speed_end                                 = 45
            

            segment.initial_battery_state_of_charge                  = 1.0
            if (day == 0) and (f_idx == 0):        
                segment.initial_battery_resistance_growth_factor     = 1
                segment.initial_battery_capacity_fade_factor         = 1             
                    
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True   
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                  
               
            mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #   Initial Climb Area Segment Flight 1  
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = 'Initial_CLimb_Area'  + "_F_" + str(flight_no) + "_D_" + str (day)   
            segment.analyses.extend( analyses.max_hex_operation )   
            segment.altitude_start                                = 50.0 * Units.feet
            segment.altitude_end                                  = 500.0 * Units.feet 
            segment.air_speed_end                                 = 50 * Units['m/s']   
            segment.climb_rate                                    = 600 * Units['ft/min']   
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                  
                  
            mission.append_segment(segment)  
           
                     
            # ------------------------------------------------------------------
            #   Climb Segment Flight 1 
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = 'Climb_1' + "_F_" + str(flight_no) + "_D_" + str (day)    
            segment.analyses.extend( analyses.hex_low_alt_climb_operation )      
            segment.altitude_start                                = 500.0 * Units.feet
            segment.altitude_end                                  = 2500 * Units.feet  
            segment.air_speed_end                                 = 120 * Units.kts 
            segment.climb_rate                                    = 500* Units['ft/min']  
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                 
                   
            mission.append_segment(segment)
            
                
            # ------------------------------------------------------------------
            #   Climb 1 : constant Speed, constant rate segment 
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment.tag = "Climb_2" + "_F_" + str(flight_no) + "_D_" + str (day)   
            segment.analyses.extend( analyses.hex_high_alt_climb_operation)
            segment.altitude_start                                = 2500.0  * Units.feet
            segment.altitude_end                                  = 5000   * Units.feet  
            segment.air_speed_end                                 = 130 * Units.kts 
            segment.climb_rate                                    = 700.034 * Units['ft/min']   
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                 
                    
            mission.append_segment(segment)
        
            # ------------------------------------------------------------------
            #   Cruise Segment: constant Speed, constant altitude
            # ------------------------------------------------------------------ 
            segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
            segment.tag = "Cruise"  + "_F_" + str(flight_no) + "_D_" + str (day)   
            segment.analyses.extend(analyses.hex_cruise_operation) 
            segment.altitude                                      = 5000   * Units.feet 
            segment.air_speed                                     = 130 * Units.kts
            segment.distance                                      = 50.   * Units.nautical_mile  
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                  
                  
            mission.append_segment(segment)    
        
        
            # ------------------------------------------------------------------
            #   Descent Segment Flight 1   
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = "Decent" + "_F_" + str(flight_no) + "_D_" + str (day)     
            segment.analyses.extend( analyses.hex_descent_operation )       
            segment.altitude_start                                = 5000   * Units.feet 
            segment.altitude_end                                  = 1000 * Units.feet  
            segment.air_speed_end                                 = 100 * Units['mph']   
            segment.climb_rate                                    = -200 * Units['ft/min']  
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                 
                  
            mission.append_segment(segment)   
                       
            # ------------------------------------------------------------------
            #  Downleg_Altitude Segment Flight 1 
            # ------------------------------------------------------------------
        
            segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
            segment.tag = 'Downleg' + "_F_" + str(flight_no) + "_D_" + str (day)   
            segment.analyses.extend(analyses.hex_cruise_operation)  
            segment.air_speed                                     = 100 * Units['mph']   
            segment.distance                                      = 6000 * Units.feet 
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                   
                    
            mission.append_segment(segment)     
            
            ## ------------------------------------------------------------------
            ##  Reserve Climb 
            ## ------------------------------------------------------------------ 
            #segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
            #segment.tag = 'Reserve_Climb' + "_F_" + str(flight_no) + "_D_" + str (day)           
            #segment.analyses.extend( analyses.hex_low_alt_climb_operation)      
            #segment.altitude_end                                  = 5000 * Units.feet
            #segment.air_speed                                     = 120 * Units['mph']
            #segment.climb_rate                                    = 500* Units['ft/min']  
            
            ## define flight dynamics to model 
            #segment.flight_dynamics.force_x                       = True  
            #segment.flight_dynamics.force_z                       = True     
            
            ## define flight controls 
            #segment.flight_controls.throttle.active               = True           
            #segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            #segment.flight_controls.body_angle.active             = True                
                
            #mission.append_segment(segment)
            
            ## ------------------------------------------------------------------
            ##  Researve Cruise Segment 
            ## ------------------------------------------------------------------ 
            #segment = Segments.Cruise.Constant_Speed_Constant_Altitude_Loiter(base_segment) 
            #segment.tag = 'Reserve_Cruise'  + "_F_" + str(flight_no) + "_D_" + str (day)    
            #segment.analyses.extend(analyses.hex_cruise_operation)  
            #segment.altitude                                      = 5000 * Units.feet
            #segment.air_speed                                     = 130 * Units.kts
            #segment.time                                          = 60*30 * Units.sec  
            
            ## define flight dynamics to model 
            #segment.flight_dynamics.force_x                       = True  
            #segment.flight_dynamics.force_z                       = True     
            
            ## define flight controls 
            #segment.flight_controls.throttle.active               = True           
            #segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            #segment.flight_controls.body_angle.active             = True                  
               
            #mission.append_segment(segment)     
            
            ## ------------------------------------------------------------------
            ##  Researve Descent
            ## ------------------------------------------------------------------ 
            #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment) 
            #segment.tag = 'Reserve_Descent' + "_F_" + str(flight_no) + "_D_" + str (day)   
            #segment.analyses.extend( analyses.hex_descent_operation)    
            #segment.altitude_end                                  = 1000 * Units.feet 
            #segment.air_speed                                     = 110 * Units['mph']
            #segment.descent_rate                                  = 300 * Units['ft/min']   
            
            ## define flight dynamics to model 
            #segment.flight_dynamics.force_x                       = True  
            #segment.flight_dynamics.force_z                       = True     
            
            ## define flight controls 
            #segment.flight_controls.throttle.active               = True           
            #segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            #segment.flight_controls.body_angle.active             = True                
            #mission.append_segment(segment)  
        
            
            ## ------------------------------------------------------------------
            ##  Baseleg Segment Flight 1  
            ## ------------------------------------------------------------------ 
            #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            #segment.tag = 'Baseleg' + "_F_" + str(flight_no) + "_D_" + str (day)   
            #segment.analyses.extend( analyses.hex_descent_operation)   
            #segment.altitude_start                                = 1000 * Units.feet
            #segment.altitude_end                                  = 500.0 * Units.feet
            #segment.air_speed_end                                 = 90 * Units['mph']  
            #segment.climb_rate                                    = -350 * Units['ft/min'] 
            
            ## define flight dynamics to model 
            #segment.flight_dynamics.force_x                       = True  
            #segment.flight_dynamics.force_z                       = True     
            
            ## define flight controls 
            #segment.flight_controls.throttle.active               = True           
            #segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            #segment.flight_controls.body_angle.active             = True                
            #mission.append_segment(segment) 
        
            ## ------------------------------------------------------------------
            ##  Final Approach Segment Flight 1  
            ## ------------------------------------------------------------------ 
            #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            #segment_name = 'Final_Approach' + "_F_" + str(flight_no) + "_D_" + str (day)   
            #segment.tag = segment_name          
            #segment.analyses.extend( analyses.hex_descent_operation)      
            #segment.altitude_start                                = 500.0 * Units.feet
            #segment.altitude_end                                  = 00.0 * Units.feet
            #segment.air_speed_end                                 = 80 * Units['mph']  
            #segment.climb_rate                                    = -300 * Units['ft/min']   
            
            ## define flight dynamics to model 
            #segment.flight_dynamics.force_x                       = True  
            #segment.flight_dynamics.force_z                       = True     
            
            ## define flight controls 
            #segment.flight_controls.throttle.active               = True           
            #segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            #segment.flight_controls.body_angle.active             = True                      
            #mission.append_segment(segment)  
        
        
            ## ------------------------------------------------------------------
            ##   Landing  
            ## ------------------------------------------------------------------  
            #segment = Segments.Ground.Landing(base_segment)
            #segment.tag = "Landing"  + "_F_" + str(flight_no) + "_D_" + str (day)     
            #segment.analyses.extend( analyses.hex_descent_operation)
        
            #segment.velocity_start                                = 150 * Units.knots# Vstall  
            #segment.velocity_end                                  = 100 * Units.knots# Vstall*0.1  
            #segment.friction_coefficient                          = 0.4   
            #segment.flight_controls.elapsed_time.active           = True  
            #segment.flight_controls.elapsed_time.initial_guess_values   = [[30.]] 
            #mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #   Mission definition complete    
            # ------------------------------------------------------------------            
  
            if recharge_battery:
                # ------------------------------------------------------------------
                #  Charge Segment: 
                # ------------------------------------------------------------------     
                # Charge Model 
                segment                                             = Segments.Ground.Battery_Recharge(base_segment)     
                segment.tag  = 'Charge_Day' + "_F_" + str(flight_no) + "_D" + str (day)  
                segment.analyses.extend(analyses.max_hex_operation)  
                segment.time                          = charging_time 
                segment.current                       = pack_charging_current          
                if flight_no  == flights_per_day:  
                    segment.increment_battery_age_by_one_day =True               
                mission.append_segment(segment)   
    
    
         
    return mission 


# ------------------------------------------------------------------
#   Constant Elevation Cruise Setup
# ------------------------------------------------------------------
def constant_elevation_in_cruise_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'constant_elevation_in_cruise_mission'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport      

    atmosphere         = RCAIDE.Analyses.Atmospheric.US_Standard_1976() 
    atmo_data          = atmosphere.compute_values(altitude = airport.altitude,temperature_deviation= 1.)     

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                             = Segments.Segment() 
    base_segment.battery_discharge                           = True  
    base_segment.state.numerics.number_control_points        = control_points
    base_segment.state.numerics.discretization_method        = linear_data
    ones_row                                                 = base_segment.state.ones_row
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability     = RCAIDE.Methods.skip


    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)    
     
      
    # ------------------------------------------------------------------
    #   Takeoff
    # ------------------------------------------------------------------      
    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff"  
    segment.analyses.extend( analyses.base )
    segment.velocity_start                                   = Vstall*0.1  
    segment.velocity_end                                     = Vstall  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 10.            
    segment.altitude                                         = 0.0  
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle   
    segment.battery_energy                                   = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)         
    mission.append_segment(segment) 
     
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Departure_End_of_Runway'       
    segment.analyses.extend( analyses.base )           
    segment.altitude_start                                   = 0.0 * Units.feet 
    segment.altitude_end                                     = 50.0 * Units.feet 
    segment.air_speed_start                                  = Vstall  
    segment.air_speed_end                                    = Vstall*1.1        
    segment.climb_rate                                       = 600 * Units['ft/min'] 
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3])    
    mission.append_segment(segment)  
                
    # ------------------------------------------------------------------
    #   Initial Climb Area Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Initial_CLimb_Area'  
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                   = 50.0 * Units.feet 
    segment.altitude_end                                     = 500.0 * Units.feet 
    segment.air_speed_start                                  = Vstall*1.1     
    segment.air_speed_end                                    = Vstall*1.2  
    segment.climb_rate                                       = 600 * Units['ft/min']  
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)    
    mission.append_segment(segment) 
            
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Climb'       
    segment.analyses.extend( analyses.base )         
    segment.altitude_start                                   = 500.0 * Units.feet 
    segment.altitude_end                                     = 1500 * Units.feet 
    segment.air_speed_start                                  = Vstall*1.2  
    segment.air_speed_end                                    = 130.* Units['mph']    
    segment.climb_rate                                       = 500* Units['ft/min']
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle       
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.85])   
    mission.append_segment(segment)  
    

    # ------------------------------------------------------------------
    #   Constant Ground Height Cruise (Climb Segment)
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Constant_Speed_Linear_Altitude(base_segment)
    segment.tag                                      = "Constant_Elevation_Cruise" 
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                           = 2500.0 * Units.ft   
    segment.altitude_end                             = 2500.0 * Units.ft     
    segment.air_speed                                = 130.  * Units['mph'] 
    cruise_distance                                  = airport_geospacial_data.flight_range - 31.5615*Units.nmi   
    segment.distance                                 = cruise_distance   
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])     
    mission.append_segment(segment)  
  
    # ------------------------------------------------------------------
    #   Descent Segment Flight 1   
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Descent'   
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                   = 1500 * Units.feet    
    segment.altitude_end                                     = 1000.0 * Units.feet   
    segment.air_speed_start                                  = 130.* Units['mph']    
    segment.air_speed_end                                    = Vstall*1.3
    segment.climb_rate                                       = -300 * Units['ft/min']
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle      
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.8])  
    mission.append_segment(segment)  
     

    # ------------------------------------------------------------------
    #  Downleg_Altitude Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment) 
    segment.tag = 'Downleg' 
    segment.analyses.extend(analyses.base) 
    segment.air_speed_start                                  = Vstall*1.3
    segment.air_speed_end                                    = Vstall*1.2             
    segment.distance                                         =  6000 * Units.feet
    segment.acceleration                                     = -0.05 * Units['m/s/s']  
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.7])  
    mission.append_segment(segment)        
    
    # ------------------------------------------------------------------
    #  Baseleg Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = 'Baseleg' 
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                                   = 1000 * Units.feet 
    segment.altitude_end                                     = 500.0 * Units.feet  
    segment.air_speed_start                                  = Vstall*1.2  
    segment.air_speed_end                                    = Vstall*1.1  
    segment.climb_rate                                       = -300 * Units['ft/min'] 
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment ,  initial_throttles = [0.7])
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #  Final Approach Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment_name = 'Final_Approach' 
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)            
    segment.altitude_start                                   = 500.0 * Units.feet 
    segment.altitude_end                                     = 00.0 * Units.feet  
    segment.air_speed_start                                  = Vstall*1.1  
    segment.air_speed_end                                    = Vstall 
    segment.climb_rate                                       = -300 * Units['ft/min'] 
    segment.state.unknowns.throttle                          =  0.8 * ones_row(1)   
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3],  initial_throttles = [0.8] )  
    mission.append_segment(segment)   
    

    # ------------------------------------------------------------------
    #   Landing  
    # ------------------------------------------------------------------  
    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"   
    segment.analyses.extend( analyses.base) 
    segment.velocity_start                                   = Vstall  
    segment.velocity_end                                     = Vstall*0.1  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 30.            
    segment.altitude                                         = 0.0  
    segment.state.unknowns.velocity_x                        = 0.1* Vstall * ones_row(1)   
    segment.true_course_angle                                = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.3])    
    mission.append_segment(segment)   
      
  
    return mission 
  

# ------------------------------------------------------------------
#   Approach and Departure Mission Setup
# ------------------------------------------------------------------
def approach_departure_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'approach_departure_mission'

    # airport
    airport           = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   =  0.0  * Units.ft
    airport.delta_isa  =  0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
   

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment() 
    base_segment.battery_discharge                                            = True  
    base_segment.state.numerics.number_control_points                         = control_points
    base_segment.state.numerics.discretization_method                         = linear_data
    ones_row                                                                  = base_segment.state.ones_row
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip


    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)    
  

    # ------------------------------------------------------------------
    #  Final Approach Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment_name = 'Final_Approach' 
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)            
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 00.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1  
    segment.air_speed_end                                    = Vstall 
    segment.climb_rate                                       = -300 * Units['ft/min'] 
    segment.state.unknowns.throttle                          =  0.8 * ones_row(1)     
    segment.battery_energy                                   = vehicle.networks.battery_electric_rotor.battery.pack.max_energy  
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3],  initial_throttles = [0.8] )  
    mission.append_segment(segment)   
    

    # ------------------------------------------------------------------
    #   Landing  
    # ------------------------------------------------------------------  
    segment = Segments.Ground.Landing(base_segment)
    segment.tag = "Landing"   
    segment.analyses.extend( analyses.base) 
    segment.velocity_start                                   = Vstall  
    segment.velocity_end                                     = Vstall*0.1  
    segment.friction_coefficient                             = 0.04 
    segment.state.unknowns.time                              = 30.            
    segment.altitude                                         = 0.0  
    segment.state.unknowns.velocity_x                        = 0.1* Vstall * ones_row(1)      
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,  initial_throttles = [0.3])    
    mission.append_segment(segment)  
    
    
    # ------------------------------------------------------------------
    #   Takeoff
    # ------------------------------------------------------------------      
    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff"  
    segment.analyses.extend( analyses.base )
    segment.velocity_start                                   = Vstall*0.01  
    segment.velocity_end                                     = Vstall  
    segment.friction_coefficient                             = 0.5
    segment.state.unknowns.time                              = 10.            
    segment.altitude                                         = 0.0   
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)         
    mission.append_segment(segment) 
     
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Departure_End_of_Runway'       
    segment.analyses.extend( analyses.base )           
    segment.altitude_start                                   = 0.0 * Units.feet
    segment.altitude_end                                     = 50.0 * Units.feet
    segment.air_speed_start                                  = Vstall  
    segment.air_speed_end                                    = Vstall*1.1        
    segment.climb_rate                                       = 600 * Units['ft/min']    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment, initial_rotor_power_coefficients  = [0.3])    
    mission.append_segment(segment)  
                
    # ------------------------------------------------------------------
    #   Initial Climb Area Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Initial_CLimb_Area'  
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                   = 50.0 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = Vstall*1.1     
    segment.air_speed_end                                    = Vstall*1.2  
    segment.climb_rate                                       = 600 * Units['ft/min']      
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)    
    mission.append_segment(segment)     
 

    return mission  
 

# ------------------------------------------------------------------
#   Range Mission Setup
# ------------------------------------------------------------------ 
def range_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
     
    mission = RCAIDE.Analyses.Mission.Variable_Range_Cruise.Given_State_of_Charge()
    mission.tag = 'Payload_Range'

    # the cruise tag to vary cruise distance
    mission.cruise_tag = 'cruise'
    mission.target_state_of_charge = 0.5

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment
    base_segment = Segments.Segment()
    ones_row                                                 = base_segment.state.ones_row    
    base_segment.state.numerics.number_control_points        = 4
    base_segment.state.numerics.discretization_method        = linear_data
    base_segment.process.iterate.conditions.stability        = RCAIDE.Methods.skip
    base_segment.process.finalize.post_process.stability     = RCAIDE.Methods.skip    
    base_segment.process.iterate.conditions.planet_position  = RCAIDE.Methods.skip    
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery

    # ------------------------------------------------------------------
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise"

    segment.analyses.extend(analyses.base)  
    segment.altitude  = 1500.0 * Units.ft
    segment.air_speed = 110.   * Units['mph']
    segment.distance  = 40.    * Units.miles     
    segment.battery_energy = base.networks.battery_electric_rotor.battery.max_energy
    segment.state.unknowns.throttle = 0.50 * ones_row(1)
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment) 

    mission.append_segment(segment)

    return mission 


def hover_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
     
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'uber_mission'
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment
    base_segment = Segments.Segment()
    base_segment.state.numerics.number_control_points    = 4
    base_segment.state.numerics.discretization_method        = linear_data
    base_segment.process.initialize.initialize_battery   = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    
    base_segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    base_segment.process.finalize.post_process.stability = RCAIDE.Methods.skip    

    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Hover.Climb(base_segment)
    segment.tag = "climb_1"

    segment.analyses.extend( analyses.hover_oei )

    segment.altitude_start  = 0.0  * Units.ft
    segment.altitude_end    = 40.  * Units.ft
    segment.climb_rate      = 200. * Units['ft/min']
    segment.battery_energy  = hover.networks.battery_electric_rotor.battery.max_energy  
    segment.process.iterate.unknowns.mission                 = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability             = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability          = RCAIDE.Methods.skip      
    
    segment = hover.networks.battery_electric_rotor.add_lift_unknowns_and_residuals_to_segment(segment)
    
    # add to misison
    mission.append_segment(segment)        

    return mission


# ------------------------------------------------------------------
#   Uber Mission Setup
# ------------------------------------------------------------------    
def uber_mission_setup(analyse,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
     
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'uber_mission'
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment
    base_segment = Segments.Segment()
    base_segment.state.numerics.number_control_points    = 4
    base_segment.state.numerics.discretization_method        = linear_data
    base_segment.process.initialize.initialize_battery   = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    
    base_segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    base_segment.process.finalize.post_process.stability = RCAIDE.Methods.skip    
    
    
    # VSTALL Calculation
    m      = base.mass_properties.max_takeoff
    g      = 9.81
    S      = base.reference_area
    atmo   = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    rho    = atmo.compute_values(1000.*Units.feet,0.).density
    CLmax  = 1.2
        
    Vstall = float(np.sqrt(2.*m*g/(rho*S*CLmax)))
    
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Hover.Climb(base_segment)
    segment.tag = "climb_1"

    segment.analyses.extend( analyses.base )

    segment.altitude_start  = 0.0  * Units.ft
    segment.altitude_end    = 40.  * Units.ft
    segment.climb_rate      = 200. * Units['ft/min']
    segment.battery_energy  = base.networks.battery_electric_rotor.battery.max_energy  
    segment.process.iterate.unknowns.mission                 = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability             = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability          = RCAIDE.Methods.skip      
    
    segment = base.networks.battery_electric_rotor.add_lift_unknowns_and_residuals_to_segment(segment)
    
    # add to misison
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2"

    segment.analyses.extend( analyses.base )

    segment.air_speed       = np.sqrt((500 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.altitude_start  = 40.0 * Units.ft
    segment.altitude_end    = 300. * Units.ft
    segment.climb_rate      = 500. * Units['ft/min']

    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Second Cruise Segment: Constant Speed, Constant Altitude
    # ------------------------------------------------------------------

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude_Loiter(base_segment)
    segment.tag = "Departure_Terminal_Procedures"

    segment.analyses.extend( analyses.base )

    segment.altitude  = 300.0 * Units.ft
    segment.time      = 60.   * Units.second
    segment.air_speed = 1.2*Vstall    
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Accelerated_Climb"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 300.0 * Units.ft
    segment.altitude_end    = 1500. * Units.ft
    segment.climb_rate      = 500.  * Units['ft/min']
    segment.air_speed_start = np.sqrt((500 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 110.  * Units['mph']              
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment,initial_prop_power_coefficient = 0.01)    
    
    # add to misison
    mission.append_segment(segment)    
    
    
    # ------------------------------------------------------------------
    #   Third Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "Cruise"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude  = 1500.0 * Units.ft
    segment.air_speed = 110.   * Units['mph']
    segment.distance  = 60.    * Units.miles              
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    
    
    # post-process aerodynamic derivatives in cruise
    #segment.process.finalize.post_process.aero_derivatives = RCAIDE.Methods.Flight_Dynamics.Static_Stability.compute_aero_derivatives    
        
    # add to misison
    mission.append_segment(segment)      
    
    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Decelerating_Descent"
    
    segment.analyses.extend( analyses.base )  
    segment.altitude_start  = 1500.0 * Units.ft
    segment.altitude_end    = 300. * Units.ft
    segment.climb_rate      = -500.  * Units['ft/min']
    segment.air_speed_start = 110.  * Units['mph']
    segment.air_speed_end   = 1.2*Vstall

    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Fourth Cruise Segment: Constant Speed, Constant Altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude_Loiter(base_segment)
    segment.tag = "Arrival_Terminal_Procedures"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude        = 300.   * Units.ft
    segment.air_speed       = 1.2*Vstall
    segment.time            = 60 * Units.seconds
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)        
        
    
    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_2"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 300.0 * Units.ft
    segment.altitude_end    = 40. * Units.ft
    segment.climb_rate      = -400.  * Units['ft/min']  # Uber has 500->300
    segment.air_speed_start = np.sqrt((400 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 1.2*Vstall                    
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)     
    
        
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Hover.Descent(base_segment)
    segment.tag = "descent_1"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 40.0  * Units.ft
    segment.altitude_end    = 40.   * Units.ft
    segment.descent_rate    = 300. * Units['ft/min']
    segment.process.iterate.unknowns.mission                 = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability             = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability          = RCAIDE.Methods.skip  
    
    segment = base.networks.battery_electric_rotor.add_lift_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #
    #   Reserve Mission
    #
    # ------------------------------------------------------------------   
    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Hover.Climb(base_segment)
    segment.tag = "reserve_climb_1"

    segment.analyses.extend( analyses.base )

    segment.altitude_start  = 0.0  * Units.ft
    segment.altitude_end    = 40.  * Units.ft
    segment.climb_rate      = 500. * Units['ft/min']
    segment.process.iterate.unknowns.mission                 = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability             = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability          = RCAIDE.Methods.skip  
    
    segment = base.networks.battery_electric_rotor.add_lift_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_climb_2"

    segment.analyses.extend( analyses.base )

    segment.air_speed       = np.sqrt((500 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.altitude_start  = 40.0 * Units.ft
    segment.altitude_end    = 300. * Units.ft
    segment.climb_rate      = 500. * Units['ft/min']

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip    

    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_Accelerated_Climb"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 300.0 * Units.ft
    segment.altitude_end    = 500. * Units.ft
    segment.climb_rate      = 500.  * Units['ft/min']
    segment.air_speed_start = np.sqrt((500 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 110.  * Units['mph']                                            

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip      
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    
    
    # add to misison
    mission.append_segment(segment)    
        
        
    # ------------------------------------------------------------------
    #   Third Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "reserve_Cruise"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude  = 500.0 * Units.ft
    segment.air_speed = 110.   * Units['mph']
    segment.distance  = 6.    * Units.miles                       

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip      
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    
    
    # add to misison
    mission.append_segment(segment)         
    
    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_Decelerating_Descent"
    
    segment.analyses.extend( analyses.base )  
    segment.altitude_start  = 500.0 * Units.ft
    segment.altitude_end    = 300. * Units.ft
    segment.climb_rate      = -500.  * Units['ft/min']
    segment.air_speed_start = 110.  * Units['mph']
    segment.air_speed_end   = 1.2*Vstall

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    
    
    # add to misison
    mission.append_segment(segment)      
    
    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_Descent_2"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 300.0 * Units.ft
    segment.altitude_end    = 40. * Units.ft
    segment.climb_rate      = -400.  * Units['ft/min']  # Uber has 500->300
    segment.air_speed_start = np.sqrt((400 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 1.2*Vstall                           

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip  
    
    segment = base.networks.battery_electric_rotor.add_cruise_unknowns_and_residuals_to_segment(segment)    
    
    # add to misison
    mission.append_segment(segment)       
    
    
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Hover.Descent(base_segment)
    segment.tag = "reserve_descent_1"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 40.0  * Units.ft
    segment.altitude_end    = 40.   * Units.ft
    segment.descent_rate    = 300. * Units['ft/min']
    segment.process.iterate.unknowns.mission                 = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability             = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability          = RCAIDE.Methods.skip      

    segment = base.networks.battery_electric_rotor.add_lift_unknowns_and_residuals_to_segment(segment)    

    
    # add to misison
    mission.append_segment(segment)        

    return mission

# ----------------------------------------------------------------------
#   Missions Setup
# ----------------------------------------------------------------------
def missions_setup(base_mission):

    # the mission container
    missions         = RCAIDE.Framework.Mission.Missions()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------

    missions.base = base_mission


    # done!
    return missions