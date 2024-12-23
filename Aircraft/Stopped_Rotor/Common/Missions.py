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
from RCAIDE.Core import Units    
import numpy as np 
from RCAIDE.Methods.Performance          import estimate_stall_speed  
from RCAIDE.Methods.Utilities.Chebyshev  import linear_data

# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def baseline_mission_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None ): 
    
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
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  

    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)      
     
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    segment     = Segments.Hover.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )  
    segment.altitude_start                          = 0.0  * Units.ft  
    segment.altitude_end                            = 200.  * Units.ft   
    segment.initial_battery_state_of_charge         = 1.0 
    segment.climb_rate                              = 500. * Units['ft/min'] 
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle  
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------
 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Vertical_Transition"  
    segment.analyses.extend( analyses.transition_flight )   
    segment.altitude                                = 200.  * Units.ft           
    segment.air_speed_start                         = 500. * Units['ft/min']
    segment.air_speed_end                           = 0.75 * Vstall
    segment.acceleration                            = 1.5
    segment.pitch_initial                           = 0.0 * Units.degrees
    segment.pitch_final                             = 2.  * Units.degrees   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle  
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Climb Transition 
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                     = "Climb_Transition" 
    segment.analyses.extend( analyses.transition_flight)    
    segment.altitude_start                          = 200.0 * Units.ft   
    segment.altitude_end                            = 500.0 * Units.ft   
    segment.air_speed_start                         = 0.75   * Vstall
    segment.climb_angle                             = 3     * Units.degrees   
    segment.acceleration                            = 0.25  * Units['m/s/s'] 
    segment.pitch_initial                           = 2.    * Units.degrees 
    segment.pitch_final                             = 7.    * Units.degrees  
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment) 
  
    # ------------------------------------------------------------------
    #   First Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_1"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 500.0 * Units.ft   
    segment.altitude_end                             = 1000. * Units.ft   
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_start                          = 104.1  * Units['mph']  
    segment.air_speed_end                            = 105.  * Units['mph']      
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #  Second Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_2"  
    segment.analyses.extend( analyses.forward_flight) 
    segment.altitude_start                           = 1000.0 * Units.ft   
    segment.altitude_end                             = 1500. * Units.ft   
    segment.climb_rate                               = 300.  * Units['ft/min']
    segment.air_speed_start                          = 105.  * Units['mph']  
    segment.air_speed_end                            = 110.  * Units['mph']     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,estimated_propulsor_group_throttles = [[0.8,0]])     
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Third Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise"  
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                 = 1500.0 * Units.ft  
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = airport_geospacial_data.flight_range  - 9.3*Units.nmi
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,estimated_propulsor_group_throttles = [[0.8,0]])      
    mission.append_segment(segment)   
    
    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                           = 1500.0 * Units.ft  
    segment.altitude_end                             = 1000. * Units.ft  
    segment.climb_rate                               = -500.  * Units['ft/min']
    segment.air_speed_start                          = 110.  * Units['mph']  
    segment.air_speed_end                            = Vstall     
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,estimated_propulsor_group_throttles = [[0.8,0]])      
    mission.append_segment(segment)  
     
            
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                      = "Descent_Transition"  
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude_start                           = 1000.0 * Units.ft   
    segment.altitude_end                             = 300.0 * Units.ft   
    segment.climb_angle                              = 3 * Units.degrees
    segment.acceleration                             = -0.25 * Units['m/s/s']    
    segment.pitch_initial                            = 4.3  * Units.degrees  
    segment.pitch_final                              = 7. * Units.degrees    
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle          
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Fifth Cuise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Decelerating_Transition"   
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude                                = 300.  * Units.ft     
    segment.air_speed_start                         = 36.5 * Units['mph'] 
    segment.air_speed_end                           = 300. * Units['ft/min'] 
    segment.acceleration                            = -0.25 * Units['m/s/s']    
    segment.pitch_initial                           = 7.  * Units.degrees  
    segment.pitch_final                             = 7. * Units.degrees     
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle  
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)       
    
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Hover.Descent(base_segment)
    segment.tag                                     = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                          = 300.0 * Units.ft   
    segment.altitude_end                            = 0.   * Units.ft  
    segment.descent_rate                            = 300. * Units['ft/min']   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)  
  

    # ------------------------------------------------------------------
    #  Charge Segment: 
    # ------------------------------------------------------------------  
    # Charge Model 
    segment                                         = Segments.Ground.Battery_Recharge(base_segment)     
    segment.tag                                     = 'Recharge' 
    segment.time                                    = 1 * Units.hr
    segment.analyses.extend(analyses.base)              
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)    
    return mission 
  
# ------------------------------------------------------------------
#   Repeated Flight Operation Setup
# ------------------------------------------------------------------
def repeated_flight_operation_setup(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'repeated_flight_operation_mission'

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
    base_segment.state.numerics.number_of_control_points        = control_points  
    base_segment.state.numerics.discretization_method        = linear_data    
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability     = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)      
    
    for day in range(simulated_days): 

        # compute daily temperature in san francisco: link: https://www.usclimatedata.com/climate/san-francisco/california/united-states/usca0987/2019/1
        daily_temp = (13.5 + (day)*(-0.00882) + (day**2)*(0.00221) + (day**3)*(-0.0000314) + (day**4)*(0.000000185)  + \
                      (day**5)*(-0.000000000483)  + (day**6)*(4.57E-13)) + 273.2
        
        base_segment.temperature_deviation = daily_temp - atmo_data.temperature[0][0]
        
                
        print(' ***********  Day ' + str(day) + ' ***********  ')
        for f_idx in range(flights_per_day):
        
            flight_no = f_idx + 1 
             
            # ------------------------------------------------------------------
            #   First Climb Segment: Constant Speed, Constant Rate
            # ------------------------------------------------------------------
            segment     = Segments.Hover.Climb(base_segment)
            segment.tag = "Vertical_Climb" + "_F_" + str(flight_no) + "_D" + str (day)   
            segment.analyses.extend( analyses.vertical_flight )  
            segment.altitude_start                          = 0.0  * Units.ft  
            segment.altitude_end                            = 200.  * Units.ft   
            segment.climb_rate                              = 500. * Units['ft/min'] 
            segment.true_course_angle                       = airport_geospacial_data.true_course_angle    
            if day == 0:        
                segment.battery_energies                    = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]   
            segment.battery_pack_temperature                = atmo_data.temperature[0,0]             
            segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
            segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
            segment.process.finalize.post_process.stability = RCAIDE.Methods.skip
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) # ,initial_rotor_power_coefficients = [0.02,0.05],initial_throttles =  [0.7,0.9] )
            mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #   First Cruise Segment: Transition
            # ------------------------------------------------------------------
         
            segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
            segment.tag                                     = "Vertical_Transition" + "_F_" + str(flight_no) + "_D" + str (day)  
            segment.analyses.extend( analyses.transition_flight ) 
            segment.altitude                                = 200.  * Units.ft           
            segment.air_speed_start                         = 500. * Units['ft/min']
            segment.air_speed_end                           = 0.8 * Vstall
            segment.acceleration                            = 2.5
            segment.pitch_initial                           = 0.0 * Units.degrees
            segment.pitch_final                             = 2.  * Units.degrees   
            segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
            segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
            segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
            segment.process.finalize.post_process.stability = RCAIDE.Methods.skip
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) #, initial_rotor_power_coefficients = [0.2,0.05],initial_throttles =  [0.9,0.9]  ) 
            mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #   Climb Transition 
            # ------------------------------------------------------------------ 
            segment                                         = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
            segment.tag                                     = "Climb_Transition" + "_F_" + str(flight_no) + "_D" + str (day) 
            segment.analyses.extend( analyses.transition_flight)    
            segment.altitude_start                          = 200.0 * Units.ft   
            segment.altitude_end                            = 500.0 * Units.ft   
            segment.air_speed                               = 0.8   * Vstall
            segment.climb_angle                             = 2     * Units.degrees   
            segment.acceleration                            = 0.25  * Units['m/s/s']    
            segment.pitch_initial                           = 7.    * Units.degrees  
            segment.pitch_final                             = 5.    * Units.degrees    
            segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
            segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
            segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
            segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) #,initial_rotor_power_coefficients = [0.1,0.05],initial_throttles =  [0.9,0.9]  ) 
            mission.append_segment(segment) 
          
            # ------------------------------------------------------------------
            #   First Climb
            # ------------------------------------------------------------------ 
            segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment.tag                                      = "Climb_1" + "_F_" + str(flight_no) + "_D" + str (day)   
            segment.analyses.extend( analyses.forward_flight ) 
            segment.altitude_start                           = 500.0 * Units.ft   
            segment.altitude_end                             = 7500. * Units.ft   
            segment.climb_rate                               = 500.  * Units['ft/min']
            segment.air_speed_start                          = 118.  * Units['mph']
            segment.air_speed_end                            = 125.  * Units['mph']      
            segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) #, initial_throttles = [0.7,0.0])   
            mission.append_segment(segment)  
        
            # ------------------------------------------------------------------
            #  Second Climb
            # ------------------------------------------------------------------ 
            segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment.tag                                      = "Climb_2" + "_F_" + str(flight_no) + "_D" + str (day)  
            segment.analyses.extend( analyses.forward_flight) 
            segment.altitude_start                           = 1000.0 * Units.ft   
            segment.altitude_end                             = 2500. * Units.ft   
            segment.climb_rate                               = 500.  * Units['ft/min']
            segment.air_speed_start                          = 125.  * Units['mph']  
            segment.air_speed_end                            = 130.  * Units['mph']     
            segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) #, initial_throttles = [0.8,0.0])     
            mission.append_segment(segment)  
        
            # ------------------------------------------------------------------
            #   Third Cruise Segment: Constant Acceleration, Constant Altitude
            # ------------------------------------------------------------------ 
            segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
            segment.tag                                      = "Cruise" + "_F_" + str(flight_no) + "_D" + str (day)  
            segment.analyses.extend( analyses.forward_flight )                  
            segment.altitude                                 = 2500.0 * Units.ft  
            segment.air_speed                                = 130.  * Units['mph'] 
            cruise_distance                                  = airport_geospacial_data.flight_range - 20.09 * Units.nmi
            segment.distance                                 = cruise_distance   
            segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
            mission.append_segment(segment)   
            
            # ------------------------------------------------------------------
            #   First Descent Segment: Constant Acceleration, Constant Rate
            # ------------------------------------------------------------------
        
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment.tag = "Descent"  + "_F_" + str(flight_no) + "_D" + str (day) 
            segment.analyses.extend(analyses.forward_flight)  
            segment.altitude_start                           = 2500.0 * Units.ft  
            segment.altitude_end                             = 1000. * Units.ft  
            segment.climb_rate                               = -500.  * Units['ft/min']
            segment.air_speed_start                          = 130.  * Units['mph']  
            segment.air_speed_end                            = 1.1*Vstall     
            segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
            mission.append_segment(segment)  
             
                    
            # ------------------------------------------------------------------
            #   First Cruise Segment: Transition
            # ------------------------------------------------------------------ 
            segment                                          = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
            segment.tag                                      = "Descent_Transition" + "_F_" + str(flight_no) + "_D" + str (day)   
            segment.analyses.extend( analyses.transition_flight ) 
            segment.altitude_start                           = 1000.0 * Units.ft   
            segment.altitude_end                             = 300.0 * Units.ft  
            segment.air_speed                                = 1.1*Vstall 
            segment.climb_angle                              = 3 * Units.degrees
            segment.acceleration                             = -0.25 * Units['m/s/s']    
            segment.pitch_initial                            = 2.3 * Units.degrees     
            segment.pitch_final                              = 7. * Units.degrees 
            segment.process.iterate.unknowns.mission         = RCAIDE.Methods.skip
            segment.process.iterate.conditions.stability     = RCAIDE.Methods.skip
            segment.process.finalize.post_process.stability  = RCAIDE.Methods.skip  
            segment.true_course_angle                        = airport_geospacial_data.true_course_angle          
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
            mission.append_segment(segment)  
        
            # ------------------------------------------------------------------
            #   Fifth Cuise Segment: Transition
            # ------------------------------------------------------------------ 
            segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
            segment.tag                                     = "Decelerating_Transition"  + "_F_" + str(flight_no) + "_D" + str (day)   
            segment.analyses.extend( analyses.transition_flight ) 
            segment.altitude                                = 300.  * Units.ft     
            segment.air_speed_start                         = 63. * Units['mph'] 
            segment.air_speed_end                           = 300. * Units['ft/min'] 
            segment.acceleration                            = -0.25 * Units['m/s/s']    
            segment.pitch_initial                           = 7.  * Units.degrees  
            segment.pitch_final                             = 7. * Units.degrees    
            segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
            segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
            segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
            segment.true_course_angle                       = airport_geospacial_data.true_course_angle    
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) #,initial_rotor_power_coefficients = [0.2,0.05],initial_throttles =  [0.3,0.9])            
            mission.append_segment(segment)       
            
            # ------------------------------------------------------------------
            #   Third Descent Segment: Constant Speed, Constant Rate
            # ------------------------------------------------------------------ 
            segment                                         = Segments.Hover.Descent(base_segment)
            segment.tag                                     = "Vertical_Descent" + "_F_" + str(flight_no) + "_D" + str (day) 
            segment.analyses.extend( analyses.vertical_flight)     
            segment.altitude_start                          = 300.0 * Units.ft   
            segment.altitude_end                            = 0.   * Units.ft  
            segment.descent_rate                            = 300. * Units['ft/min']   
            segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
            segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
            segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
            segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) #, initial_rotor_power_coefficients = [0.0,0.02], initial_throttles = [0.0,0.7])
            mission.append_segment(segment)  
           
            # ------------------------------------------------------------------
            #  Charge Segment: 
            # ------------------------------------------------------------------  
            # Charge Model 
            segment                                                  = Segments.Ground.Battery_Charge_Discharge(base_segment)     
            segment.tag                                              = 'Charge_Day_' + "_F_" + str(flight_no) + "_D" + str (day)  
            segment.analyses.extend(analyses.base)           
            segment.battery_discharge                                = False   
            if flight_no  == flights_per_day:  
                segment.increment_battery_cycle_day=True         
            segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
  
    return mission 
 
# ------------------------------------------------------------------
#   Cruise at 1000 ft 
# ------------------------------------------------------------------
def direct_mission_setup_at_1000ft(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'direct_mission_setup_at_1000ft'

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
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)      
     
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    segment     = Segments.Hover.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )    
    segment.altitude_start                                      = 0.0  * Units.ft    
    segment.altitude_end                                        = 200.  * Units.ft     
    segment.climb_rate                                          = 300. * Units['ft/min'] 
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle  
    segment.battery_energies                                      = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------
 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                                 = "Vertical_Transition"  
    segment.analyses.extend( analyses.transition_flight )  
    segment.altitude                                            = 200.  * Units.ft             
    segment.air_speed_start                                     = 500. * Units['ft/min']
    segment.air_speed_end                                       = 0.75 * Vstall
    segment.acceleration                                        = 1.5
    segment.pitch_initial                                       = 0.0 * Units.degrees
    segment.pitch_final                                         = 2.  * Units.degrees   
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Climb Transition 
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                                 = "Climb_Transition" 
    segment.analyses.extend( analyses.transition_flight)     
    segment.altitude_start                                      = 200.0 * Units.ft     
    segment.altitude_end                                        = 500.0 * Units.ft     
    segment.air_speed_start                                     = 0.75   * Vstall
    segment.climb_angle                                         = 3     * Units.degrees   
    segment.acceleration                                        = 0.25  * Units['m/s/s'] 
    segment.pitch_initial                                       = 2.    * Units.degrees 
    segment.pitch_final                                         = 7.    * Units.degrees  
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle      
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment) 
   
    # ------------------------------------------------------------------
    #   First Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_1"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 500.0 * Units.ft  
    segment.altitude_end                             = 750. * Units.ft      
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_end                            = 107   * Units['mph']     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #   Second Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_2"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 750.0 * Units.ft 
    segment.altitude_end                             = 1000. * Units.ft     
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_end                            = 110  * Units['mph']     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)       

    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------  
    cruise_distance                                  = airport_geospacial_data.flight_range - 11690.893101175745
    num_cruise_legs = int(np.ceil(cruise_distance/5000))
    for i in range(num_cruise_legs): 
        if i == num_cruise_legs-1:
            distance =  cruise_distance%5000
        else:
            distance = 5000
            
        segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
        segment.tag                                      = "Cruise_" + str(i+1)
        segment.analyses.extend( analyses.forward_flight )                   
        segment.altitude                                 = 1000. * Units.ft     
        segment.air_speed                                = 110.  * Units['mph']  
        segment.distance                                 = distance   
        segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
        segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
        mission.append_segment(segment)     
                 
            
    # ------------------------------------------------------------------
    #   Descent 1 
    # ------------------------------------------------------------------

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent_1"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                                      = 1000.0 * Units.ft   
    segment.altitude_end                                        = 750.   * Units.ft   
    segment.climb_rate                                          = -300.  * Units['ft/min']
    segment.air_speed_end                                       = 105.  * Units['mph']       
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)      
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #   Descent 2
    # ------------------------------------------------------------------            

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent_2"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                                      = 750.0 * Units.ft   
    segment.altitude_end                                        = 500.   * Units.ft    
    segment.climb_rate                                          = -300.  * Units['ft/min'] 
    segment.air_speed_end                                       = 100.  * Units['mph']       
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)      
    mission.append_segment(segment)   
                        
                        
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                                 = "Descent_Transition"  
    segment.analyses.extend( analyses.transition_flight )  
    segment.altitude_start                                      = 500.0 * Units.ft    
    segment.altitude_end                                        = 300.0 * Units.ft     
    segment.climb_angle                                         = 3 * Units.degrees
    segment.acceleration                                        = -0.5 * Units['m/s/s']    
    segment.pitch_initial                                       = 4.3  * Units.degrees  
    segment.pitch_final                                         = 7. * Units.degrees   
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip  
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle          
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Fifth Cuise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                                 = "Decelerating_Transition"   
    segment.analyses.extend( analyses.transition_flight )  
    segment.altitude                                            = 300.  * Units.ft       
    segment.air_speed_start                                     = 64.5 * Units['mph'] 
    segment.air_speed_end                                       = 300. * Units['ft/min'] 
    segment.acceleration                                        = -0.5 * Units['m/s/s']    
    segment.pitch_initial                                       = 7.  * Units.degrees  
    segment.pitch_final                                         = 7. * Units.degrees    
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip   
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle  
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)       
    
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Hover.Descent(base_segment)
    segment.tag                                                 = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_flight)      
    segment.altitude_start                                      = 300.0 * Units.ft      
    segment.altitude_end                                        = 0.   * Units.ft    
    segment.descent_rate                                        = 300. * Units['ft/min']   
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle    
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)  
     
    return mission 
  

# ------------------------------------------------------------------
#   Cruise at 1500 ft 
# ------------------------------------------------------------------
def direct_mission_setup_at_1500ft(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'direct_mission_setup_at_1500ft'

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
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)      
     
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    segment     = Segments.Hover.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )  
    segment.altitude_start                          = 0.0  * Units.ft    
    segment.altitude_end                            = 200.  * Units.ft     
    segment.climb_rate                              = 300. * Units['ft/min'] 
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle    
    segment.battery_energies                          = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------
 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Vertical_Transition"  
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude                                = 200.  * Units.ft             
    segment.air_speed_start                         = 500. * Units['ft/min']
    segment.air_speed_end                           = 0.75 * Vstall
    segment.acceleration                            = 1.5
    segment.pitch_initial                           = 0.0 * Units.degrees
    segment.pitch_final                             = 2.  * Units.degrees   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Climb Transition 
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                     = "Climb_Transition" 
    segment.analyses.extend( analyses.transition_flight)    
    segment.altitude_start                          = 200.0 * Units.ft     
    segment.altitude_end                            = 500.0 * Units.ft     
    segment.air_speed_start                         = 0.75   * Vstall
    segment.climb_angle                             = 3     * Units.degrees   
    segment.acceleration                            = 0.25  * Units['m/s/s'] 
    segment.pitch_initial                           = 2.    * Units.degrees 
    segment.pitch_final                             = 7.    * Units.degrees  
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment) 
    

    # ------------------------------------------------------------------
    #   First Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_1"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 500.0 * Units.ft  
    segment.altitude_end                             = 1000. * Units.ft      
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_end                            = 107   * Units['mph']     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #   Second Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_2"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 1000. * Units.ft 
    segment.altitude_end                             = 1500. * Units.ft     
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_end                            = 110  * Units['mph']     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)      

    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------  
    cruise_distance                                  = airport_geospacial_data.flight_range - 19243.78416114785  
    num_cruise_legs = int(np.ceil(cruise_distance/5000))
    for i in range(num_cruise_legs): 
        if i == num_cruise_legs-1:
            distance =  cruise_distance%5000
        else:
            distance = 5000
        segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
        segment.tag                                      = "Cruise_" + str(i+1)
        segment.analyses.extend( analyses.forward_flight )                  
        segment.altitude                                 = 1500. * Units.ft     
        segment.air_speed                                = 110.  * Units['mph']  
        segment.distance                                 = distance
        segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
        segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
        mission.append_segment(segment)     
            
            
    # ------------------------------------------------------------------
    #   Descent 1 
    # ------------------------------------------------------------------

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent_1"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                                      = 1500.0 * Units.ft   
    segment.altitude_end                                        = 1000.   * Units.ft   
    segment.climb_rate                                          = -300.  * Units['ft/min']
    segment.air_speed_end                                       = 105.  * Units['mph']       
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)      
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #   Descent 2
    # ------------------------------------------------------------------            

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent_2"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                                      = 1000.0 * Units.ft   
    segment.altitude_end                                        = 500.   * Units.ft    
    segment.climb_rate                                          = -300.  * Units['ft/min'] 
    segment.air_speed_end                                       = 100.  * Units['mph']       
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)      
    mission.append_segment(segment)   
                        
            
                        
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                      = "Descent_Transition"  
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude_start                           = 500.0 * Units.ft    
    segment.altitude_end                             = 300.0 * Units.ft     
    segment.climb_angle                              = 3 * Units.degrees
    segment.acceleration                             = -0.5 * Units['m/s/s']    
    segment.pitch_initial                            = 4.3  * Units.degrees  
    segment.pitch_final                              = 7. * Units.degrees   
    segment.process.iterate.unknowns.mission         = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability     = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability  = RCAIDE.Methods.skip  
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle          
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Fifth Cuise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Decelerating_Transition"   
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude                                = 300.  * Units.ft       
    segment.air_speed_start                         = 64.5 * Units['mph'] 
    segment.air_speed_end                           = 300. * Units['ft/min'] 
    segment.acceleration                            = -0.5 * Units['m/s/s']    
    segment.pitch_initial                           = 7.  * Units.degrees  
    segment.pitch_final                             = 7. * Units.degrees    
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle  
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)       
    
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Hover.Descent(base_segment)
    segment.tag                                     = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                          = 300.0 * Units.ft      
    segment.altitude_end                            = 0.   * Units.ft    
    segment.descent_rate                            = 300. * Units['ft/min']   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)

    return mission 



# ------------------------------------------------------------------
#   Cruise at 2000 ft 
# ------------------------------------------------------------------
def direct_mission_setup_at_2000ft(analyses,vehicle,simulated_days = 1,flights_per_day = 1,control_points = 10,recharge_battery = True,microphone_terrain_data = None ,airport_geospacial_data = None):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'direct_mission_setup_at_2000ft'

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
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)      
     
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    segment     = Segments.Hover.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )  
    segment.altitude_start                          = 0.0  * Units.ft    
    segment.altitude_end                            = 200.  * Units.ft     
    segment.climb_rate                              = 300. * Units['ft/min'] 
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle    
    segment.battery_energies                          = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------
 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Vertical_Transition"  
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude                                = 200.  * Units.ft             
    segment.air_speed_start                         = 500. * Units['ft/min']
    segment.air_speed_end                           = 0.75 * Vstall
    segment.acceleration                            = 1.5
    segment.pitch_initial                           = 0.0 * Units.degrees
    segment.pitch_final                             = 2.  * Units.degrees   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Climb Transition 
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                     = "Climb_Transition" 
    segment.analyses.extend( analyses.transition_flight)    
    segment.altitude_start                          = 200.0 * Units.ft     
    segment.altitude_end                            = 500.0 * Units.ft     
    segment.air_speed_start                         = 0.75   * Vstall
    segment.climb_angle                             = 3     * Units.degrees   
    segment.acceleration                            = 0.25  * Units['m/s/s'] 
    segment.pitch_initial                           = 2.    * Units.degrees 
    segment.pitch_final                             = 7.    * Units.degrees  
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   First Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_1"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 500.0 * Units.ft  
    segment.altitude_end                             = 1250. * Units.ft      
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_end                            = 107   * Units['mph']     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #   Second Climb
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                                      = "Climb_2"   
    segment.analyses.extend( analyses.forward_flight ) 
    segment.altitude_start                           = 1250. * Units.ft 
    segment.altitude_end                             = 2000. * Units.ft     
    segment.climb_rate                               = 500.  * Units['ft/min']
    segment.air_speed_end                            = 110  * Units['mph']     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)      

    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------  
    cruise_distance                                  = airport_geospacial_data.flight_range - 14.469047095637153 
    num_cruise_legs = int(np.ceil(cruise_distance/5000))
    for i in range(num_cruise_legs): 
        if i == num_cruise_legs-1:
            distance =  cruise_distance%5000
        else:
            distance = 5000
        segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
        segment.tag                                      = "Cruise_" + str(i+1)
        segment.analyses.extend( analyses.forward_flight )                  
        segment.altitude                                 = 2000. * Units.ft     
        segment.air_speed                                = 110.  * Units['mph']  
        segment.distance                                 = distance
        segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
        segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
        mission.append_segment(segment)     
            
            
    # ------------------------------------------------------------------
    #   Descent 1 
    # ------------------------------------------------------------------

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent_1"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                                      = 2000.0 * Units.ft   
    segment.altitude_end                                        = 1250.   * Units.ft   
    segment.climb_rate                                          = -300.  * Units['ft/min']
    segment.air_speed_end                                       = 105.  * Units['mph']       
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)      
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #   Descent 2
    # ------------------------------------------------------------------            

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Descent_2"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude_start                                      = 1250. * Units.ft   
    segment.altitude_end                                        = 500.   * Units.ft    
    segment.climb_rate                                          = -300.  * Units['ft/min'] 
    segment.air_speed_end                                       = 100.  * Units['mph']       
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)      
    mission.append_segment(segment)   
                         
                        
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                          = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                      = "Descent_Transition"  
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude_start                           = 500.0 * Units.ft    
    segment.altitude_end                             = 300.0 * Units.ft     
    segment.climb_angle                              = 3 * Units.degrees
    segment.acceleration                             = -0.5 * Units['m/s/s']    
    segment.pitch_initial                            = 4.3  * Units.degrees  
    segment.pitch_final                              = 7. * Units.degrees   
    segment.process.iterate.unknowns.mission         = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability     = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability  = RCAIDE.Methods.skip  
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle          
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Fifth Cuise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                     = "Decelerating_Transition"   
    segment.analyses.extend( analyses.transition_flight ) 
    segment.altitude                                = 300.  * Units.ft       
    segment.air_speed_start                         = 64.5 * Units['mph'] 
    segment.air_speed_end                           = 300. * Units['ft/min'] 
    segment.acceleration                            = -0.5 * Units['m/s/s']    
    segment.pitch_initial                           = 7.  * Units.degrees  
    segment.pitch_final                             = 7. * Units.degrees    
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle  
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)       
    
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                         = Segments.Hover.Descent(base_segment)
    segment.tag                                     = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_flight)     
    segment.altitude_start                          = 300.0 * Units.ft      
    segment.altitude_end                            = 0.   * Units.ft    
    segment.descent_rate                            = 300. * Units['ft/min']   
    segment.true_course_angle                       = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission        = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
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
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                 = 200. * Units.ft 
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energies                           = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
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
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                 = 500. * Units.ft 
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energies                           = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
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
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                 = 1000. * Units.ft 
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energies                           = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
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
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                 = 1500. * Units.ft  
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energies                           = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
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
    airport.altitude   = 0.0  * Units.ft
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 
 
    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                      = "Cruise" 
    segment.analyses.extend( analyses.forward_flight )                  
    segment.altitude                                 = 2000. * Units.ft  
    segment.air_speed                                = 110.  * Units['mph']  
    segment.distance                                 = 1000     
    segment.true_course_angle                        = airport_geospacial_data.true_course_angle   
    segment.battery_energies                           = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment, initial_throttles = [0.8,0.0])  
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
    airport            = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   = 0 
    airport.delta_isa  = 0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport       
    
    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment           
    base_segment                                                              = Segments.Segment()  
    base_segment.state.numerics.number_of_control_points                         = control_points 
    base_segment.state.numerics.discretization_method                         = linear_data
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health  
    base_segment.process.finalize.post_process.stability                      = RCAIDE.Methods.skip 

    # VSTALL Calculation  
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)     
                        
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                                 = "Descent_Transition"  
    segment.analyses.extend( analyses.transition_flight )  
    segment.battery_energies                                      = [vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.pack.maximum_energy]  
    segment.air_speed_start                                     = 100.  * Units['mph']    
    segment.altitude_start                                      = 500.0 * Units.ft    
    segment.altitude_end                                        = 300.0 * Units.ft     
    segment.climb_angle                                         = 3 * Units.degrees
    segment.acceleration                                        = -0.5 * Units['m/s/s']    
    segment.pitch_initial                                       = 4.3  * Units.degrees  
    segment.pitch_final                                         = 7. * Units.degrees   
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip          
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)  

    # ------------------------------------------------------------------
    #   Fifth Cuise Segment: Transition
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                                 = "Decelerating_Transition"   
    segment.analyses.extend( analyses.transition_flight )  
    segment.altitude                                            = 300.  * Units.ft       
    segment.air_speed_start                                     = 64.5 * Units['mph'] 
    segment.air_speed_end                                       = 300. * Units['ft/min'] 
    segment.acceleration                                        = -0.5 * Units['m/s/s']    
    segment.pitch_initial                                       = 7.  * Units.degrees  
    segment.pitch_final                                         = 7. * Units.degrees    
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip    
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)       
    
    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Hover.Descent(base_segment)
    segment.tag                                                 = "Vertical_Descent" 
    segment.analyses.extend( analyses.vertical_flight)      
    segment.altitude_start                                      = 300.0 * Units.ft      
    segment.altitude_end                                        = 0.   * Units.ft    
    segment.descent_rate                                        = 300. * Units['ft/min']      
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)  
    

    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    segment     = Segments.Hover.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )    
    segment.altitude_start                                      = 0.0  * Units.ft    
    segment.altitude_end                                        = 200.  * Units.ft     
    segment.climb_rate                                          = 300. * Units['ft/min'] 
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle    
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment) 
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Transition
    # ------------------------------------------------------------------
 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                                 = "Vertical_Transition"  
    segment.analyses.extend( analyses.transition_flight )  
    segment.altitude                                            = 200.  * Units.ft             
    segment.air_speed_start                                     = 500. * Units['ft/min']
    segment.air_speed_end                                       = 0.75 * Vstall
    segment.acceleration                                        = 1.5
    segment.pitch_initial                                       = 0.0 * Units.degrees
    segment.pitch_final                                         = 2.  * Units.degrees   
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle   
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Climb Transition 
    # ------------------------------------------------------------------ 
    segment                                                     = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    segment.tag                                                 = "Climb_Transition" 
    segment.analyses.extend( analyses.transition_flight)     
    segment.altitude_start                                      = 200.0 * Units.ft     
    segment.altitude_end                                        = 500.0 * Units.ft     
    segment.air_speed_start                                     = 0.75   * Vstall
    segment.climb_angle                                         = 3     * Units.degrees   
    segment.acceleration                                        = 0.25  * Units['m/s/s'] 
    segment.pitch_initial                                       = 2.    * Units.degrees 
    segment.pitch_final                                         = 7.    * Units.degrees  
    segment.true_course_angle                                   = airport_geospacial_data.true_course_angle      
    segment.process.iterate.unknowns.mission                    = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability                = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability             = RCAIDE.Methods.skip 
    segment = vehicle.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
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
    base_segment.state.numerics.number_of_control_points    = 4
    base_segment.state.numerics.discretization_method    = linear_data
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
    segment.battery_energies  = [hover.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.max_energy]  
    segment.process.iterate.unknowns.mission                 = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability             = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability          = RCAIDE.Methods.skip      
    
    segment = hover.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
    
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
    base_segment.state.numerics.number_of_control_points    = 4
    base_segment.state.numerics.discretization_method    = linear_data
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
    segment.battery_energies  = base.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.max_energy  
    segment.process.iterate.unknowns.mission                 = RCAIDE.Methods.skip
    segment.process.iterate.conditions.stability             = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability          = RCAIDE.Methods.skip      
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)
    
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
    segment.altitude_end    = 200.  * Units.ft
    segment.climb_rate      = 500. * Units['ft/min']

    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Second Cruise Segment: Constant Speed, Constant Altitude
    # ------------------------------------------------------------------

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude_Loiter(base_segment)
    segment.tag = "Departure_Terminal_Procedures"

    segment.analyses.extend( analyses.base )

    segment.altitude  = 200.0 * Units.ft
    segment.time      = 60.   * Units.second
    segment.air_speed = 1.2*Vstall    
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "Accelerated_Climb"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 200.0 * Units.ft
    segment.altitude_end    = 1500. * Units.ft
    segment.climb_rate      = 500.  * Units['ft/min']
    segment.air_speed_start = np.sqrt((500 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 110.  * Units['mph']              
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment,initial_prop_power_coefficient = 0.01)    
    
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
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
    
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
    segment.altitude_end    = 200.  * Units.ft
    segment.climb_rate      = -500.  * Units['ft/min']
    segment.air_speed_start = 110.  * Units['mph']
    segment.air_speed_end   = 1.2*Vstall

    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Fourth Cruise Segment: Constant Speed, Constant Altitude
    # ------------------------------------------------------------------
    
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude_Loiter(base_segment)
    segment.tag = "Arrival_Terminal_Procedures"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude        = 200.  * Units.ft
    segment.air_speed       = 1.2*Vstall
    segment.time            = 60 * Units.seconds
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)        
        
    
    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_2"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 200.0 * Units.ft
    segment.altitude_end    = 40. * Units.ft
    segment.climb_rate      = -400.  * Units['ft/min']  # Uber has 500->300
    segment.air_speed_start = np.sqrt((400 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 1.2*Vstall                    
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

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
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

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
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

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
    segment.altitude_end    = 200.  * Units.ft
    segment.climb_rate      = 500. * Units['ft/min']

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip    

    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

    # add to misison
    mission.append_segment(segment)    
    
    # ------------------------------------------------------------------
    #   Third Climb Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_Accelerated_Climb"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 200.0 * Units.ft
    segment.altitude_end    = 500. * Units.ft
    segment.climb_rate      = 500.  * Units['ft/min']
    segment.air_speed_start = np.sqrt((500 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 110.  * Units['mph']                                            

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip      
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
    
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
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
    
    # add to misison
    mission.append_segment(segment)         
    
    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_Decelerating_Descent"
    
    segment.analyses.extend( analyses.base )  
    segment.altitude_start  = 500.0 * Units.ft
    segment.altitude_end    = 200.  * Units.ft
    segment.climb_rate      = -500.  * Units['ft/min']
    segment.air_speed_start = 110.  * Units['mph']
    segment.air_speed_end   = 1.2*Vstall

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip 
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
    
    # add to misison
    mission.append_segment(segment)      
    
    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------
    
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = "reserve_Descent_2"
    
    segment.analyses.extend( analyses.base )
    
    segment.altitude_start  = 200.0 * Units.ft
    segment.altitude_end    = 40. * Units.ft
    segment.climb_rate      = -400.  * Units['ft/min']  # Uber has 500->300
    segment.air_speed_start = np.sqrt((400 * Units['ft/min'])**2 + (1.2*Vstall)**2)
    segment.air_speed_end   = 1.2*Vstall                           

    segment.process.iterate.conditions.stability    = RCAIDE.Methods.skip
    segment.process.finalize.post_process.stability = RCAIDE.Methods.skip  
    
    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    
    
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

    segment = base.networks.all_electric.add_unknowns_and_residuals_to_segment(segment)    

    
    # add to misison
    mission.append_segment(segment)        

    return mission



# ----------------------------------------------------------------------
#   Missions Setup
# ----------------------------------------------------------------------
def missions_setup(base_mission):

    # the mission container
    missions = RCAIDE.Analyses.Mission.Mission.Container()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------

    missions.base = base_mission


    # done!
    return missions  