# Missions.py

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Methods.Performance.estimate_stall_speed        import estimate_stall_speed 

# ------------------------------------------------------------------
#   Repeated Flight Operation Setup
# ------------------------------------------------------------------
def repeated_flight_operation_setup(configs,analyses,day_group,g_idx,group,days_per_group, flights_per_day ,charge_throughput_nmc,cycle_day_nmc,resistance_growth_nmc,capacity_fade_nmc,charge_throughput_lfp,cycle_day_lfp,resistance_growth_lfp,capacity_fade_lfp):
 

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission' 

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment()
    base_segment.temperature_deviation  = 10
    base_segment.state.numerics.number_of_control_points  = 16
     # VSTALL Calculation  
    vehicle        = analyses.base.aerodynamics.vehicle
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)
  
 
    for d_idx in  range(days_per_group):
        day =  ((group -1) * days_per_group) +  d_idx + 1
        print(' ***********  Day ' + str(day) + ' ***********  ')
        for f_idx in range(flights_per_day):
            current_flight_no = f_idx + 1
            print(' ***********  Flight ' + str(current_flight_no) + ' ***********  ')
               
            # ------------------------------------------------------------------
            #   Departure End of Runway Segment Flight 1 : 
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = 'Departure_End_of_Runway'+ "_F_" + str(current_flight_no) + "_D_" + str (day)           
            segment.analyses.extend( analyses.lfp )  
            segment.altitude_start                                = 0.0 * Units.feet
            segment.altitude_end                                  = 50.0 * Units.feet
            segment.air_speed_start                               = Vstall *1.2  
            segment.air_speed_end                                 = Vstall *1.25
            segment.initial_battery_state_of_charge               = 1.0
            if group != 1 and d_idx == 0:
                try:
                    segment.charge_throughput            = Data()
                    segment.charge_throughput.power_bus  = charge_throughput_lfp[str(group-1)][0]
                    segment.charge_throughput.cruise_bus = charge_throughput_nmc[str(group-1)][0]
                    segment.resistance_growth            = Data()
                    segment.resistance_growth.power_bus  = resistance_growth_lfp[str(group-1)]
                    segment.resistance_growth.cruise_bus = resistance_growth_nmc[str(group-1)]
                    segment.capacity_fade                = Data()
                    segment.capacity_fade.power_bus      = capacity_fade_lfp[str(group-1)]
                    segment.capacity_fade.cruise_bus     = capacity_fade_nmc[str(group-1)]
                    segment.cycle_day                    = Data()
                    segment.cycle_day.power_bus          = cycle_day_lfp[str(group-1)]
                    segment.cycle_day.cruise_bus         = cycle_day_nmc[str(group-1)]
                except KeyError:
                    raise Exception(f"Error: The key '{group-1}' was not found in charge_throughput or cycle_day. Run the simulation for the previous day group.")
                        
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                  
               
            mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #   Initial Climb Area Segment Flight 1  
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = 'Initial_CLimb_Area'+ "_F_" + str(current_flight_no) + "_D_" + str (day)     
            segment.analyses.extend( analyses.lfp )     
            segment.altitude_start                                = 50.0 * Units.feet
            segment.altitude_end                                  = 500.0 * Units.feet 
            segment.air_speed_end                                 = Vstall *1.3 
            segment.climb_rate                                    = 600 * Units['ft/min']   
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                  
                  
            mission.append_segment(segment)  
           
                     
            # ------------------------------------------------------------------
            #   Climb Segment Flight 1 
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = 'Climb_1'+ "_F_" + str(current_flight_no) + "_D_" + str (day)            
            segment.analyses.extend( analyses.lfp )   
            segment.altitude_start                                = 500.0 * Units.feet
            segment.altitude_end                                  = 2500 * Units.feet   
            segment.air_speed_end                                 = 120 * Units.kts  
            segment.climb_rate                                    = 500* Units['ft/min']  
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                 
                   
            mission.append_segment(segment)
            
                
            # ------------------------------------------------------------------
            #   Climb 1 : constant Speed, constant rate segment 
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment.tag = "Climb_2"+ "_F_" + str(current_flight_no) + "_D_" + str (day)    
            segment.analyses.extend( analyses.nmc )  
            segment.altitude_start                                = 2500.0  * Units.feet
            segment.altitude_end                                  = 5000   * Units.feet  
            segment.air_speed_end                                 = 130 * Units.kts 
            segment.climb_rate                                    = 700.034 * Units['ft/min']   
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                 
                    
            mission.append_segment(segment)
        
            # ------------------------------------------------------------------
            #   Cruise Segment: constant Speed, constant altitude
            # ------------------------------------------------------------------ 
            segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
            segment.tag = "Cruise"+ "_F_" + str(current_flight_no) + "_D_" + str (day)    
            segment.analyses.extend( analyses.nmc )  
            segment.altitude                                      = 5000   * Units.feet 
            segment.air_speed                                     = 130 * Units.kts
            segment.distance                                      = 13.   * Units.nautical_mile  
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                  
                  
            mission.append_segment(segment)    
        
        
            # ------------------------------------------------------------------
            #   Descent Segment Flight 1   
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = "Decent" + "_F_" + str(current_flight_no) + "_D_" + str (day)     
            segment.analyses.extend( analyses.nmc )       
            segment.altitude_start                                = 5000   * Units.feet 
            segment.altitude_end                                  = 1000 * Units.feet  
            segment.air_speed_end                                 = 100 * Units['mph']   
            segment.climb_rate                                    = -200 * Units['ft/min']  
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                 
                  
            mission.append_segment(segment)   
                       
            # ------------------------------------------------------------------
            #  Downleg_Altitude Segment Flight 1 
            # ------------------------------------------------------------------
        
            segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
            segment.tag = 'Downleg'+ "_F_" + str(current_flight_no) + "_D_" + str (day)    
            segment.analyses.extend( analyses.nmc )  
            segment.air_speed                                     = 100 * Units['mph']   
            segment.distance                                      = 6000 * Units.feet 
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                   
                    
            mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #  Baseleg Segment Flight 1  
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment.tag = 'Baseleg'+ "_F_" + str(current_flight_no) + "_D_" + str (day)    
            segment.analyses.extend( analyses.lfp ) 
            segment.altitude_start                                = 1000 * Units.feet
            segment.altitude_end                                  = 500.0 * Units.feet
            segment.air_speed_end                                 = 90 * Units['mph']  
            segment.climb_rate                                    = -350 * Units['ft/min'] 
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                
            mission.append_segment(segment) 
        
            # ------------------------------------------------------------------
            #  Final Approach Segment Flight 1  
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = 'Final_Approach'+ "_F_" + str(current_flight_no) + "_D_" + str (day)    
            segment.analyses.extend( analyses.lfp )        
            segment.altitude_start                                = 500.0 * Units.feet
            segment.altitude_end                                  = 00.0 * Units.feet
            segment.air_speed_end                                 = 80 * Units['mph']  
            segment.climb_rate                                    = -300 * Units['ft/min']   
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.assigned_control_variables.throttle.active               = True           
            segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.assigned_control_variables.body_angle.active             = True                      
            mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #  Charge Segment: 
            # ------------------------------------------------------------------     
            # Charge Model 
            segment      = Segments.Ground.Battery_Recharge(base_segment)     
            segment.tag  = 'Charge_Day' + "_F_" + str(current_flight_no) + "_D_" + str (day)    
            segment.analyses.extend( analyses.base) 
            segment.cooling_time = 30 * Units.minutes
            segment.state.numerics.number_of_control_points = 64
            if f_idx ==  (flights_per_day - 1): 
                segment.increment_battery_age_by_one_day =  True 
                #segment.increment_battery_cycle_day      =  day-1
            mission.append_segment(segment)             
           
            
            # ------------------------------------------------------------------
            #   Mission definition complete    
            # ------------------------------------------------------------------       
            
         
    return mission 


# ----------------------------------------------------------------------
#   Missions Setup
# ----------------------------------------------------------------------
def missions_setup(mission):

    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 


    # done!
    return missions