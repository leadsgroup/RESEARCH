Vstall =  45
# ------------------------------------------------------------------
#   Takeoff
# ------------------------------------------------------------------      
#segment = Segments.Ground.Takeoff(base_segment)
#segment.tag = "Takeoff"+ "_F_" + str(current_flight_no) + "_D" + str (day)    
##segment.analyses.extend( analyses.max_hex_operation )
#segment.analyses.extend( analyses.base )
#segment.velocity_end                                     = Vstall*1.2  
#segment.friction_coefficient                             = 0.04   
#segment.throttle                                         = 0.8   
#segment.velocity_start                                   = 5* Units.kts
##segment.flight_dynamics.force_x                           = True 
##segment.assigned_control_variables.elapsed_time.active   = True


#if f_idx == 0 and  d_idx == 0:
    #segment.day               = day
    #segment.charge_throughput = cumulative_Q
    #segment.initial_battery_state_of_charge               = 1.0
#mission.append_segment(segment) 


# ------------------------------------------------------------------
#   Departure End of Runway Segment Flight 1 : 
# ------------------------------------------------------------------ 
segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
segment.tag = 'Departure_End_of_Runway'+ "_F_" + str(current_flight_no) + "_D" + str (day)         
#segment.analyses.extend( analyses.max_hex_operation )
segment.analyses.extend( analyses.base )
segment.altitude_start                                = 0.0 * Units.feet
segment.altitude_end                                  = 50.0 * Units.feet
segment.air_speed_start                               = Vstall *1.2  
segment.air_speed_end                                 = Vstall *1.25  
        
# define flight dynamics to model 
segment.flight_dynamics.force_x                       = True  
segment.flight_dynamics.force_z                       = True     

# define flight controls 
segment.assigned_control_variables.throttle.active               = True           
segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
segment.assigned_control_variables.body_angle.active             = True                  
   
if f_idx == 0 and  d_idx == 0:
    segment.day               = day
    segment.charge_throughput = cumulative_Q
    segment.initial_battery_state_of_charge               = 1.0
    
mission.append_segment(segment)

# ------------------------------------------------------------------
#   Initial Climb Area Segment Flight 1  
# ------------------------------------------------------------------ 
segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
segment.tag = 'Initial_CLimb_Area' + "_F_" + str(current_flight_no) + "_D" + str (day)  
#segment.analyses.extend( analyses.max_hex_operation )
segment.analyses.extend( analyses.base )
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
segment.tag = 'Climb_1' + "_F_" + str(current_flight_no) + "_D" + str (day)         
#segment.analyses.extend( analyses.hex_low_alt_climb_operation )
segment.analyses.extend( analyses.base )
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
segment.tag = "Climb_2"+ "_F_" + str(current_flight_no) + "_D" + str (day)  
#segment.analyses.extend( analyses.hex_high_alt_climb_operation)
segment.analyses.extend( analyses.base )
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
segment.tag = "Cruise" + "_F_" + str(current_flight_no) + "_D" + str (day)  
#segment.analyses.extend(analyses.hex_cruise_operation)
segment.analyses.extend( analyses.base )
segment.altitude                                      = 5000   * Units.feet 
segment.air_speed                                     = 130 * Units.kts
segment.distance                                      = 20.   * Units.nautical_mile 

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
segment.tag = "Decent" + "_F_" + str(current_flight_no) + "_D" + str (day)   
#segment.analyses.extend( analyses.hex_descent_operation )
segment.analyses.extend( analyses.base )
segment.altitude_start                                = 5000 * Units.feet 
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
segment.tag = 'Downleg'+ "_F_" + str(current_flight_no) + "_D" + str (day)  
segment.air_speed                                     = 100 * Units['mph']   
segment.distance                                      = 6000 * Units.feet 
# define flight dynamics to model 
segment.flight_dynamics.force_x                       = True  
segment.flight_dynamics.force_z                       = True     
segment.analyses.extend( analyses.base )
# define flight controls 
segment.assigned_control_variables.throttle.active               = True           
segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
segment.assigned_control_variables.body_angle.active             = True                   
        
mission.append_segment(segment)     

 
# ------------------------------------------------------------------
#  Baseleg Segment Flight 1  
# ------------------------------------------------------------------ 
segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
segment.tag = 'Baseleg'+ "_F_" + str(current_flight_no) + "_D" + str (day)    
segment.altitude_start                                = 1000 * Units.feet
segment.altitude_end                                  = 500.0 * Units.feet
segment.air_speed_end                                 = 90 * Units['mph']  
segment.climb_rate                                    = -350 * Units['ft/min'] 
segment.analyses.extend( analyses.base )
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
segment_name = 'Final_Approach'+ "_F_" + str(current_flight_no) + "_D" + str (day)  
segment.tag = segment_name          
#segment.analyses.extend( analyses.hex_descent_operation)
segment.analyses.extend( analyses.base )
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
#   Landing  
# ------------------------------------------------------------------  
segment = Segments.Ground.Landing(base_segment)
segment.tag = "Landing"+ "_F_" + str(current_flight_no) + "_D" + str (day)     
#segment.analyses.extend( analyses.hex_descent_operation)
segment.analyses.extend( analyses.base )
segment.velocity_end                                     = Vstall*0.1  

segment.flight_dynamics.force_x                           = True 
segment.assigned_control_variables.elapsed_time.active               = True

mission.append_segment(segment)
# ------------------------------------------------------------------
#  Charge Segment: 
# ------------------------------------------------------------------     
# Charge Model 
segment                               = Segments.Ground.Battery_Recharge(base_segment)     
segment.tag  = 'Charge_Day'+ "_F_" + str(current_flight_no) + "_D" + str (day)   
segment.analyses.extend( analyses.base )
if f_idx ==  (flights_per_day - 1): 
    segment.increment_battery_age_by_one_day =  True 
#segment.time    =  vehicle.networks.electric.busses.bus.charging_time
#segment.current = vehicle.networks.electric.busses.bus.charging_current

mission.append_segment(segment)   
