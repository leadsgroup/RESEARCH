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

import sys
sys.path.append('Common')  

import Vehicle
import Analyses

# ------------------------------------------------------------------
#   Repeated Flight Operation Setup
# ------------------------------------------------------------------
def repeated_flight_operation_setup(vehicle,simulated_days,flight_no, idle_time, taxi_time, cruise_distance, mean_temperature,ceiling_altitude, tms_operation, recharge_battery, airport, month):
    
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
    
    
    bat                   = vehicle.networks.all_electric.busses.bus.batteries.lithium_ion_nmc
    Charging_C_Rate       = 1
    pack_charging_current = bat.cell.nominal_capacity * Charging_C_Rate *  bat.pack.electrical_configuration.series
    pack_capacity_Ah      = bat.cell.nominal_capacity * bat.pack.electrical_configuration.parallel
    charging_time         = (pack_capacity_Ah)/pack_charging_current * Units.hrs
    
    
    for day in range(simulated_days):
                
        base_segment.temperature_deviation = (mean_temperature[0] + 274) - atmo_data.temperature[0][0]
        
        print(' *********** '+ str(airport) +' ' + str(month) + ' ***********  ')
        
        for f_idx in range(len(flight_no)): 
            current_flight_no = f_idx + 1
            print(' ***********  Flight ' + str(current_flight_no) + ' ***********  ')
            
            
            # -------------------------------------------------------------------------------------------    
            # SET UP CONFIGURATIONS 
            # -------------------------------------------------------------------------------------------
            configs           = Vehicle.configs_setup(vehicle, tms_operation, f_idx)  
            analyses          = Analyses.analyses_setup(configs)
            
            # VSTALL Calculation  
            vehicle        = analyses.base.aerodynamics.geometry
            vehicle_mass   = vehicle.mass_properties.max_takeoff
            reference_area = vehicle.reference_area 
            Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)            
            
            
            # ------------------------------------------------------------------
            #   Idle  cruise_distance[f_idx]
            # ------------------------------------------------------------------      
            segment = Segments.Ground.Idle(base_segment)
            segment.tag = "Idle"+ "_F_" + str(current_flight_no) + "_D" + str (day)    
            segment.analyses.extend(analyses.no_hex_operation)  
            segment.time                               = idle_time[f_idx]* Units.hr
            segment.initial_battery_state_of_charge    = 1.0  
            mission.append_segment(segment)
             
            # ------------------------------------------------------------------
            #   Taxi 
            # ------------------------------------------------------------------      
            segment = Segments.Ground.Taxi(base_segment)
            segment.tag = "Taxi"+ "_F_" + str(current_flight_no) + "_D" + str (day)    
            segment.analyses.extend( analyses.max_hex_operation)  
            segment.velocity                                      = 25 * Units.knots
            segment.time                                          = taxi_time[f_idx]* Units.mins
        
            segment.flight_dynamics.force_x                       = True   
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']]
            mission.append_segment(segment)  
            
            # ------------------------------------------------------------------
            #   Takeoff
            # ------------------------------------------------------------------      
            segment = Segments.Ground.Takeoff(base_segment)
            segment.tag = "Takeoff"+ "_F_" + str(current_flight_no) + "_D" + str (day)    
            segment.analyses.extend( analyses.max_hex_operation ) 
            segment.velocity_end                                     = Vstall*1.2  
            segment.friction_coefficient                             = 0.04   
            segment.throttle                                         = 0.8   
            
            segment.flight_dynamics.force_x                           = True 
            segment.flight_controls.elapsed_time.active               = True         
           
            mission.append_segment(segment) 
          
            
            # ------------------------------------------------------------------
            #   Departure End of Runway Segment Flight 1 : 
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
            segment.tag = 'Departure_End_of_Runway'+ "_F_" + str(current_flight_no) + "_D" + str (day)         
            segment.analyses.extend( analyses.max_hex_operation )  
            segment.altitude_start                                = 0.0 * Units.feet
            segment.altitude_end                                  = 50.0 * Units.feet
            segment.air_speed_start                               = Vstall *1.2  
            segment.air_speed_end                                 = Vstall *1.25  
                    
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
            segment.tag = 'Initial_CLimb_Area' + "_F_" + str(current_flight_no) + "_D" + str (day)  
            segment.analyses.extend( analyses.max_hex_operation )   
            segment.altitude_start                                = 50.0 * Units.feet
            segment.altitude_end                                  = 500.0 * Units.feet 
            segment.air_speed_end                                 = Vstall *1.3 
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
            segment.tag = 'Climb_1' + "_F_" + str(current_flight_no) + "_D" + str (day)         
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
            segment.tag = "Climb_2"+ "_F_" + str(current_flight_no) + "_D" + str (day)  
            segment.analyses.extend( analyses.hex_high_alt_climb_operation)
            segment.altitude_start                                = 2500.0  * Units.feet
            segment.altitude_end                                  = ceiling_altitude[f_idx]   * Units.feet  
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
            segment.tag = "Cruise" + "_F_" + str(current_flight_no) + "_D" + str (day)  
            segment.analyses.extend(analyses.hex_cruise_operation) 
            segment.altitude                                      = ceiling_altitude[f_idx]    * Units.feet 
            segment.air_speed                                     = 130 * Units.kts
            segment.distance                                      = cruise_distance[f_idx]* Units.nautical_mile  
            
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
            segment.tag = "Decent" + "_F_" + str(current_flight_no) + "_D" + str (day)   
            segment.analyses.extend( analyses.hex_descent_operation )       
            segment.altitude_start                                = ceiling_altitude[f_idx]    * Units.feet 
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
            segment.tag = 'Downleg'+ "_F_" + str(current_flight_no) + "_D" + str (day)  
            segment.analyses.extend(analyses.hex_descent_operation)  
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
            
             
            # ------------------------------------------------------------------
            #  Baseleg Segment Flight 1  
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment.tag = 'Baseleg'+ "_F_" + str(current_flight_no) + "_D" + str (day)  
            segment.analyses.extend( analyses.hex_descent_operation)   
            segment.altitude_start                                = 1000 * Units.feet
            segment.altitude_end                                  = 500.0 * Units.feet
            segment.air_speed_end                                 = 90 * Units['mph']  
            segment.climb_rate                                    = -350 * Units['ft/min'] 
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                
            mission.append_segment(segment) 
        
            # ------------------------------------------------------------------
            #  Final Approach Segment Flight 1  
            # ------------------------------------------------------------------ 
            segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
            segment_name = 'Final_Approach'+ "_F_" + str(current_flight_no) + "_D" + str (day)  
            segment.tag = segment_name          
            segment.analyses.extend( analyses.hex_descent_operation)      
            segment.altitude_start                                = 500.0 * Units.feet
            segment.altitude_end                                  = 00.0 * Units.feet
            segment.air_speed_end                                 = 80 * Units['mph']  
            segment.climb_rate                                    = -300 * Units['ft/min']   
            
            # define flight dynamics to model 
            segment.flight_dynamics.force_x                       = True  
            segment.flight_dynamics.force_z                       = True     
            
            # define flight controls 
            segment.flight_controls.throttle.active               = True           
            segment.flight_controls.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
            segment.flight_controls.body_angle.active             = True                      
            mission.append_segment(segment)  
        
        
            # ------------------------------------------------------------------
            #   Landing  
            # ------------------------------------------------------------------  
            segment = Segments.Ground.Landing(base_segment)
            segment.tag = "Landing"+ "_F_" + str(current_flight_no) + "_D" + str (day)     
            segment.analyses.extend( analyses.hex_descent_operation)  
            segment.velocity_end                                     = Vstall*0.1  
            
            segment.flight_dynamics.force_x                           = True 
            segment.flight_controls.elapsed_time.active               = True
            
            mission.append_segment(segment)
            
            # ------------------------------------------------------------------
            #  Charge Segment: 
            # ------------------------------------------------------------------     
            # Charge Model 
            segment                               = Segments.Ground.Battery_Recharge(base_segment)     
            segment.tag  = 'Charge_Day'+ "_F_" + str(current_flight_no) + "_D" + str (day)  
            segment.state.numerics.number_of_control_points  = 64 
            segment.analyses.extend(analyses.recharge)  
            segment.time                          = charging_time 
            segment.current                       = pack_charging_current        
            mission.append_segment(segment)   
        
            
            # ------------------------------------------------------------------
            #   Mission definition complete    
            # ------------------------------------------------------------------       
    
         
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