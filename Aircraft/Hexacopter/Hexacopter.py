
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
import RCAIDE
from RCAIDE.Framework.Core import Units, Data    
from RCAIDE.Library.Methods.Energy.Sources.Batteries.Common                    import initialize_from_circuit_configuration 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion            import compute_motor_weight
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor                     import design_motor
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor                        import design_lift_rotor  
from RCAIDE.Library.Methods.Weights.Physics_Based_Buildups.Electric            import converge_physics_based_weight_buildup 
from RCAIDE.Library.Plots                                                      import *      
# python imports 
import os
import numpy as np 
from copy import deepcopy
import pickle
import  pandas as pd
import matplotlib.pyplot as plt  
 
# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def main():           
         
    # vehicle data
    new_geometry =  True 
    if new_geometry :
        vehicle  = vehicle_setup()
        save_aircraft_geometry(vehicle , 'Hexacopter')
    else: 
        vehicle = load_aircraft_geometry('Hexacopter')

    # plot vehicle 
    plot_3d_vehicle(vehicle, 
                    min_x_axis_limit            = -5,
                    max_x_axis_limit            = 15,
                    min_y_axis_limit            = -10,
                    max_y_axis_limit            = 10,
                    min_z_axis_limit            = -10,
                    max_z_axis_limit            = 10,
                    show_figure                 = False 
                    )           

    # Set up configs
    configs  = configs_setup(vehicle)

    # vehicle analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(analyses)
    missions = missions_setup(mission) 
     
    results = missions.base_mission.evaluate() 
     
    # plot the results 
    plot_results(results)    
     
    return
# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup() :
     
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    vehicle                                     = RCAIDE.Vehicle()
    vehicle.tag                                 = 'Hexacopter'
    vehicle.configuration                       = 'eVTOL'
    
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    
    # mass properties 
    vehicle.mass_properties.max_takeoff         = 3300 
    vehicle.mass_properties.takeoff             = vehicle.mass_properties.max_takeoff
    vehicle.mass_properties.operating_empty     = vehicle.mass_properties.max_takeoff 
    vehicle.mass_properties.center_of_gravity   = [[2.6, 0., 0. ] ] 
                                                
    # This needs updating                       
    vehicle.passengers                          = 6
    vehicle.reference_area                      = 73  * Units.feet**2 
    vehicle.envelope.ultimate_load              = 5.7   
    vehicle.envelope.limit_load                 = 3.  
                                                
    wing                                        = RCAIDE.Library.Components.Wings.Main_Wing()  # this is the body of the vehicle 
    wing.tag                                    = 'main_wing'   
    wing.aspect_ratio                           = 0.5 
    wing.sweeps.quarter_chord                   = 0.  
    wing.thickness_to_chord                     = 0.01   
    wing.spans.projected                        = 0.01  
    wing.chords.root                            = 0.01
    wing.total_length                           = 0.01
    wing.chords.tip                             = 0.01
    wing.chords.mean_aerodynamic                = 0.01
    wing.dihedral                               = 0.0  
    wing.areas.reference                        = 0.0001 
    wing.areas.wetted                           = 0.01
    wing.areas.exposed                          = 0.01  
    wing.symbolic                               = True 
    wing.symmetric                              = True 
    
    vehicle.append_component(wing)
    
    # ------------------------------------------------------    
    # FUSELAGE    
    # ------------------------------------------------------    
    # FUSELAGE PROPERTIES
    fuselage                                    = RCAIDE.Library.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage' 
    fuselage.seats_abreast                      = 2.  
    fuselage.seat_pitch                         = 3.  
    fuselage.fineness.nose                      = 0.88   
    fuselage.fineness.tail                      = 1.13   
    fuselage.lengths.nose                       = 0.5 
    fuselage.lengths.tail                       = 0.5
    fuselage.lengths.cabin                      = 4.
    fuselage.lengths.total                      = 5.
    fuselage.width                              = 1.8
    fuselage.heights.maximum                    = 1.8
    fuselage.heights.at_quarter_length          = 1.8
    fuselage.heights.at_wing_root_quarter_chord = 1.8
    fuselage.heights.at_three_quarters_length   = 1.8
    fuselage.areas.wetted                       = 19.829265
    fuselage.areas.front_projected              = 1.4294246 
    fuselage.effective_diameter                 = 1.4
    fuselage.differential_pressure              = 1. 
    
    # Segment  
    segment                          = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                      = 'segment_0'   
    segment.percent_x_location       = 0.  
    segment.percent_z_location       = 0.0 
    segment.height                   = 0.1   
    segment.width                    = 0.1   
    fuselage.append_segment(segment)            
                                                
    # Segment                                   
    segment                         = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                     = 'segment_1'   
    segment.percent_x_location      = 0.200/4.
    segment.percent_z_location      = 0.1713/4.
    segment.height                  = 0.737
    segment.width                   = 1.2
    segment.vsp_data.top_angle      = 53.79 * Units.degrees 
    segment.vsp_data.bottom_angle   = 28.28 * Units.degrees     
    fuselage.append_segment(segment)            
                                                
    # Segment                                   
    segment                         = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                     = 'segment_2'   
    segment.percent_x_location      = 0.8251/4.
    segment.percent_z_location      = 0.2840/4.
    segment.height                  = 1.40 
    segment.width                   = 1.8
    segment.vsp_data.top_angle      = 0 * Units.degrees 
    segment.vsp_data.bottom_angle   = 0 * Units.degrees     
    fuselage.append_segment(segment)            
                                                
    # Segment                                  
    segment                         = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                     = 'segment_3'   
    segment.percent_x_location      = 3.342/4.
    segment.percent_z_location      = 0.356/4.
    segment.height                  = 1.40
    segment.width                   = 1.8
    #segment.vsp_data.top_angle      = 0 * Units.degrees 
    #segment.vsp_data.bottom_angle   = 0 * Units.degrees     
    fuselage.append_segment(segment)  
                                                
    # Segment                                   
    segment                         = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                     = 'segment_4'   
    segment.percent_x_location      = 3.70004/4.
    segment.percent_z_location      = 0.4636/4.
    segment.height                  = 0.9444
    segment.width                   = 1.2
    segment.vsp_data.top_angle      = -36.59 * Units.degrees 
    segment.vsp_data.bottom_angle   = -57.94 * Units.degrees 
    fuselage.append_segment(segment)             
    
    # Segment                                   
    segment                         = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                     = 'segment_5'   
    segment.percent_x_location      = 1.
    segment.percent_z_location      = 0.6320/4.
    segment.height                  = 0.1    
    segment.width                   = 0.1    
    fuselage.append_segment(segment)             
    
                                                 
    # add to vehicle
    vehicle.append_component(fuselage)   
    
    
    

    #==================================================================================================================================== 
    # Lift Bus 
    #====================================================================================================================================          
    lift_bus                                               = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus()
    lift_bus.tag                                           = 'lift_bus' 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus Battery
    #------------------------------------------------------------------------------------------------------------------------------------ 
    bat                                                    = RCAIDE.Library.Components.Energy.Sources.Battery_Modules.Lithium_Ion_NMC()
    number_of_modules                                      = 1 
    bat.tag                                                = 'lift_bus_battery'
    bat.electrical_configuration.series                    = 140   
    bat.electrical_configuration.parallel                  = 20
    initialize_from_circuit_configuration(bat)  
    bat.geometrtic_configuration.total                      = bat.electrical_configuration.total
    bat.voltage                                             = bat.maximum_voltage 
    bat.geometrtic_configuration.normal_count               = 25
    bat.geometrtic_configuration.parallel_count             = 40 

    for _ in range(number_of_modules):
        lift_bus.battery_modules.append(deepcopy(bat))
        
    lift_bus.initialize_bus_electrical_properties()

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Lift Propulsors 
    #------------------------------------------------------------------------------------------------------------------------------------    
     
    # Define Lift Propulsor Container 
    lift_propulsor_1                                       = RCAIDE.Library.Components.Propulsors.Electric_Rotor()
    lift_propulsor_1.tag                                   = 'lift_propulsor_1'
    lift_propulsor_1.wing_mounted                          = True         
              
    # Electronic Speed Controller           
    lift_rotor_esc                                         = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller() 
    lift_rotor_esc.efficiency                              = 0.95    
    lift_rotor_esc.tag                                     = 'lift_rotor_esc_1' 
    lift_rotor_esc.origin                                  = [[-0.073 ,  1.950 , 1.2]] 
    lift_propulsor_1.electronic_speed_controller           = lift_rotor_esc 
           
    # Lift Rotor Design              
    lift_rotor                                             = RCAIDE.Library.Components.Propulsors.Converters.Lift_Rotor()   
    lift_rotor.tag                                         = 'lift_rotor_1'  
    lift_rotor.origin                                      = [[-0.073 ,  1.950 , 1.2]] 
    lift_rotor.active                                      = True          
    lift_rotor.orientation_euler_angles                    = [10.0*Units.degrees,np.pi/2.,0.]   # vector of angles defining default orientation of rotor
    lift_rotor.tip_radius                                  = 2.8/2
    lift_rotor.hub_radius                                  = 0.1 
    lift_rotor.number_of_blades                            = 3     
    lift_rotor.hover.design_altitude                       = 40 * Units.feet  
    lift_rotor.hover.design_thrust                         = Hover_Load/8
    lift_rotor.hover.design_freestream_velocity            = np.sqrt(lift_rotor.hover.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
    lift_rotor.oei.design_altitude                         = 40 * Units.feet  
    lift_rotor.oei.design_thrust                           = Hover_Load/7  
    lift_rotor.oei.design_freestream_velocity              = np.sqrt(lift_rotor.oei.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
    airfoil                                                = RCAIDE.Library.Components.Airfoils.Airfoil()   
    airfoil.coordinate_file                                = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                                    = [rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
    lift_rotor.append_airfoil(airfoil)                         
    lift_rotor.airfoil_polar_stations                      = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]    
    test_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../Tests/analysis_weights'))
    
    if new_regression:
        design_lift_rotor(lift_rotor)
        save_rotor(lift_rotor, os.path.join(test_dir, 'stopped_rotor_geometry.res'))
    else:
        regression_lift_rotor = deepcopy(lift_rotor)
        design_lift_rotor(regression_lift_rotor, iterations=2)
        loaded_lift_rotor = load_rotor(os.path.join(test_dir, 'stopped_rotor_geometry.res'))
        
        for key,item in lift_rotor.items():
            lift_rotor[key] = loaded_lift_rotor[key] 
        lift_rotor.Wake   = RCAIDE.Framework.Analyses.Propulsion.Rotor_Wake_Fidelity_Zero()         
            
    lift_propulsor_1.rotor =  lift_rotor          
    
    #------------------------------------------------------------------------------------------------------------------------------------               
    # Lift Rotor Motor  
    #------------------------------------------------------------------------------------------------------------------------------------    
    lift_rotor_motor                                       = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    lift_rotor_motor.efficiency                            = 0.9
    lift_rotor_motor.nominal_voltage                       = lift_bus.voltage*3/4  
    lift_rotor_motor.origin                                = [[-0.073 ,  1.950 , 1.2]]
    lift_rotor_motor.propeller_radius                      = lift_rotor.tip_radius
    lift_rotor_motor.tag                                   = 'lift_rotor_motor_1' 
    lift_rotor_motor.no_load_current                       = 0.01  
    lift_rotor_motor.wing_tag                              = 'main_wing'
    lift_rotor_motor.rotor_radius                          = lift_rotor.tip_radius
    lift_rotor_motor.design_torque                         = lift_rotor.hover.design_torque
    lift_rotor_motor.angular_velocity                      = lift_rotor.hover.design_angular_velocity/lift_rotor_motor.gear_ratio  
    design_motor(lift_rotor_motor)
    lift_rotor_motor.mass_properties.mass                  = compute_motor_weight(lift_rotor_motor.design_torque)     
    lift_propulsor_1.motor                                 = lift_rotor_motor
    


    #------------------------------------------------------------------------------------------------------------------------------------               
    # Lift Rotor Nacelle
    #------------------------------------------------------------------------------------------------------------------------------------     
    nacelle                           = RCAIDE.Library.Components.Nacelles.Nacelle()
    nacelle.tag                       = 'rotor_nacelle' 
    nacelle.length                    = 0.4
    nacelle.diameter                  = 2.6*2
    nacelle.inlet_diameter            = 2.55*2     
    nacelle.orientation_euler_angles  = [0,-90*Units.degrees,0.]    
    nacelle.flow_through   = True  
    lift_propulsor_1.nacelle          =  nacelle 
    lift_bus.propulsors.append(lift_propulsor_1)   


    #------------------------------------------------------------------------------------------------------------------------------------   
    # make and append copy of lift propulsor (efficient coding)    
    #------------------------------------------------------------------------------------------------------------------------------------   
    lift_propulsor_2                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_2.tag                                   = 'lift_propulsor_2' 
    lift_propulsor_2.rotor.origin                          = [[-0.073  , -1.950  , 1.2]]  
    lift_propulsor_2.rotor.orientation_euler_angle         = [-10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_2.motor.origin                          = [[-0.073  , -1.950  , 1.2]] 
    lift_propulsor_2.rotor.origin                          = [[-0.073  , -1.950  , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_2' 
    rotor_nacelle.origin                                   = [[ -0.073, -1.950, 1.2]]
    lift_propulsor_2.nacelle                               = rotor_nacelle  
    lift_bus.propulsors.append(lift_propulsor_2)    
        

    lift_propulsor_3                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_3.tag                                   = 'lift_propulsor_3' 
    lift_propulsor_3.rotor.origin                          = [[ 4.440 ,  1.950 , 1.2]]  
    lift_propulsor_3.rotor.orientation_euler_angle         = [10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_3.motor.origin                          = [[ 4.440 ,  1.950 , 1.2]] 
    lift_propulsor_3.rotor.origin                          = [[ 4.440 ,  1.950 , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_3'  
    rotor_nacelle.origin                                   = [[   4.413,   1.950 ,1.2]]
    lift_propulsor_3.nacelle                               = rotor_nacelle  
    lift_bus.propulsors.append(lift_propulsor_3)    
    

    lift_propulsor_4                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_4.tag                                   = 'lift_propulsor_4' 
    lift_propulsor_4.rotor.origin                          = [[ 4.440  , -1.950  , 1.2]]  
    lift_propulsor_4.rotor.orientation_euler_angle         = [-10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_4.motor.origin                          = [[ 4.440  , -1.950  , 1.2]] 
    lift_propulsor_4.rotor.origin                          = [[ 4.440  , -1.950  , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_4' 
    rotor_nacelle.origin                                   = [[   4.413, -1.950, 1.2]]
    lift_propulsor_4.nacelle                               = rotor_nacelle   
    lift_bus.propulsors.append(lift_propulsor_4)    
    

    lift_propulsor_5                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_5.tag                                   = 'lift_propulsor_5' 
    lift_propulsor_5.rotor.origin                          = [[ 0.219 ,  4.891 , 1.2]]  
    lift_propulsor_5.rotor.orientation_euler_angle         = [10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_5.motor.origin                          = [[ 0.219 ,  4.891 , 1.2]] 
    lift_propulsor_5.rotor.origin                          = [[ 0.219 ,  4.891 , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_5'  
    rotor_nacelle.origin                                   = [[   0.219 ,   4.891 , 1.2]] 
    lift_propulsor_5.nacelle                               = rotor_nacelle   
    lift_bus.propulsors.append(lift_propulsor_5)    
    

    lift_propulsor_6                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_6.tag                                   = 'lift_propulsor_6' 
    lift_propulsor_6.rotor.origin                          = [[ 0.219  , - 4.891 , 1.2]]  
    lift_propulsor_6.rotor.orientation_euler_angle         = [-10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_6.motor.origin                          = [[ 0.219  , - 4.891 , 1.2]] 
    lift_propulsor_6.rotor.origin                          = [[ 0.219  , - 4.891 , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_6'  
    rotor_nacelle.origin                                   = [[   0.219 , -  4.891 ,1.2]]
    lift_propulsor_6.nacelle                               = rotor_nacelle   
    lift_bus.propulsors.append(lift_propulsor_6)    
     
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Additional Bus Loads
    #------------------------------------------------------------------------------------------------------------------------------------            
    # Payload   
    payload                                                 = RCAIDE.Library.Components.Systems.Avionics()
    payload.power_draw                                      = 10. # Watts 
    payload.mass_properties.mass                            = 1.0 * Units.kg
    lift_bus.payload                                        = payload 
                             
    # Avionics                            
    avionics                                                = RCAIDE.Library.Components.Systems.Avionics()
    avionics.power_draw                                     = 10. # Watts  
    avionics.mass_properties.mass                           = 1.0 * Units.kg
    lift_bus.avionics                                       = avionics    

   
    network.busses.append(lift_bus)       
        
    # append energy network 
    vehicle.append_energy_network(network) 
     
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################   Determine Vehicle Mass Properties Using Physic Based Methods  ################################ 
    #------------------------------------------------------------------------------------------------------------------------------------   
    converged_vehicle, breakdown = converge_physics_based_weight_buildup(vehicle)  
    print(breakdown)
    
     
 
 
    cowling_origins   =  [[ -1.5,2.6,1.8],[ -1.5,-2.6,1.8],
                                    [2.5,6.0,1.8] ,[2.5,-6.,1.8],
                                    [6.5,2.6,1.8] ,[6.5,-2.6,1.8]]  
    
    
    
    #------------------------------------------------------------------
    # Design Battery 
    #------------------------------------------------------------------
    bat                                   = RCAIDE.Library.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion_LiNiMnCoO2_18650() 
    bat.pack.electrical_configuration.series                = 150   
    bat.pack.electrical_configuration.parallel              = 270
    initialize_from_circuit_configuration(bat,module_weight_factor = 1.05)    
    net.voltage                                             = bat.pack.max_voltage  
    bat.module_config.number_of_modules                     = 20  
    bat.module.geometrtic_configuration.total               = bat.pack.electrical_configuration.total
    bat.module_config.voltage                               = bat.pack.max_voltage/bat.module_config.number_of_modules # must be less than max_module_voltage (~50) /safety_factor (~ 1.5) 
    bat.module_config.layout_ratio                          = 0.5  
    bat.module.geometrtic_configuration.normal_count        = int(bat.module.geometrtic_configuration.total**(bat.module_config.layout_ratio))
    bat.module.geometrtic_configuration.parallel_count      = int(bat.module.geometrtic_configuration.total**(1-bat.module_config.layout_ratio)) 
    net.battery                                             = bat        
 
    # --------------------------------------------------------------
    # Lift Propulsor System 
    # -------------------------------------------------------------- 
    # 1. Electronic Speed Controller  
    lift_rotor_esc              = RCAIDE.Library.Components.Energy.Distributors.Electronic_Speed_Controller()
    lift_rotor_esc.efficiency   = 0.95
    for i in range(6):
        lift_rotor_ESC          = deepcopy(lift_rotor_esc)
        lift_rotor_ESC.tag      = 'lift_rotor_esc' + str(i + 1)  
        net.electronic_speed_controllers.append(lift_rotor_ESC) 
         
    # 2. Rotors 
    # atmosphere and flight conditions for propeller/lift_rotor design
    g                                       = 9.81                                   # gravitational acceleration  
    speed_of_sound                          = 340                                    # speed of sound 
    Hover_Load                              = vehicle.mass_properties.takeoff*g      # hover load   
    design_tip_mach                         = 0.7                                    # design tip mach number 
                  
    rotor                                   = RCAIDE.Library.Components.Energy.Converters.Lift_Rotor() 
    rotor.tip_radius                        = 2.5
    rotor.hub_radius                        = 0.15 * rotor.tip_radius 
    rotor.number_of_blades                  = 3

    rotor.hover.design_angular_velocity     = (design_tip_mach*speed_of_sound)/rotor.tip_radius    
    rotor.hover.design_altitude             = 500 * Units.feet                   
    rotor.hover.design_thrust               = Hover_Load/(6) 
    rotor.hover.design_freestream_velocity  = np.sqrt(rotor.hover.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2)))  

    rotor.oei.design_angular_velocity       = (design_tip_mach*speed_of_sound)/rotor.tip_radius    
    rotor.oei.design_altitude               = 500 * Units.feet                   
    rotor.oei.design_thrust                 = Hover_Load/(6-1)  
    rotor.oei.design_freestream_velocity    = np.sqrt(rotor.oei.design_thrust/(2*1.2*np.pi*(rotor.tip_radius**2)))
    
    rotor.variable_pitch                    = True  
    
    rotor                                  = lift_rotor_design(rotor)     
    
    # Appending rotors with different origins
    origins                 = [[ -1.5,2.6,1.7],[ -1.5,-2.6,1.7],
                                [2.5,6.0,1.7] ,[2.5,-6.,1.7],
                                [6.5,2.6,1.7] ,[6.5,-2.6,1.7]]   
    
    for ii in range(6):
        lift_rotor          = deepcopy(rotor)
        lift_rotor.tag      = 'mr_lift_rotor_' + str(ii+1)
        lift_rotor.origin   = [origins[ii]]
        net.rotors.append(lift_rotor)
    

    # Component 7: Motors
    lift_rotor_motor                         = RCAIDE.Library.Components.Energy.Converters.Motor() 
    lift_rotor_motor.efficiency              = 0.95
    lift_rotor_motor.nominal_voltage         = bat.pack.max_voltage*3/4   
    lift_rotor_motor.mass_properties.mass    = 3. * Units.kg  
    lift_rotor_motor.propeller_radius        = rotor.tip_radius   
    lift_rotor_motor.no_load_current         = 2.0      
    lift_rotor_motor.rotor_radius            = lift_rotor.tip_radius
    lift_rotor_motor.design_torque           = lift_rotor.hover.design_torque
    lift_rotor_motor.angular_velocity        = lift_rotor.hover.design_angular_velocity/lift_rotor_motor.gear_ratio  
    lift_rotor_motor                         = size_optimal_motor(lift_rotor_motor)
    lift_rotor_motor.mass_properties.mass    = nasa_motor(lift_rotor_motor.design_torque)  
 
     
    # Appending motors with different origins    
    for ii in range(6):
        lift_motor     = deepcopy(lift_rotor_motor)
        lift_motor.tag = 'motor_' + str(ii+1)
        lift_motor.origin    = [origins[ii]]
        net.motors.append(lift_motor)       
            





    
    return vehicle 
  
# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------
def configs_setup(vehicle):
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs = RCAIDE.Library.Components.Configs.Config.Container()

    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)
    
    # ------------------------------------------------------------------
    #   Hover Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'hover'  
    configs.append(config)
    
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'vertical_flight'    
    configs.append(config)
    
  
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'vertical_transition'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     = 3.  * Units.degrees
    configs.append(config)
      
      
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'descent_transition'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     = 3.  * Units.degrees
    configs.append(config)
  
    
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'climb'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     =2.  * Units.degrees
    configs.append(config) 
    
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'forward_flight'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     = 5.  * Units.degrees
    configs.append(config)     
    
    return configs

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE  

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ---------------------------------------------------------------------- 
def analyses_setup(configs,run_noise_analysis_flag,use_topology_flag,microphone_terrain_data,airport_geospacial_data):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config,run_noise_analysis_flag,use_topology_flag,microphone_terrain_data,airport_geospacial_data)
        analyses[tag] = analysis

    return analyses

# ------------------------------------------------------------------
# Base Analysis
# ------------------------------------------------------------------
def base_analysis(vehicle,run_noise_analysis_flag,use_topology_flag,microphone_terrain_data,airport_geospacial_data):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = RCAIDE.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)

    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle = vehicle  
    aerodynamics.settings.model_fuselage = True 
    aerodynamics.settings.number_spanwise_vortices           = 25
    aerodynamics.settings.number_chordwise_vortices          = 5    
    analyses.append(aerodynamics)  
        
    if run_noise_analysis_flag:  
        # ------------------------------------------------------------------
        #  Noise Analysis
        noise = RCAIDE.Analyses.Noise.Fidelity_Zero()   
        noise.geometry = vehicle

        # ------------------------------------------------------------------
        #  Noise Analysis
        noise = RCAIDE.Analyses.Noise.Fidelity_Zero()   
        noise.geometry = vehicle  
        noise.settings.mean_sea_level_altitude           = False 
        noise.settings.ground_microphone_x_resolution    = microphone_terrain_data.ground_microphone_x_resolution           
        noise.settings.ground_microphone_y_resolution    = microphone_terrain_data.ground_microphone_y_resolution             
        noise.settings.ground_microphone_min_x           = microphone_terrain_data.ground_microphone_min_x                 
        noise.settings.ground_microphone_max_x           = microphone_terrain_data.ground_microphone_max_x                 
        noise.settings.ground_microphone_min_y           = microphone_terrain_data.ground_microphone_min_y                 
        noise.settings.ground_microphone_max_y           = microphone_terrain_data.ground_microphone_max_y    
        noise.settings.ground_microphone_x_stencil       = microphone_terrain_data.ground_microphone_x_stencil             
        noise.settings.ground_microphone_y_stencil       = microphone_terrain_data.ground_microphone_y_stencil        
        
        if use_topology_flag:                             
            noise.settings.ground_microphone_locations   = microphone_terrain_data.cartesian_microphone_locations 
            noise.settings.aircraft_departure_location   = airport_geospacial_data.departure_location   
            noise.settings.aircraft_destination_location = airport_geospacial_data.destination_location   
        
        analyses.append(noise)                                                       
                                                                              
    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Analyses.Energy.Energy()
    energy.network = vehicle.networks
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    return analyses   


# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def mission_setup(analyses,vehicle,): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    starting_elevation  = 0 * Units.ft
    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'baseline_mission'

    # airport
    airport            = RCAIDE.Attributes.Airports.Airport() 
    airport.delta_isa  =  0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976() 
    mission.airport    = airport           

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment 
    base_segment                                             = Segments.Segment() 
    base_segment.state.numerics.number_of_control_points        = control_points    
    ones_row                                                 = base_segment.state.ones_row    
    base_segment.process.initialize.initialize_battery = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health   
    
           
            
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 

    segment                                            = Segments.Hover.Climb(base_segment)
    segment.tag                                        = "Vertical_Climb" 
    segment.analyses.extend( analyses.vertical_flight) 
    segment.altitude_start                             = 0.0  * Units.ft 
    segment.altitude_end                               = 200.  * Units.ft  
    segment.climb_rate                                 = 300. * Units['ft/min']  
    segment.battery_energy                             = vehicle.networks.battery_electric_rotor.battery.pack.max_energy    
    segment.true_course_angle                          = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #  First Transition Segment
    # ------------------------------------------------------------------  

    segment                                  = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag                              = "Vertical_Transition"  
    segment.analyses.extend( analyses.vertical_transition) 
    segment.altitude                         = 200.  * Units.ft       
    segment.air_speed_start                  = 300. * Units['ft/min'] 
    segment.air_speed_end                    = 35 * Units['mph']    
    segment.acceleration                     = 0.5
    segment.pitch_initial                    = 0. * Units.degrees
    segment.pitch_final                      = 0. * Units.degrees 
    segment.state.unknowns.throttle          = 0.8  * ones_row(1)  
    segment.true_course_angle                = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)

    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Climb_1"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 600. * Units['ft/min']
    segment.air_speed_start                  = 35.   * Units['mph']
    segment.air_speed_end                    = 55.  * Units['mph']       
    segment.altitude_start                   = 200.0 * Units.ft  
    segment.altitude_end                     = 500.0 * Units.ft    
    segment.true_course_angle                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Climb_2"  
    segment.analyses.extend(analyses.climb) 
    segment.climb_rate                       = 500. * Units['ft/min']
    segment.air_speed_start                  = 55.   * Units['mph']
    segment.air_speed_end                    = 75.  * Units['mph']       
    segment.altitude_start                   = 500.0 * Units.ft     
    segment.altitude_end                     = 2500.0 * Units.ft    
    segment.true_course_angle                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)                

    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                              = "Cruise"  
    segment.analyses.extend(analyses.forward_flight)  
    segment.altitude                         = 2500.0 * Units.ft      
    segment.air_speed                        = 75. * Units['mph']      
    segment.distance                         = airport_geospacial_data.flight_range  - 9.3*Units.nmi 
    segment.true_course_angle                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)      
    
                
    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag                              = "Descent_1"  
    segment.analyses.extend(analyses.forward_flight)
    segment.climb_rate                       = -200. * Units['ft/min']
    segment.air_speed_start                  = 75. * Units['mph']      
    segment.air_speed_end                    = 35. * Units['mph']      
    segment.altitude_start                   = 2500.0 * Units.ft 
    segment.altitude_end                     = 200.0 * Units.ft   
    segment.true_course_angle                = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)    
    mission.append_segment(segment)     
        
     
    # ------------------------------------------------------------------
    #  Third Transition Segment
    # ------------------------------------------------------------------

    segment                           = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment)
    segment.tag                       = "Decent_Transition" 
    segment.analyses.extend( analyses.descent_transition) 
    segment.altitude                  = 200.  * Units.ft 
    segment.air_speed_start           = 35.  * Units['mph'] 
    segment.air_speed_end             = 300. * Units['ft/min']
    segment.acceleration              = -0.5307 
    segment.pitch_initial             = 1. * Units.degrees
    segment.pitch_final               = 2. * Units.degrees        
    segment.true_course_angle         = airport_geospacial_data.true_course_angle    
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)

    
    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------ 
    segment                           = Segments.Hover.Descent(base_segment)
    segment.tag                       = "Vertical_Descent"  
    segment.analyses.extend( analyses.vertical_flight) 
    segment.altitude_start            = 200.0  * Units.ft  
    segment.altitude_end              = 0.  * Units.ft  
    segment.descent_rate              = 300. * Units['ft/min']  
    segment.true_course_angle         = airport_geospacial_data.true_course_angle     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)
             
  
    return mission 



def save_aircraft_geometry(geometry,filename): 
    pickle_file  = filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(geometry, file) 
    return 


def load_aircraft_geometry(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results