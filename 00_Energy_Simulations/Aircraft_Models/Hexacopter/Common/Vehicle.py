''' 
# Vehicle.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import MARC
from MARC.Core import Units, Data   
import pickle
from MARC.Visualization.Performance.Aerodynamics.Vehicle import *  
from MARC.Visualization.Performance.Mission import *     
from MARC.Visualization.Performance.Energy.Battery import *   
from MARC.Visualization.Performance.Noise import *  
from MARC.Visualization.Geometry import *
from MARC.Components.Energy.Networks.Battery_Electric_Rotor                 import Battery_Electric_Rotor
from MARC.Methods.Power.Battery.Sizing                                      import initialize_from_mass 
from MARC.Methods.Geometry.Two_Dimensional.Planform                         import segment_properties
from MARC.Methods.Power.Battery.Sizing                                      import initialize_from_circuit_configuration 
from MARC.Methods.Weights.Correlations.Propulsion                           import nasa_motor
from MARC.Methods.Propulsion.electric_motor_sizing                          import size_optimal_motor
from MARC.Methods.Propulsion                                                import lift_rotor_design 
from MARC.Methods.Weights.Buildups.eVTOL.empty                              import empty
from MARC.Methods.Center_of_Gravity.compute_component_centers_of_gravity    import compute_component_centers_of_gravity
from MARC.Methods.Geometry.Two_Dimensional.Planform.wing_segmented_planform import wing_segmented_planform 
from MARC.Methods.Weights.Buildups.eVTOL.converge_evtol_weight              import converge_evtol_weight  
 
import os
import numpy as np 
from copy import deepcopy 

# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup(resize_aircraft,vehicle_name = 'Hexacopter_CRM') :
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------  
    
    if resize_aircraft:
         
         # ------------------------------------------------------------------
        #   Initialize the Vehicle
        # ------------------------------------------------------------------    
        vehicle                                     = MARC.Vehicle()
        vehicle.tag                                 = vehicle_name
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
                                                    
        wing                                        = MARC.Components.Wings.Main_Wing()  # this is the body of the vehicle 
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
        fuselage                                    = MARC.Components.Fuselages.Fuselage()
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
        segment                          = MARC.Components.Lofted_Body_Segment.Segment() 
        segment.tag                      = 'segment_0'   
        segment.percent_x_location       = 0.  
        segment.percent_z_location       = 0.0 
        segment.height                   = 0.1   
        segment.width                    = 0.1   
        fuselage.append_segment(segment)            
                                                    
        # Segment                                   
        segment                         = MARC.Components.Lofted_Body_Segment.Segment()
        segment.tag                     = 'segment_1'   
        segment.percent_x_location      = 0.200/4.
        segment.percent_z_location      = 0.1713/4.
        segment.height                  = 0.737
        segment.width                   = 1.2
        segment.vsp_data.top_angle      = 53.79 * Units.degrees 
        segment.vsp_data.bottom_angle   = 28.28 * Units.degrees     
        fuselage.append_segment(segment)            
                                                    
        # Segment                                   
        segment                         = MARC.Components.Lofted_Body_Segment.Segment()
        segment.tag                     = 'segment_2'   
        segment.percent_x_location      = 0.8251/4.
        segment.percent_z_location      = 0.2840/4.
        segment.height                  = 1.40 
        segment.width                   = 1.8
        segment.vsp_data.top_angle      = 0 * Units.degrees 
        segment.vsp_data.bottom_angle   = 0 * Units.degrees     
        fuselage.append_segment(segment)            
                                                    
        # Segment                                  
        segment                         = MARC.Components.Lofted_Body_Segment.Segment()
        segment.tag                     = 'segment_3'   
        segment.percent_x_location      = 3.342/4.
        segment.percent_z_location      = 0.356/4.
        segment.height                  = 1.40
        segment.width                   = 1.8
        #segment.vsp_data.top_angle      = 0 * Units.degrees 
        #segment.vsp_data.bottom_angle   = 0 * Units.degrees     
        fuselage.append_segment(segment)  
                                                    
        # Segment                                   
        segment                         = MARC.Components.Lofted_Body_Segment.Segment()
        segment.tag                     = 'segment_4'   
        segment.percent_x_location      = 3.70004/4.
        segment.percent_z_location      = 0.4636/4.
        segment.height                  = 0.9444
        segment.width                   = 1.2
        segment.vsp_data.top_angle      = -36.59 * Units.degrees 
        segment.vsp_data.bottom_angle   = -57.94 * Units.degrees 
        fuselage.append_segment(segment)             
        
        # Segment                                   
        segment                         = MARC.Components.Lofted_Body_Segment.Segment()
        segment.tag                     = 'segment_5'   
        segment.percent_x_location      = 1.
        segment.percent_z_location      = 0.6320/4.
        segment.height                  = 0.1    
        segment.width                   = 0.1    
        fuselage.append_segment(segment)             
        
                                                     
        # add to vehicle
        vehicle.append_component(fuselage)   
        

        # ------------------------------------------------------------------
        #   Nacelles
        # ------------------------------------------------------------------ 
        nacelle                = MARC.Components.Nacelles.Nacelle()
        nacelle.tag            = 'lift_rotor_nacelle'
        nacelle.length         = 0.5
        nacelle.diameter       = 0.8
        nacelle.inlet_diameter = 0.6
        nacelle.orientation_euler_angles  = [0,-90*Units.degrees,0.]    
        nacelle.flow_through   = True 
     
        lift_rotor_nacelle_origins   =  [[ -1.5,2.6,1.9],[ -1.5,-2.6,1.9],
                                    [2.5,6.0,1.9] ,[2.5,-6.,1.9],
                                    [6.5,2.6,1.9] ,[6.5,-2.6,1.9]]   
     
        for ii in range(6):
            rotor_nacelle          = deepcopy(nacelle)
            rotor_nacelle.tag      = 'engine_nacelle_' + str(ii+1) 
            rotor_nacelle.origin   = [lift_rotor_nacelle_origins[ii]]
            vehicle.append_component(rotor_nacelle)   
           
    
        # ------------------------------------------------------------------
        #   Nacelles
        # ------------------------------------------------------------------ 
        cowling                = MARC.Components.Nacelles.Nacelle()
        cowling.tag            = 'cowling'
        cowling.length         = 0.4
        cowling.diameter       = 2.6*2
        cowling.inlet_diameter = 2.55*2
        cowling.orientation_euler_angles  = [0,-90*Units.degrees,0.]    
        cowling.flow_through   = True 
     
        lift_rotor_nacelle_origins   =  [[ -1.5,2.6,1.8],[ -1.5,-2.6,1.8],
                                        [2.5,6.0,1.8] ,[2.5,-6.,1.8],
                                        [6.5,2.6,1.8] ,[6.5,-2.6,1.8]]  
        
     
        for ii in range(6):
            rotor_cowling          = deepcopy(cowling)
            rotor_cowling.tag      = 'cowling_' + str(ii+1) 
            rotor_cowling.origin   = [lift_rotor_nacelle_origins[ii]]
            vehicle.append_component(rotor_cowling)    
        
        #------------------------------------------------------------------
        # Network
        #------------------------------------------------------------------
        net                              = MARC.Components.Energy.Networks.Battery_Electric_Rotor() 
        net.rotor_group_indexes          = [0,0,0,0,0,0]
        net.motor_group_indexes          = [0,0,0,0,0,0]  
        net.esc_group_indexes            = [0,0,0,0,0,0]     
        net.active_propulsor_groups      = [True] 
        
        #------------------------------------------------------------------
        # Design Battery 
        #------------------------------------------------------------------
        bat                                   = MARC.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion_LiNiMnCoO2_18650() 
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
        lift_rotor_esc              = MARC.Components.Energy.Distributors.Electronic_Speed_Controller()
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
                      
        rotor                                   = MARC.Components.Energy.Converters.Lift_Rotor() 
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

        airfoil                                = MARC.Components.Airfoils.Airfoil()
        ospath    = os.path.abspath(__file__)
        separator = os.path.sep
        rel_path  = ospath.split( 'Hexacopter' + separator + 'Common')[0]  
        airfoil.coordinate_file                     =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
        airfoil.polar_files                         = [rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
        rotor.append_airfoil(airfoil)          
        rotor.airfoil_polar_stations            = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]    
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
        lift_rotor_motor                         = MARC.Components.Energy.Converters.Motor() 
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
                
    
        #------------------------------------------------------------------
        # Design Payload
        #------------------------------------------------------------------
        payload                       = MARC.Components.Energy.Peripherals.Avionics()
        payload.power_draw            = 0. 
        net.payload                   = payload
    
        #------------------------------------------------------------------
        # Design Avionics
        #------------------------------------------------------------------
        avionics            = MARC.Components.Energy.Peripherals.Avionics()
        avionics.power_draw = 200. * Units.watts
        net.avionics        = avionics
                                                    
        vehicle.append_component(net) 
        vehicle.wings['main_wing'].motor_spanwise_locations   = np.array(origins)[:,1]/ vehicle.wings['main_wing'].spans.projected
     

    
        settings = Data()
        converge_evtol_weight(vehicle,settings,contingency_factor = 1.0) 
        breakdown = empty(vehicle,settings,contingency_factor = 1.0 )
        print(breakdown)
        
        vehicle.weight_breakdown  = breakdown
        compute_component_centers_of_gravity(vehicle)
        vehicle.center_of_gravity() 
        
        save_aircraft_geometry(vehicle,vehicle.tag)        
        
    else: 
        vehicle = load_aircraft_geometry(vehicle_name) 
    
    return vehicle 
  
# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------
def configs_setup(vehicle):
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs = MARC.Components.Configs.Config.Container()

    base_config = MARC.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)
    
    # ------------------------------------------------------------------
    #   Hover Configuration
    # ------------------------------------------------------------------
    config = MARC.Components.Configs.Config(base_config)
    config.tag = 'hover'  
    configs.append(config)
    
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = MARC.Components.Configs.Config(base_config)
    config.tag = 'vertical_flight'    
    configs.append(config)
    
  
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = MARC.Components.Configs.Config(base_config)
    config.tag = 'vertical_transition'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     = 3.  * Units.degrees
    configs.append(config)
      
      
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = MARC.Components.Configs.Config(base_config)
    config.tag = 'descent_transition'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     = 3.  * Units.degrees
    configs.append(config)
  
    
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = MARC.Components.Configs.Config(base_config)
    config.tag = 'climb'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     =2.  * Units.degrees
    configs.append(config) 
    
    # ------------------------------------------------------------------
    #    Configuration
    # ------------------------------------------------------------------
    config = MARC.Components.Configs.Config(base_config)
    config.tag = 'forward_flight'   
    for rotor in config.networks.battery_electric_rotor.rotors:  
        rotor.inputs.pitch_command                     = 5.  * Units.degrees
    configs.append(config)     
    
    return configs


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