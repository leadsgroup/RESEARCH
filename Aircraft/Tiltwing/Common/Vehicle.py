''' 
# Vehicle.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Core import Units, Data   
import pickle 
from RCAIDE.Visualization.Performance.Aerodynamics.Vehicle import *  
from RCAIDE.Visualization.Performance.Mission import *   
from RCAIDE.Visualization.Performance.Energy.Battery import *   
from RCAIDE.Visualization.Performance.Noise import *  
from RCAIDE.Visualization.Geometry import *
from RCAIDE.Components.Energy.Networks.Battery_Electric_Rotor                 import Battery_Electric_Rotor 
from RCAIDE.Methods.Performance.estimate_cruise_drag                          import estimate_cruise_drag
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform                         import segment_properties
from RCAIDE.Methods.Power.Battery.Sizing                                      import initialize_from_circuit_configuration 
from RCAIDE.Methods.Weights.Correlations.Propulsion                           import nasa_motor
from RCAIDE.Methods.Propulsion.electric_motor_sizing                          import size_optimal_motor
from RCAIDE.Methods.Propulsion                                                import prop_rotor_design
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform.wing_planform           import wing_planform 
from RCAIDE.Methods.Weights.Buildups.eVTOL.empty                              import empty
from RCAIDE.Methods.Center_of_Gravity.compute_component_centers_of_gravity    import compute_component_centers_of_gravity
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform.wing_segmented_planform import wing_segmented_planform 
from RCAIDE.Methods.Weights.Buildups.eVTOL.converge_evtol_weight              import converge_evtol_weight  
 
import os
import numpy as np 
from copy import deepcopy 

# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup(resize_aircraft,vehicle_name = 'Tiltwing_CRM') :
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------  
    
    if resize_aircraft:
         
        vehicle               = RCAIDE.Vehicle()
        vehicle.tag           = vehicle_name
        vehicle.configuration = 'eVTOL'
         
        vehicle.mass_properties.max_takeoff         = 2200
        vehicle.mass_properties.takeoff             = vehicle.mass_properties.max_takeoff
        vehicle.mass_properties.operating_empty     = vehicle.mass_properties.max_takeoff
        vehicle.mass_properties.center_of_gravity   = [[ 2.0144,   0.  ,  0.]] 
        vehicle.passengers                          = 6
        vehicle.flight_envelope.ultimate_load              = 5.7
        vehicle.flight_envelope.positive_limit_load                 = 3.     
        
        # ------------------------------------------------------    
        # WINGS    
        # ------------------------------------------------------  
        wing                                        = RCAIDE.Components.Wings.Main_Wing()
        wing.tag                                    = 'canard_wing'  
        wing.aspect_ratio                           = 8.51507 
        wing.sweeps.quarter_chord                   = 0.0
        wing.thickness_to_chord                     = 0.16 
        wing.taper                                  = 1.  
        span                                        = 12
        wing.spans.projected                        = span 
        chord                                       = 0.85
        wing.chords.root                            = chord
        wing.total_length                           = chord
        wing.chords.tip                             = chord
        wing.chords.mean_aerodynamic                = chord
        wing.dihedral                               = 0.0  
        wing.areas.reference                        = wing.chords.root*wing.spans.projected 
        wing.areas.wetted                           = 2*wing.chords.root*wing.spans.projected*0.95  
        wing.areas.exposed                          = 2*wing.chords.root*wing.spans.projected*0.95 
        wing.twists.root                            = 4. * Units.degrees 
        wing.twists.tip                             = 2. * Units.degrees  
        wing.origin                                 = [[0.1,  0.0 , 0.0]]  
        wing.aerodynamic_center                     = [0., 0., 0.]     
        wing.winglet_fraction                       = 0.0  
        wing.symmetric                              = True    
        airfoil                                   = RCAIDE.Components.Airfoils.Airfoil() 
        ospath    = os.path.abspath(__file__)
        separator = os.path.sep
        rel_path  = os.path.dirname(ospath) + separator             
        airfoil.coordinate_file                   = rel_path + 'Airfoils' + separator + 'NACA_63_412.txt' 
        wing.append_airfoil(airfoil)
                                                    
        # add to vehicle   
        vehicle.append_component(wing)                            
                                                    
        wing                                        = RCAIDE.Components.Wings.Main_Wing()
        wing.tag                                    = 'main_wing'  
        wing.aspect_ratio                           = 8.51507 
        wing.sweeps.quarter_chord                   = 0.0
        wing.thickness_to_chord                     = 0.16 
        wing.taper                                  = 1.  
        wing.spans.projected                        = span 
        wing.chords.root                            = chord
        wing.total_length                           = chord 
        wing.chords.tip                             = chord
        wing.chords.mean_aerodynamic                = chord
        wing.dihedral                               = 0.0  
        wing.areas.reference                        = wing.chords.root*wing.spans.projected 
        wing.areas.wetted                           = 2*wing.chords.root*wing.spans.projected*0.95  
        wing.areas.exposed                          = 2*wing.chords.root*wing.spans.projected*0.95 
        wing.twists.root                            = 4. * Units.degrees 
        wing.twists.tip                             = 2. * Units.degrees  
        wing.origin                                 = [[ 5.138, 0.0  ,  1.323 ]]  # for images 1.54
        wing.aerodynamic_center                     = [0., 0., 0.]     
        wing.winglet_fraction                       = 0.0  
        wing.symmetric                              = True 
        wing.append_airfoil(airfoil)
    
        # compute reference properties 
        wing_planform(wing) 
        wing = wing_planform(wing)
        vehicle.reference_area = wing.areas.reference*2   
         
        # add to vehicle 
        vehicle.append_component(wing)      
        
        
        # ------------------------------------------------------    
        # FUSELAGE    
        # ------------------------------------------------------    
        # FUSELAGE PROPERTIES                       
        fuselage                                    = RCAIDE.Components.Fuselages.Fuselage()
        fuselage.tag                                = 'fuselage' 
        fuselage.seats_abreast                      = 0.  
        fuselage.seat_pitch                         = 1.  
        fuselage.fineness.nose                      = 1.5 
        fuselage.fineness.tail                      = 4.0 
        fuselage.lengths.nose                       = 1.7   
        fuselage.lengths.tail                       = 2.7 
        fuselage.lengths.cabin                      = 1.7  
        fuselage.lengths.total                      = 6.3  
        fuselage.width                              = 1.15  
        fuselage.heights.maximum                    = 1.7 
        fuselage.heights.at_quarter_length          = 1.2  
        fuselage.heights.at_wing_root_quarter_chord = 1.7  
        fuselage.heights.at_three_quarters_length   = 0.75 
        fuselage.areas.wetted                       = 12.97989862  
        fuselage.areas.front_projected              = 1.365211404  
        fuselage.effective_diameter                 = 1.318423736  
        fuselage.differential_pressure              = 0.  
        
        # Segment  
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment() 
        segment.tag                                 = 'segment_0'   
        segment.percent_x_location                  = 0.  
        segment.percent_z_location                  = 0.  
        segment.height                              = 0.09  
        segment.width                               = 0.23473  
        segment.length                              = 0.  
        segment.effective_diameter                  = 0. 
        fuselage.Segments.append(segment)             
        
        # Segment                                   
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_1'   
        segment.percent_x_location                  = 0.97675/6.1 
        segment.percent_z_location                  = 0.21977/6.1
        segment.height                              = 0.9027  
        segment.width                               = 1.01709  
        fuselage.Segments.append(segment)             
        
        
        # Segment                                   
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_2'    
        segment.percent_x_location                  = 1.93556/6.1 
        segment.percent_z_location                  = 0.39371/6.1
        segment.height                              = 1.30558   
        segment.width                               = 1.38871  
        fuselage.Segments.append(segment)             
        
        
        # Segment                                   
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_3'    
        segment.percent_x_location                  = 3.44137/6.1 
        segment.percent_z_location                  = 0.57143/6.1
        segment.height                              = 1.52588 
        segment.width                               = 1.47074 
        fuselage.Segments.append(segment)             
        
        # Segment                                   
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_4'   
        segment.percent_x_location                  = 4.61031/6.1
        segment.percent_z_location                  = 0.81577/6.1
        segment.height                              = 1.14788 
        segment.width                               = 1.11463  
        fuselage.Segments.append(segment)              
        
        # Segment                                   
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_5'   
        segment.percent_x_location                  = 0.9827
        segment.percent_z_location                  = 0.180
        segment.height                              = 0.6145
        segment.width                               = 0.3838
        fuselage.Segments.append(segment)            
        
        
        # Segment                                   
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_6'   
        segment.percent_x_location                  = 1. 
        segment.percent_z_location                  = 0.2058
        segment.height                              = 0.4
        segment.width                               = 0.25
        fuselage.Segments.append(segment)                 
        
        # add to vehicle
        vehicle.append_component(fuselage)    

        
        # Nacelles 
        nacelle                = RCAIDE.Components.Nacelles.Nacelle()
        nacelle.tag            = 'nacelle'
        nacelle.length         = 1.5
        nacelle.diameter       = 0.5
        nacelle.orientation_euler_angles  = [0.,0.,0.]    
        nacelle.flow_through   = False  
        
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_1'
        nac_segment.percent_x_location = 0.0  
        nac_segment.height             = 0.0
        nac_segment.width              = 0.0
        nacelle.append_segment(nac_segment)    
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_2'
        nac_segment.percent_x_location = 0.10  
        nac_segment.height             = 0.3
        nac_segment.width              = 0.3
        nacelle.append_segment(nac_segment)    
        
        
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_3'
        nac_segment.percent_x_location = 0.25  
        nac_segment.height             = 0.45
        nac_segment.width              = 0.45
        nacelle.append_segment(nac_segment)    
        
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_4'
        nac_segment.percent_x_location = 0.5 
        nac_segment.height             = 0.5
        nac_segment.width              = 0.5
        nacelle.append_segment(nac_segment)    
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_5'
        nac_segment.percent_x_location = 0.75
        nac_segment.height             = 0.45
        nac_segment.width              = 0.45
        nacelle.append_segment(nac_segment)        
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_5'
        nac_segment.percent_x_location = 0.9
        nac_segment.height             = 0.3
        nac_segment.width              = 0.3
        nacelle.append_segment(nac_segment)    
        
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_7'
        nac_segment.percent_x_location = 1.0  
        nac_segment.height             = 0.0
        nac_segment.width              = 0.0
        nacelle.append_segment(nac_segment)     
    
    
        prop_nacelle_origins = [[-0.5, 2.0, 0.0], [-0.5, 4.8, 0.0],[-0.5, -2.0, 0.0], [-0.5, -4.8, 0.0],\
                   [4.5, 2.0 ,1.4], [4.5, 4.8, 1.4],[4.5, -2.0, 1.4], [4.5, -4.8, 1.4]] 
        
        for ii in range(8):
            prop_nacelle          = deepcopy(nacelle)
            prop_nacelle.tag      = 'propeller_nacelle_' + str(ii+1) 
            prop_nacelle.origin   = [prop_nacelle_origins[ii]]
            vehicle.append_component(prop_nacelle)        
        
        #------------------------------------------------------------------
        # PROPULSOR
        #------------------------------------------------------------------ 
        net                              = Battery_Electric_Rotor()
        net.rotor_group_indexes          = [0,0,0,0,0,0,0,0]
        net.motor_group_indexes          = [0,0,0,0,0,0,0,0]  
        net.esc_group_indexes            = [0,0,0,0,0,0,0,0]     
        net.active_propulsor_groups      = [True]
        
        
        # Battery         
        bat                                                    = RCAIDE.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion_LiNiMnCoO2_18650()  
        bat.pack.electrical_configuration.series               = 140   
        bat.pack.electrical_configuration.parallel             = 100
        bat.module_config.number_of_modules                    = 14 
        bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
        bat.module_config.voltage                              = bat.pack.max_voltage/bat.module_config.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
        bat.module.geometrtic_configuration.normal_count       = 25
        bat.module.geometrtic_configuration.parallel_count     = 40
        net.battery                                            = bat    
        net.voltage                                            = bat.pack.max_voltage 
     
        # 1. Electronic Speed Controller  
        prop_rotor_esc              = RCAIDE.Components.Energy.Distributors.Electronic_Speed_Controller()
        prop_rotor_esc.efficiency   = 0.95
        for i in range(8):
            prop_rotor_ESC          = deepcopy(prop_rotor_esc)
            prop_rotor_ESC.tag      = 'prop_rotor_esc' + str(i + 1)  
            net.electronic_speed_controllers.append(prop_rotor_ESC)    
        
        # 2. Prop Rotors     
        g               = 9.81                                   # gravitational acceleration  
        Drag            = estimate_cruise_drag(vehicle,altitude = 1500. * Units.ft,speed= 130.* Units['mph'] ,lift_coefficient = 0.5 ,profile_drag = 0.06)
        Hover_Load      = vehicle.mass_properties.takeoff*g      # hover load          
     
        prop_rotor                                   = RCAIDE.Components.Energy.Converters.Prop_Rotor() 
        prop_rotor.tag                               = 'prop_rotor'     
        prop_rotor.tip_radius                        = 1.25
        prop_rotor.hub_radius                        = 0.15 * prop_rotor.tip_radius
        prop_rotor.number_of_blades                  = 3 
           
        prop_rotor.hover.design_altitude             = 40 * Units.feet  
        prop_rotor.hover.design_thrust               = Hover_Load/8 
        prop_rotor.hover.design_freestream_velocity  = np.sqrt(prop_rotor.hover.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2))) 
    
        prop_rotor.oei.design_altitude               = 40 * Units.feet  
        prop_rotor.oei.design_thrust                 = Hover_Load/6
        prop_rotor.oei.design_freestream_velocity    = np.sqrt(prop_rotor.oei.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2)))  
    
        # CRUISE                   
        prop_rotor.cruise.design_altitude             = 1500 * Units.feet                      
        prop_rotor.cruise.design_thrust               = Drag*1.1/8
        prop_rotor.cruise.design_freestream_velocity  = 130.* Units['mph']  
 
        ospath                                      = os.path.abspath(__file__)
        separator                                   = os.path.sep
        rel_path                                    = ospath.split( 'Tiltwing' + separator + 'Common')[0]
        airfoil                                     = RCAIDE.Components.Airfoils.Airfoil() 
        airfoil.coordinate_file                     = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
        airfoil.polar_files                         = [rel_path+ 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
        prop_rotor.append_airfoil(airfoil)          
        prop_rotor.airfoil_polar_stations            = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]    
        prop_rotor.variable_pitch                    = True  
          
        prop_rotor                                   = prop_rotor_design(prop_rotor)   
                
        # Prop Rotors                                          
        prop_rotor_origins           =  [[-0.3, 2.0, 0.0], [-0.3, 4.8, 0.0],[-0.3, -2.0, 0.0], [-0.3, -4.8, 0.0],\
                   [4.7, 2.0 ,1.4], [4.7, 4.8, 1.4],[4.7, -2.0, 1.4], [4.7, -4.8, 1.4]]  
        for ii in range(8):
            pr                        = deepcopy(prop_rotor)
            pr.tag                    = 'prop_rotor_' + str(ii+1) 
            pr.origin                 = [prop_rotor_origins[ii]] 
            net.rotors.append(pr)  
    
        # Component 7 the Motors  
        prop_rotor_motor                         = RCAIDE.Components.Energy.Converters.Motor()
        prop_rotor_motor.efficiency              = 0.9
        prop_rotor_motor.nominal_voltage         = bat.pack.max_voltage*0.75 
        prop_rotor_motor.propeller_radius        = prop_rotor.tip_radius   
        prop_rotor_motor.no_load_current         = 0.1   
        prop_rotor_motor.rotor_radius            = prop_rotor.tip_radius
        prop_rotor_motor.design_torque           = prop_rotor.hover.design_torque
        prop_rotor_motor.angular_velocity        = prop_rotor.hover.design_angular_velocity/prop_rotor_motor.gear_ratio  
        prop_rotor_motor                         = size_optimal_motor(prop_rotor_motor)
        prop_rotor_motor.mass_properties.mass    = nasa_motor(prop_rotor_motor.design_torque)     
     
        # Appending motors with different origins
        for ii in range(8):
            motor        = deepcopy(prop_rotor_motor)
            motor.tag    = 'prop_rotor_motor_' + str(ii+1)
            motor.origin =  [prop_rotor_origins[ii]]                     
            net.motors.append(motor)   
                
                
        # Component 2: Payload
        payload                        = RCAIDE.Components.Energy.Peripherals.Payload()
        payload.power_draw             = 10. # Watts 
        payload.mass_properties.mass   = 1.0 * Units.kg
        net.payload                    = payload
        
        # Component 3: Avionics    
        avionics                       = RCAIDE.Components.Energy.Peripherals.Avionics()
        avionics.power_draw            = 20. # Watts  
        net.avionics                   = avionics   
        
        # Component 4: Miscellaneous Systems 
        sys                            = RCAIDE.Components.Systems.System()
        sys.mass_properties.mass       = 5 # kg       
    
        # append motor origin spanwise locations onto wing data structure
        motor_origins_front                                   = np.array(prop_rotor_origins[:4])
        motor_origins_rear                                    = np.array(prop_rotor_origins[5:])
        vehicle.wings['canard_wing'].motor_spanwise_locations = motor_origins_front[:,1]/ vehicle.wings['canard_wing'].spans.projected
        vehicle.wings['main_wing'].motor_spanwise_locations   = motor_origins_rear[:,1]/ vehicle.wings['main_wing'].spans.projected
    
        vehicle.append_component(net)  
    
     
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
    '''
    The configration set up below the scheduling of the nacelle angle and vehicle speed.
    Since one propeller operates at varying flight conditions, one must perscribe  the 
    pitch command of the propeller which us used in the variable pitch model in the analyses
    Note: low pitch at take off & low speeds, high pitch at cruise
    '''
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs                                     = RCAIDE.Components.Configs.Config.Container() 
    base_config                                 = RCAIDE.Components.Configs.Config(vehicle)
    base_config.tag                             = 'base'
    configs.append(base_config)
 
    # ------------------------------------------------------------------
    #   Hover Configuration
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Components.Configs.Config(base_config)
    config.tag                                        = 'vertical_flight'
    vector_angle                                      = 90*Units.degrees
    for rotor in config.networks.battery_electric_rotor.rotors: 
        rotor.orientation_euler_angles                 = [0,vector_angle,0] 
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle
    configs.append(config) 

    # ------------------------------------------------------------------
    #  Vertical Transition 1 
    # ------------------------------------------------------------------  
    config                                            = RCAIDE.Components.Configs.Config(base_config)
    vector_angle                                      = 60.0  * Units.degrees
    config.tag                                        = 'low_speed_transition_flight'
    for rotor in config.networks.battery_electric_rotor.rotors: 
        rotor.orientation_euler_angles                 = [0,vector_angle,0]
        rotor.inputs.pitch_command                     = rotor.cruise.design_collective_pitch*1/3
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle 
    configs.append(config)    


    # ------------------------------------------------------------------
    # Vertical Transition 2    
    # ------------------------------------------------------------------ 
    config                                            = RCAIDE.Components.Configs.Config(base_config)
    vector_angle                                      = 30.0  * Units.degrees
    config.tag                                        = 'high_speed_transition_flight'
    for rotor in config.networks.battery_electric_rotor.rotors: 
        rotor.orientation_euler_angles                 = [0,vector_angle,0]
        rotor.inputs.pitch_command                     = rotor.cruise.design_collective_pitch*2/3
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle 
    configs.append(config)    
     
     
    # ------------------------------------------------------------------
    #   Forward Flight Transition 
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Components.Configs.Config(base_config)
    config.tag                                        = 'forward_flight'
    vector_angle                                      = 0.0 * Units.degrees 
    for rotor in config.networks.battery_electric_rotor.rotors: 
        rotor.orientation_euler_angles                 = [0,vector_angle,0]
        rotor.inputs.pitch_command                     = rotor.cruise.design_collective_pitch 
    config.wings.main_wing.twists.root                = vector_angle
    config.wings.main_wing.twists.tip                 = vector_angle
    config.wings.canard_wing.twists.root              = vector_angle
    config.wings.canard_wing.twists.tip               = vector_angle 
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