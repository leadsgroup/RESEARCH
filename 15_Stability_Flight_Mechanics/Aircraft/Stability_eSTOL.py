# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units   
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor             import design_propeller 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor          import design_motor 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion import compute_motor_weight
from RCAIDE.Library.Methods.Stability.Center_of_Gravity            import compute_component_centers_of_gravity
from RCAIDE.Library.Methods.Geometry.Planform                      import segment_properties
from RCAIDE.Library.Plots                                          import *     

# python imports 
import numpy as np  
from copy import deepcopy
import matplotlib.pyplot as plt  
import os   


def vehicle_setup(): 
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Electric_STOL'    
    vehicle.mass_properties.max_takeoff               = 4090 * Units.kilogram  
    vehicle.mass_properties.takeoff                   = 4090 * Units.kilogram    
    vehicle.mass_properties.operating_empty           = 2950 * Units.kilogram  
    vehicle.mass_properties.max_zero_fuel             = 4090 * Units.kilogram # CHECK if this includes battery or not!
    vehicle.mass_properties.cargo                     = 1140  * Units.kilogram  
    vehicle.envelope.ultimate_load                    = 3.75
    vehicle.envelope.limit_load                       = 2.5 
    vehicle.reference_area                            = 38.25 * Units['meters**2']   
    vehicle.passengers                                = 9
    vehicle.systems.control                           = "fully powered" 
    vehicle.systems.accessories                       = "commuter" 
 
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################################### Landing Gear ################################################################    
    #------------------------------------------------------------------------------------------------------------------------------------
    landing_gear                    = RCAIDE.Library.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag                = "main_landing_gear" 
    landing_gear.main_tire_diameter = 0.6 * Units.m
    landing_gear.nose_tire_diameter = 0.5 * Units.m
    landing_gear.main_strut_length  = 0.75 * Units.m
    landing_gear.nose_strut_length  = 0.75 * Units.m
    landing_gear.main_units         = 2 # Number of main landing gear
    landing_gear.nose_units         = 1 # Number of nose landing gear
    landing_gear.main_wheels        = 1 # Number of wheels on the main landing gear
    landing_gear.nose_wheels        = 1 # Number of wheels on the nose landing gear      
    vehicle.landing_gear            = landing_gear

    #------------------------------------------------------------------------------------------------------------------------------------
    # ######################################################## Wings ####################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------
    #   Main Wing. Airfoil file comes from UIUC prop database, Prof. Selig
    # ------------------------------------------------------------------

    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.aspect_ratio                     = 8.47
    wing.sweeps.quarter_chord             = -3.5 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.taper                            = 1 / 3.25
    wing.spans.projected                  = 18
    wing.chords.root                      = 3.25 * Units.meter
    wing.chords.tip                       = 1 * Units.meter
    wing.chords.mean_aerodynamic          = 2.125 * Units.meter 
    wing.areas.reference                  = 38.25 * Units['meters**2']  
    wing.areas.wetted                     = 70 * Units['meters**2']  
    wing.twists.root                      = 0 * Units.degrees
    wing.twists.tip                       = 0 * Units.degrees 
    wing.origin                           = [[4.673, 0, 0.902]] 
    wing.aerodynamic_center               = [0,0,0] 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.dynamic_pressure_ratio           = 1.0

    # Wing Segments
    root_airfoil                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator
    root_airfoil.coordinate_file          = rel_path  + 'Airfoils' + separator + 'NACA_2412.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 0 * Units.deg
    segment.root_chord_percent            = 1.
    segment.thickness_to_chord            = 0.12 
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = -3.5 * Units.degrees
    segment.append_airfoil(root_airfoil)
    wing.append_segment(segment)

    tip_airfoil                           = RCAIDE.Library.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = rel_path  + 'Airfoils' + separator + 'NACA_2412.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0. * Units.degrees
    segment.root_chord_percent            = 1 / 3.25
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 0.
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12 
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)    

    # control surfaces -------------------------------------------
    flap                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.6
    flap.deflection               = 0.0 * Units.degrees
    flap.configuration_type       = 'single_slotted'
    flap.chord_fraction           = 0.50
    wing.append_control_surface(flap)

    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.6
    aileron.span_fraction_end     = 0.9
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.2
    wing.append_control_surface(aileron)

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

    wing     = RCAIDE.Library.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'

    wing.aspect_ratio            = 6.70
    wing.sweeps.quarter_chord    = 9.44 * Units.deg  
    wing.thickness_to_chord      = 0.1
    wing.taper                   = 0.5  
    wing.spans.projected         = 9.0 * Units.meter 
    wing.chords.root             = 1.8 * Units.meter 
    wing.chords.tip              = 0.9 * Units.meter 
    wing.chords.mean_aerodynamic = 1.35 * Units.meter 
    wing.areas.reference         = 12.2 * Units['meters**2']  
    wing.areas.exposed           = 25.0   # Exposed area of the horizontal tail
    wing.areas.wetted            = 25.0   # Wetted area of the horizontal tail
    wing.twists.root             = 0 * Units.degrees
    wing.twists.tip              = 0 * Units.degrees 
    wing.origin                  = [[14.918, 0, 2.6]]
    wing.aerodynamic_center      = [0,0,0] 
    wing.vertical                = False
    wing.symmetric               = True 
    wing.dynamic_pressure_ratio  = 0.9



    # Wing Segments
    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'root_segment'
    segment.percent_span_location  = 0.0
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 1.0
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 9.44  * Units.degrees 
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)

    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'tip_segment'
    segment.percent_span_location  = 1.
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 0.5               
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 0 * Units.degrees  
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)
    
    # control surfaces -------------------------------------------
    elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.05
    elevator.span_fraction_end     = 0.95
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    wing.append_control_surface(elevator)

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'

    wing.aspect_ratio            = 0.90298
    wing.sweeps.quarter_chord    = 31.2  * Units.deg   
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.4

    wing.spans.projected         = 2.27 * Units.meter 
    wing.total_length            = wing.spans.projected 
    
    wing.chords.root             = 3.75 * Units.meter 
    wing.chords.tip              = 1.5 * Units.meter  
    wing.chords.mean_aerodynamic = 2.975 * Units.meter 

    wing.areas.reference         = 5.74689 * Units['meters**2']  
    wing.areas.wetted            = 11.0  * Units['meters**2']  
    
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees

    wing.origin                  = [[11.394,0,0.328]]
    wing.aerodynamic_center      = [0,0,0]

    wing.vertical                = True
    wing.symmetric               = False
    wing.t_tail                  = True

    wing.dynamic_pressure_ratio  = 1.0


    # Wing Segments
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 1.
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = 20 * Units.degrees  
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'segment_1'
    segment.percent_span_location         = 0.195
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 1.0
    segment.dihedral_outboard             = 0. * Units.degrees
    segment.sweeps.quarter_chord          = 72.6 * Units.degrees   
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'segment_2'
    segment.percent_span_location         = 0.35
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.667 
    segment.dihedral_outboard             = 0.0 * Units.degrees
    segment.sweeps.quarter_chord          = 48.95 * Units.degrees 
    segment.thickness_to_chord            = .1  
    wing.append_segment(segment)
    
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'segment_3'
    segment.percent_span_location         = 1.00
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.4 
    segment.dihedral_outboard             = 0.0 * Units.degrees
    segment.sweeps.quarter_chord          = 0.0    
    segment.thickness_to_chord            = .1  
    wing.append_segment(segment)    
    
    # control surfaces -------------------------------------------
    rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
    rudder.tag                   = 'rudder'
    rudder.span_fraction_start   = 0.05
    rudder.span_fraction_end     = 0.95
    rudder.deflection            = 0.0  * Units.deg
    rudder.chord_fraction        = 0.4
    wing.append_control_surface(rudder)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)        

    # add to vehicle
    vehicle.append_component(wing)

    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################  Fuselage ############################################################## 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    
    fuselage                                    = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.number_coach_seats                 = vehicle.passengers 
    fuselage.seats_abreast                      = 2
    fuselage.seat_pitch                         = 1 * Units.meter 
    fuselage.fineness.nose                      = 1.6
    fuselage.fineness.tail                      = 3.4
    fuselage.lengths.nose                       = 4.0 * Units.meter
    fuselage.lengths.tail                       = 7.5 * Units.meter
    fuselage.lengths.total                      = 15.0 * Units.meter  
    fuselage.lengths.fore_space                 = 1.3 * Units.meter
    fuselage.lengths.aft_space                  = 3.5 * Units.meter
    fuselage.width                              = 2.2 * Units.meter
    fuselage.heights.maximum                    = 2.45 * Units.meter
    fuselage.effective_diameter                 = 2.2 * Units.meter
    fuselage.areas.side_projected               = 24.55 * Units['meters**2'] 
    fuselage.areas.wetted                       = 72.3 * Units['meters**2'] 
    fuselage.areas.front_projected              = 18.85 * Units['meters**2']  
    fuselage.differential_pressure              = 5e03 * Units.pascal #CHECK THIS 
    fuselage.heights.at_quarter_length          = 2.5 * Units.meter
    fuselage.heights.at_three_quarters_length   = 2.5 * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = 2.5 * Units.meter

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_0'    
    segment.percent_x_location                  = 0.0000
    segment.percent_z_location                  = -0.066
    segment.height                              = 0.0 
    segment.width                               = 0.0  
    fuselage.Segments.append(segment)   
    
    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_1'    
    segment.percent_x_location                  = 0.02561 
    segment.percent_z_location                  = -0.06567 
    segment.height                              = 0.816
    segment.width                               = 0.65
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'   
    segment.percent_x_location                  = 0.06424 
    segment.percent_z_location                  = -0.05775 
    segment.height                              = 1.200 
    segment.width                               = 1.32653 
    fuselage.Segments.append(segment)      
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'   
    segment.percent_x_location                  = 0.116 
    segment.percent_z_location                  = -0.045 
    segment.height                              = 1.7 
    segment.width                               = 1.73469 
    fuselage.Segments.append(segment)   

    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'   
    segment.percent_x_location                  = 0.1561	
    segment.percent_z_location                  = -0.037 
    segment.height                              = 2.00 
    segment.width                               = 1.9898 
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'   
    segment.percent_x_location                  = 0.19574 
    segment.percent_z_location                  = -0.029 
    segment.height                              = 2.25 
    segment.width                               = 2.14286 
    fuselage.Segments.append(segment)     
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  = 0.2531 
    segment.percent_z_location                  = -0.022 
    segment.height                              = 2.5 
    segment.width                               = 2.39796 
    fuselage.Segments.append(segment)             
     
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'   
    segment.percent_x_location                  = 0.38869 
    segment.percent_z_location                  = -0.01639 
    segment.height                              = 2.5 
    segment.width                               = 2.44898 
    fuselage.Segments.append(segment)    
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_8'   
    segment.percent_x_location                  = 0.55051
    segment.percent_z_location                  = -0.0082 
    segment.height                              = 2.12245 
    segment.width                               = 2.14286 
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_9'     
    segment.percent_x_location                  = 0.75 
    segment.percent_z_location                  = 0.0082
    segment.height                              = 1.27347
    segment.width                               = 1.17347
    fuselage.Segments.append(segment)     
        
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_10'     
    segment.percent_x_location                  = 0.8627 
    segment.percent_z_location                  = 0.01639
    segment.height                              = 0.8
    segment.width                               = 0.71429
    fuselage.Segments.append(segment)   
        
    # Segment                                   
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_11'     
    segment.percent_x_location                  = 1.0
    segment.percent_z_location                  = 0.02459
    segment.height                              = 0.0
    segment.width                               = 0.0
    fuselage.Segments.append(segment)    
    
    # add to vehicle
    vehicle.append_component(fuselage)
     

    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################################### Energy Network ##############################################################    
    #------------------------------------------------------------------------------------------------------------------------------------ 
   
    #initialize the electric network
    net                              = RCAIDE.Framework.Networks.Electric()
    
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                              = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus()


    #------------------------------------------------------------------------------------------------------------------------------------           
    # Battery
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat_module                                             = RCAIDE.Library.Components.Energy.Sources.Battery_Modules.Lithium_Ion_NMC()
    bat_module.electrical_configuration.series             = 10
    bat_module.electrical_configuration.parallel           = 210
    bat_module.cell.nominal_capacity                       = 3.8
    bat_module.geometrtic_configuration.total              = bat_module.electrical_configuration.parallel*bat_module.electrical_configuration.series  
    bat_module.voltage                                     = bat_module.maximum_voltage 
    bat_module.geometrtic_configuration.normal_count       = 42
    bat_module.geometrtic_configuration.parallel_count     = 50
    bat_module.nominal_capacity                            = bat_module.cell.nominal_capacity* bat_module.electrical_configuration.parallel

    for _ in range(14):
        bat_copy = deepcopy(bat_module)
        bus.battery_modules.append(bat_copy)

    bus.initialize_bus_properties()    

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Starboard Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    starboard_propulsor                              = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    starboard_propulsor.tag                          = 'inner_starboard_cruise_propulsor'
    starboard_propulsor.active_batteries             = ['li_ion_battery']   
  
    # Electronic Speed Controller       
    esc1                                              = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc1.tag                                          = 'esc1'
    esc1.efficiency                                   = 0.95 
    starboard_propulsor.electronic_speed_controller  = esc1   
     
    # Propeller. Optimized for cruise performance
    propeller1                                        = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    propeller1.tag                                    = 'propeller1'  
    propeller1.tip_radius                             = 1.  * Units.m
    propeller1.number_of_blades                       = 3 
    propeller1.hub_radius                             = 0.2 * Units.m
    propeller1.cruise.design_freestream_velocity      = 180.*Units['mph']   
    propeller1.cruise.design_angular_velocity         = 3500. * Units.rpm # Verify ?
    propeller1.cruise.design_Cl                       = 0.7 # CHANGE???
    propeller1.cruise.design_altitude                 = 2500. * Units.feet 
    propeller1.cruise.design_thrust                   = 800 # Based on a total 650 kw power system  
    propeller1.clockwise_rotation                     = False
    propeller1.variable_pitch                         = True  
    propeller1.origin                                 = [[3.85,2.4,0.652]]   
    airfoil                                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                                      = 'NACA_4412' 
    airfoil.coordinate_file                          =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                              =[ rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt']   
    propeller1.append_airfoil(airfoil)                       
    propeller1.airfoil_polar_stations                 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    design_propeller(propeller1)    
    starboard_propulsor.rotor                        = propeller1   
              
    # DC_Motor       
    motor1                                            = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    motor1.efficiency                                 = 0.98
    motor1.origin                                     = [[3.95,2.4,0.652]]   
    motor1.nominal_voltage                            = bus.voltage *0.5
    motor1.no_load_current                            = 1
    motor1.rotor_radius                               = propeller1.tip_radius
    motor1.design_torque                              = propeller1.cruise.design_torque
    motor1.angular_velocity                           = propeller1.cruise.design_angular_velocity 
    design_motor(motor1)  
    motor1.mass_properties.mass                       = compute_motor_weight(motor1.design_torque) 
    starboard_propulsor.motor                         = motor1
    
    nacelle1                    = RCAIDE.Library.Components.Nacelles.Stack_Nacelle()
    nacelle1.tag                = 'nacelle_1'
    nacelle1.length             = 2
    nacelle1.diameter           = 0.4 
    nacelle1.areas.wetted       = 1.8 * Units['meters**2'] 
    nacelle1.origin             = [[3.85,2.4,0.652]]
    nacelle1.flow_through       = False  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle1.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.001 
    nac_segment.percent_z_location = 0.
    nac_segment.height             = 0.3 
    nac_segment.width              = 0.3 
    nacelle1.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.4 
    nac_segment.percent_z_location = 0
    nac_segment.height             = 0.4	 
    nac_segment.width              = 0.4 
    nacelle1.append_segment(nac_segment)  
     
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.6  
    nac_segment.percent_z_location = 0
    nac_segment.height             = 0.35	 
    nac_segment.width              = 0.35 
    nacelle1.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.8  
    nac_segment.percent_z_location = 0.0
    nac_segment.height             = 0.2 
    nac_segment.width              = 0.2
    nacelle1.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 1.0   
    nac_segment.percent_z_location = 0
    nac_segment.height             = 0. 
    nac_segment.width              = 0. 
    nacelle1.append_segment(nac_segment)  
    
    starboard_propulsor.nacelle = nacelle1      
   
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Starboard takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    propulsor5                              = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    propulsor5.tag                          = 'inner_starboard_takeoff_propulsor'
    propulsor5.active_batteries             = ['li_ion_battery']   
  
    # Electronic Speed Controller       
    esc5                                              = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc5.tag                                          = 'esc_5'
    esc5.efficiency                                   = 0.95 
    propulsor5.electronic_speed_controller  = esc5 
     
    # Propeller. Optimized for low-speed flight              
    propeller5                                       = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    propeller5.tag                                   = 'propeller_5'  
    propeller5.tip_radius                            = 0.6  * Units.m
    propeller5.number_of_blades                      = 3 
    propeller5.hub_radius                            = 0.15 * Units.m
    propeller5.cruise.design_freestream_velocity     = 60.*Units['mph']   # Low speed
    propeller5.cruise.design_angular_velocity        = 4000. * Units.rpm 
    propeller5.cruise.design_Cl                      = 0.7 # CHANGE???
    propeller5.cruise.design_altitude                = 1000. * Units.feet 
    propeller5.cruise.design_thrust                  = 2400 #    
    propeller5.clockwise_rotation                    = False
    propeller5.variable_pitch                        = True  
    propeller5.origin                                = [[4.07,5.87,0.652]]   
    airfoil                                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                                      = 'NACA_4412' 
    airfoil.coordinate_file                          =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                              =[ rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt']   
    propeller5.append_airfoil(airfoil)                       
    propeller5.airfoil_polar_stations                 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    design_propeller(propeller5)    
    propulsor5.rotor                        = propeller5   
              
    # DC_Motor       
    motor5                                           = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    motor5.efficiency                                = 0.98
    motor5.origin                                    = [[4.2,5.87,0.64]]   
    motor5.nominal_voltage                           = bus.voltage*0.5
    motor5.no_load_current                           = 1
    motor5.rotor_radius                              = propeller5.tip_radius
    motor5.design_torque                             = propeller5.cruise.design_torque
    motor5.angular_velocity                          = propeller5.cruise.design_angular_velocity 
    design_motor(motor5)  
    motor5.mass_properties.mass                      = compute_motor_weight(motor5.design_torque) 
    propulsor5.motor                                 = motor5 
   
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Starboard Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    propulsor2                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    propulsor2.tag                         = "outer_starboard_cruise_propulsor"
    propulsor2.active_batteries            = ['li_ion_battery']   
            
    esc2                                      = deepcopy(esc1)
    esc2.origin                               = [[4.5, 4.2, 0.652]]      
    propulsor2.electronic_speed_controller    = esc2  

    propeller2                                = deepcopy(propeller1)
    propeller2.tag                            = 'propeller2' 
    propeller2.origin                         = [[3.85,4.25,0.652]]
    propeller2.clockwise_rotation             = False        
    propulsor2.rotor                          = propeller2  
              
    motor2                                    = deepcopy(motor1)
    motor2.origin                             = [[3.95, 4.25, 0.652]]      
    propulsor2.motor                       = motor2   

    nacelle2                                  = deepcopy(nacelle1)
    nacelle2.tag                              = 'nacelle2'
    nacelle2.origin                           = [[3.95,4.25,0.652]]
    propulsor2.nacelle                        = nacelle2
     
    # append propulsor to distribution line 
    bus.propulsors.append(propulsor2)
    '''
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Inner Port Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    propulsor3                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    propulsor3.tag                         = "inner_port_cruise_propulsor"
    propulsor3.active_batteries            = ['li_ion_battery']   
        
    esc3                                      = deepcopy(esc1)
    esc3.origin                               = [[4.50, -2.4,0.652]]    
    propulsor3.electronic_speed_controller = esc3 

    propeller3                                = deepcopy(propeller1)
    propeller3.tag                            = 'propeller3' 
    propeller3.origin                         = [[3.85, -2.4,0.652]]
    propeller3.clockwise_rotation             = False        
    propulsor3.rotor                          = propeller3  
              
    motor3                                    = deepcopy(motor1)
    motor3.origin                             = [[3.95, -2.4,0.652]]    
    propulsor3.motor                          = motor3   

    nacelle3                                  = deepcopy(nacelle1)
    nacelle3.tag                              = 'nacelle3'
    nacelle3.origin                           = [[3.95, -2.4,0.652]]
    propulsor3.nacelle                        = nacelle3
     
    # append propulsor to distribution line 
    bus.propulsors.append(propulsor3)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Port Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    propulsor4                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    propulsor4.tag                         = "outer_port_cruise_propulsor"
    propulsor4.active_batteries            = ['li_ion_battery']   
            
    esc4                                      = deepcopy(esc1)
    esc4.origin                               = [[4.5, -4.25,0.652]]     
    propulsor4.electronic_speed_controller    = esc4 

    propeller4                                = deepcopy(propeller1)
    propeller4.tag                            = 'propeller4' 
    propeller4.origin                         = [[3.85, -4.25,0.652]]
    propeller4.clockwise_rotation             = False        
    propulsor4.rotor                          = propeller4 
              
    motor4                                    = deepcopy(motor1)
    motor4.origin                             = [[3.95, -4.25,0.652]]     
    propulsor4.motor                          = motor4   

    nacelle4                                  = deepcopy(nacelle1)
    nacelle4.tag                              = 'nacelle4'
    nacelle4.origin                           = [[3.95, -4.25,0.652]]
    propulsor4.nacelle                        = nacelle4
     
    # append propulsor to distribution line 
    bus.propulsors.append(propulsor4) 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Starboard Takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    propulsor6                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    propulsor6.tag                         = "outer_starboard_takeoff_propulsor"
    propulsor6.active_batteries            = ['li_ion_battery']   
            
    esc6                                      = deepcopy(esc5)
    esc6.origin                               = [[5.0, 7.4, 0.64]]      
    propulsor6.electronic_speed_controller = esc6  

    propeller6                                = deepcopy(propeller5)
    propeller6.tag                            = 'propeller6' 
    propeller6.origin                         = [[4.07,7.4,0.64]]
    propeller6.clockwise_rotation             = False        
    propulsor6.rotor                       = propeller6  
              
    motor6                                    = deepcopy(motor5)
    motor6.origin                             = [[4.2, 7.4, 0.64]]      
    propulsor6.motor                       = motor6   

    nacelle6                                  = deepcopy(nacelle1)
    nacelle6.tag                              = 'nacelle6'
    nacelle6.origin                           = [[4.2, 7.4,0.64]]
    propulsor6.nacelle                     = nacelle6
     
    # append propulsor to distribution line 
    bus.propulsors.append(propulsor6) 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Inner Port Takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    propulsor7                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    propulsor7.tag                         = "inner_port_takeoff_propulsor"
    propulsor7.active_batteries            = ['li_ion_battery']   
            
    esc7                                      = deepcopy(esc5)
    esc7.origin                               = [[5.0, -5.87,0.64]]       
    propulsor7.electronic_speed_controller = esc7  

    propeller7                                = deepcopy(propeller5)
    propeller7.tag                            = 'propeller7' 
    propeller7.origin                         = [[4.07, -5.87,0.64]] 
    propeller7.clockwise_rotation             = False        
    propulsor7.rotor                       = propeller7 
              
    motor7                                    = deepcopy(motor5)
    motor7.origin                             = [[4.2, -5.87,0.64]]      
    propulsor7.motor                       = motor7   

    nacelle7                                  = deepcopy(nacelle1)
    nacelle7.tag                              = 'nacelle7'
    nacelle7.origin                           = [[4.2, -5.87,0.64]] 
    propulsor7.nacelle                     = nacelle7
     
    # append propulsor to distribution line 
    bus.propulsors.append(propulsor7) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Port Takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    propulsor8                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    propulsor8.tag                         = "port_propulsor"
    propulsor8.active_batteries            = ['li_ion_battery']   
            
    esc8                                      = deepcopy(esc5)
    esc8.origin                               = [[5.0, -7.4, 0.64]]      
    propulsor8.electronic_speed_controller = esc8  

    propeller8                                = deepcopy(propeller5)
    propeller8.tag                            = 'propeller8' 
    propeller8.origin                         = [[4.07, -7.4, 0.64]] 
    propeller8.clockwise_rotation             = False        
    propulsor8.rotor                       = propeller8  
              
    motor8                                    = deepcopy(motor5)
    motor8.origin                             = [[4.20, -7.4, 0.64]]      
    propulsor8.motor                       = motor8   

    nacelle8                                  = deepcopy(nacelle1)
    nacelle8.tag                              = 'nacelle8'
    nacelle8.origin                           = [[4.20, -7.4, 0.64]] 
    propulsor8.nacelle                     = nacelle8
    
    # append propulsor to distribution line 
    bus.propulsors.append(propulsor8) 
    '''
    # append bus   
    net.busses.append(bus)
    # Append energy network to aircraft 
    vehicle.append_energy_network(net)    
    
    #------------------------------------------------------------------------------------------------------------------------- 
    # Compute Center of Gravity of aircraft (Optional)
    #------------------------------------------------------------------------------------------------------------------------- 
   
    vehicle.center_of_gravity()    
    compute_component_centers_of_gravity(vehicle)
    
    #------------------------------------------------------------------------------------------------------------------------- 
    # Done ! 
    #------------------------------------------------------------------------------------------------------------------------- 
      
    return vehicle
 
# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    """This function sets up vehicle configurations for use in different parts of the mission.
    Here, this is mostly in terms of high lift settings."""
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base' 
    base_config.landing_gear.gear_condition                      = 'up'
    configs.append(base_config)

    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'cruise'
    configs.append(config)


    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'takeoff'
    config.wings['main_wing'].control_surfaces.flap.deflection  = 20. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection  = 25. * Units.deg 
    config.networks.fuel.fuel_lines['fuel_line'].propulsors['starboard_propulsor'].fan.angular_velocity =  3470. * Units.rpm
    config.networks.fuel.fuel_lines['fuel_line'].propulsors['port_propulsor'].fan.angular_velocity      =  3470. * Units.rpm
    config.landing_gear.gear_condition                          = 'up'     
    configs.append(config)

    
    # ------------------------------------------------------------------
    #   Cutback Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'cutback'
    config.wings['main_wing'].control_surfaces.flap.deflection  = 20. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection  = 20. * Units.deg
    config.networks.fuel.fuel_lines['fuel_line'].propulsors['starboard_propulsor'].fan.angular_velocity =  2780. * Units.rpm
    config.networks.fuel.fuel_lines['fuel_line'].propulsors['port_propulsor'].fan.angular_velocity      =  2780. * Units.rpm
    config.landing_gear.gear_condition                          = 'up'       
    configs.append(config)   
    
        
    
    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'landing'
    config.wings['main_wing'].control_surfaces.flap.deflection  = 30. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection  = 25. * Units.deg
    config.networks.fuel.fuel_lines['fuel_line'].propulsors['starboard_propulsor'].fan.angular_velocity =  2030. * Units.rpm
    config.networks.fuel.fuel_lines['fuel_line'].propulsors['port_propulsor'].fan.angular_velocity      =  2030. * Units.rpm
    config.landing_gear.gear_condition                          = 'down'    
    configs.append(config)   
     
    # ------------------------------------------------------------------
    #   Short Field Takeoff Configuration
    # ------------------------------------------------------------------  

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'reverse_thrust'
    config.wings['main_wing'].control_surfaces.flap.deflection  = 30. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection  = 25. * Units.deg 
    config.landing_gear.gear_condition                          = 'down'    
    configs.append(config)    


    return configs

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------


# Step 1 design a vehicle
vehicle  = vehicle_setup()    

# plot vehicle 
plot_3d_vehicle(vehicle,
                min_x_axis_limit            = -5,
                max_x_axis_limit            = 40,
                min_y_axis_limit            = -20,
                max_y_axis_limit            = 20,
                min_z_axis_limit            = -20,
                max_z_axis_limit            = 20)          
