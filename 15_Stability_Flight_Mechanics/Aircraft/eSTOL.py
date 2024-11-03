# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units   
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor             import design_propeller 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor          import design_motor 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion import nasa_motor
from RCAIDE.Library.Methods.Energy.Sources.Batteries.Common         import initialize_from_circuit_configuration
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

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                              = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus() 

    #------------------------------------------------------------------------------------------------------------------------------------           
    # Battery
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat                                                    = RCAIDE.Library.Components.Energy.Sources.Batteries.Lithium_Ion_NMC()
    bat.tag                                                = 'li_ion_battery'
    bat.pack.electrical_configuration.series               = 140  # CHANGE  
    bat.pack.electrical_configuration.parallel             = 100  # CHANGE
    initialize_from_circuit_configuration(bat)  
    bat.module.number_of_modules                           = 14   # CHANGE
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module.voltage                                     = bat.pack.maximum_voltage/bat.module.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
    bat.module.geometrtic_configuration.normal_count       = 24 # CHANGE
    bat.module.geometrtic_configuration.parallel_count     = 40 # CHANGE
    bat.thermal_management_system.heat_acquisition_system  = RCAIDE.Library.Components.Thermal_Management.Batteries.Heat_Acquisition_Systems.Direct_Air()    # DO we need this??  
    bus.voltage                                            = bat.pack.maximum_voltage  
    bus.batteries.append(bat)            
    

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Starboard Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    starboard_propulsor                              = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    starboard_propulsor.tag                          = 'inner_starboard_cruise_propulsor'
    starboard_propulsor.active_batteries             = ['li_ion_battery']   
  
    # Electronic Speed Controller       
    esc                                              = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc.tag                                          = 'esc_1'
    esc.efficiency                                   = 0.95 
    starboard_propulsor.electronic_speed_controller  = esc   
     
    # Propeller              
    propeller                                        = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    propeller.tag                                    = 'propeller_1'  
    propeller.tip_radius                             = 1.2  * Units.m
    propeller.number_of_blades                       = 6 
    propeller.hub_radius                             = 0.2 * Units.m
    propeller.cruise.design_freestream_velocity      = 175.*Units['mph']   
    propeller.cruise.design_angular_velocity         = 2700. * Units.rpm # CHANGE
    propeller.cruise.design_Cl                       = 0.7 # CHANGE???
    propeller.cruise.design_altitude                 = 2500. * Units.feet 
    propeller.cruise.design_thrust                   = # CHANGE???   
    propeller.clockwise_rotation                     = False
    propeller.variable_pitch                         = True  
    propeller.origin                                 = [[4.068,2.4,0.64]]   
    airfoil                                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                                      = 'NACA_4412' 
    airfoil.coordinate_file                          =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                              =[ rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt']   
    propeller.append_airfoil(airfoil)                       
    propeller.airfoil_polar_stations                 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    design_propeller(propeller)    
    starboard_propulsor.rotor                        = propeller   
              
    # DC_Motor       
    motor                                            = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    motor.efficiency                                 = 0.98
    motor.origin                                     = [[4.0,2.4,0.64]]   
    motor.nominal_voltage                            = bat.pack.maximum_voltage*0.5
    motor.no_load_current                            = 1
    motor.rotor_radius                               = propeller.tip_radius
    motor.design_torque                              = propeller.cruise.design_torque
    motor.angular_velocity                           = propeller.cruise.design_angular_velocity 
    design_motor(motor)  
    motor.mass_properties.mass                       = nasa_motor(motor.design_torque) 
    starboard_propulsor.motor                        = motor 
   
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Inner Starboard takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    starboard_propulsor                              = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    starboard_propulsor.tag                          = 'inner_starboard_takeoff_propulsor'
    starboard_propulsor.active_batteries             = ['li_ion_battery']   
  
    # Electronic Speed Controller       
    esc                                              = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc.tag                                          = 'esc_5'
    esc.efficiency                                   = 0.95 
    starboard_propulsor.electronic_speed_controller  = esc   
     
    # Propeller              
    propeller                                        = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    propeller.tag                                    = 'propeller_5'  
    propeller.tip_radius                             = 0.8  * Units.m # CHECK THIS
    propeller.number_of_blades                       = 6 
    propeller.hub_radius                             = 0.15 * Units.m
    propeller.cruise.design_freestream_velocity      = 60.*Units['mph']   # CHANGE?? Cruise speed
    propeller.cruise.design_angular_velocity         = 2700. * Units.rpm # CHANGE
    propeller.cruise.design_Cl                       = 0.7 # CHANGE???
    propeller.cruise.design_altitude                 = 1000. * Units.feet 
    propeller.cruise.design_thrust                   = # CHANGE???   
    propeller.clockwise_rotation                     = False
    propeller.variable_pitch                         = True  
    propeller.origin                                 = [[4.068,5.87,0.64]]   
    airfoil                                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                                      = 'NACA_4412' 
    airfoil.coordinate_file                          =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                              =[ rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt']   
    propeller.append_airfoil(airfoil)                       
    propeller.airfoil_polar_stations                 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    design_propeller(propeller)    
    starboard_propulsor.rotor                        = propeller   
              
    # DC_Motor       
    motor                                            = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    motor.efficiency                                 = 0.98
    motor.origin                                     = [[4.0,5.87,0.64]]   
    motor.nominal_voltage                            = bat.pack.maximum_voltage*0.5
    motor.no_load_current                            = 1
    motor.rotor_radius                               = propeller.tip_radius
    motor.design_torque                              = propeller.cruise.design_torque
    motor.angular_velocity                           = propeller.cruise.design_angular_velocity 
    design_motor(motor)  
    motor.mass_properties.mass                       = nasa_motor(motor.design_torque) 
    starboard_propulsor.motor                        = motor 
   
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Starboard Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_batteries            = ['li_ion_battery']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.origin                               = [[2., -2.5, 0.95]]      
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         = [[2.,-2.5,0.95]]
    propeller_2.clockwise_rotation             = False        
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.origin                             = [[2., -2.5, 0.95]]      
    port_propulsor.motor                       = motor_2   

    nacelle_2                                  = deepcopy(nacelle)
    nacelle_2.tag                              = 'nacelle_2'
    nacelle_2.origin                           = [[2.5,-2.5,1.0]]
    port_propulsor.nacelle                     = nacelle_2
     
    # append propulsor to distribution line 
    bus.propulsors.append(port_propulsor)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Inner Port Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_batteries            = ['li_ion_battery']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.origin                               = [[2., -2.5, 0.95]]      
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         = [[2.,-2.5,0.95]]
    propeller_2.clockwise_rotation             = False        
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.origin                             = [[2., -2.5, 0.95]]      
    port_propulsor.motor                       = motor_2   

    nacelle_2                                  = deepcopy(nacelle)
    nacelle_2.tag                              = 'nacelle_2'
    nacelle_2.origin                           = [[2.5,-2.5,1.0]]
    port_propulsor.nacelle                     = nacelle_2
     
    # append propulsor to distribution line 
    bus.propulsors.append(port_propulsor)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Port Cruise Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_batteries            = ['li_ion_battery']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.origin                               = [[2., -2.5, 0.95]]      
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         = [[2.,-2.5,0.95]]
    propeller_2.clockwise_rotation             = False        
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.origin                             = [[2., -2.5, 0.95]]      
    port_propulsor.motor                       = motor_2   

    nacelle_2                                  = deepcopy(nacelle)
    nacelle_2.tag                              = 'nacelle_2'
    nacelle_2.origin                           = [[2.5,-2.5,1.0]]
    port_propulsor.nacelle                     = nacelle_2
     
    # append propulsor to distribution line 
    bus.propulsors.append(port_propulsor) 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Starboard Takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_batteries            = ['li_ion_battery']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.origin                               = [[2., -2.5, 0.95]]      
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         = [[2.,-2.5,0.95]]
    propeller_2.clockwise_rotation             = False        
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.origin                             = [[2., -2.5, 0.95]]      
    port_propulsor.motor                       = motor_2   

    nacelle_2                                  = deepcopy(nacelle)
    nacelle_2.tag                              = 'nacelle_2'
    nacelle_2.origin                           = [[2.5,-2.5,1.0]]
    port_propulsor.nacelle                     = nacelle_2
     
    # append propulsor to distribution line 
    bus.propulsors.append(port_propulsor) 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Inner Port Takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_batteries            = ['li_ion_battery']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.origin                               = [[2., -2.5, 0.95]]      
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         = [[2.,-2.5,0.95]]
    propeller_2.clockwise_rotation             = False        
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.origin                             = [[2., -2.5, 0.95]]      
    port_propulsor.motor                       = motor_2   

    nacelle_2                                  = deepcopy(nacelle)
    nacelle_2.tag                              = 'nacelle_2'
    nacelle_2.origin                           = [[2.5,-2.5,1.0]]
    port_propulsor.nacelle                     = nacelle_2
     
    # append propulsor to distribution line 
    bus.propulsors.append(port_propulsor) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Outer Port Takeoff Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_batteries            = ['li_ion_battery']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.origin                               = [[2., -2.5, 0.95]]      
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         = [[2.,-2.5,0.95]]
    propeller_2.clockwise_rotation             = False        
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.origin                             = [[2., -2.5, 0.95]]      
    port_propulsor.motor                       = motor_2   

    nacelle_2                                  = deepcopy(nacelle)
    nacelle_2.tag                              = 'nacelle_2'
    nacelle_2.origin                           = [[2.5,-2.5,1.0]]
    port_propulsor.nacelle                     = nacelle_2
     
    # append propulsor to distribution line 
    bus.propulsors.append(port_propulsor) 
    
    
    
   
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
