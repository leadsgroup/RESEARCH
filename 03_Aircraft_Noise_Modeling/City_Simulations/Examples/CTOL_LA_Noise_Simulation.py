# digital_elevation_and_noise_hemispheres_test.py
#
# Created: Dec 2023 M. Clarke  

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE Imports 
import RCAIDE
from RCAIDE.Framework.Core import Units , Data 
from RCAIDE.Library.Plots import *     
from RCAIDE.Library.Methods.Noise.Metrics import *  
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor             import design_propeller 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor          import design_motor 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion import nasa_motor
from RCAIDE.Library.Methods.Energy.Sources.Batteries.Common         import initialize_from_circuit_configuration
from RCAIDE.Library.Methods.Geometry.Planform                       import wing_segmented_planform 
from RCAIDE.Library.Methods.Noise.Common.generate_microphone_locations        import generate_terrain_elevated_microphone_locations
from RCAIDE.Library.Mission.Common.compute_point_to_point_geospacial_data     import compute_point_to_point_geospacial_data

# Python imports
import matplotlib.pyplot as plt  
import sys 
import numpy as np     
from copy import deepcopy
import os
  

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():    
    microphone_terrain_data =  generate_terrain_elevated_microphone_locations(topography_file   ='LA_Metropolitan_Area.txt',
                                                           ground_microphone_x_resolution    = 300,  
                                                           ground_microphone_y_resolution    = 90, 
                                                           ground_microphone_x_stencil       = 1,   
                                                           ground_microphone_y_stencil       = 1)    
    

    geospacial_data =  compute_point_to_point_geospacial_data(topography_file  = 'LA_Metropolitan_Area.txt',
                                                                                 departure_tag                         = 'A',
                                                                                 destination_tag                       = 'B',
                                                                                 departure_coordinates                 = [33.94067953101678, -118.40513722978149],
                                                                                 destination_coordinates               = [33.81713622114423, -117.92111163722772] )    
  
    
    plot_elevation_contours(topography_file   ='LA_Metropolitan_Area.txt',use_lat_long_coordinates = True, save_filename = "Elevation_Contours_Lat_Long")

    plot_elevation_contours(topography_file   ='LA_Metropolitan_Area.txt',use_lat_long_coordinates = False, save_filename = "Elevation_Contours_XY")  
      
    vehicle  = vehicle_setup()      
    vehicle.networks.electric.busses.bus.identical_propulsors     = False # only for regression     
    configs  = configs_setup(vehicle) 
    analyses = analyses_setup(configs,microphone_terrain_data,geospacial_data)  
    mission  = mission_setup(analyses,geospacial_data)
    missions = missions_setup(mission)  
    results  = missions.base_mission.evaluate()   
    
    plot_results(results)   
     
    return      

# ----------------------------------------------------------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------------------------------------------------------
def vehicle_setup():

    #------------------------------------------------------------------------------------------------------------------------------------
    #   Initialize the Vehicle
    #------------------------------------------------------------------------------------------------------------------------------------

    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'X57_Maxwell_Mod2'

 
    # ################################################# Vehicle-level Properties ########################################################  

    # mass properties
    vehicle.mass_properties.max_takeoff   = 2550. * Units.pounds
    vehicle.mass_properties.takeoff       = 2550. * Units.pounds
    vehicle.mass_properties.max_zero_fuel = 2550. * Units.pounds 
    vehicle.flight_envelope.ultimate_load        = 5.7
    vehicle.flight_envelope.positive_limit_load           = 3.8 
    vehicle.reference_area                = 14.76
    vehicle.passengers                    = 4
    vehicle.systems.control               = "fully powered"
    vehicle.systems.accessories           = "commuter"    
    
    cruise_speed                          = 135.*Units['mph']    
    altitude                              = 2500. * Units.ft
    atmo                                  = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    freestream                            = atmo.compute_values (0.)
    freestream0                           = atmo.compute_values (altitude)
    mach_number                           = (cruise_speed/freestream.speed_of_sound)[0][0] 
    vehicle.design_dynamic_pressure       = ( .5 *freestream0.density*(cruise_speed*cruise_speed))[0][0]
    vehicle.design_mach_number            =  mach_number

         
    # ##########################################################  Wings ################################################################    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Main Wing
    #------------------------------------------------------------------------------------------------------------------------------------
    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.sweeps.quarter_chord             = 0.0 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 14.76
    wing.spans.projected                  = 11.4 
    wing.chords.root                      = 1.46
    wing.chords.tip                       = 0.92
    wing.chords.mean_aerodynamic          = 1.19
    wing.taper                            = wing.chords.root/wing.chords.tip 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 3.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[2.93, 0., 1.01]]
    wing.aerodynamic_center               = [3., 0., 1.01] 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0  
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator + '..'  + separator + '..' + separator + '..' + separator  
    airfoil                               = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                           = 'NACA_63_412.txt' 
    airfoil.coordinate_file               = rel_path + 'Airfoils' + separator + 'NACA_63_412.txt'   # absolute path     
    cg_x                                  = wing.origin[0][0] + 0.25*wing.chords.mean_aerodynamic
    cg_z                                  = wing.origin[0][2] - 0.2*wing.chords.mean_aerodynamic
    vehicle.mass_properties.center_of_gravity = [[cg_x,   0.  ,  cg_z ]]  # SOURCE: Design and aerodynamic analysis of a twin-engine commuter aircraft

    # Wing Segments
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'inboard'
    segment.percent_span_location         = 0.0 
    segment.twist                         = 3. * Units.degrees   
    segment.root_chord_percent            = 1. 
    segment.dihedral_outboard             = 0.  
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)

    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'outboard'
    segment.percent_span_location         = 0.5438
    segment.twist                         = 2.* Units.degrees 
    segment.root_chord_percent            = 1. 
    segment.dihedral_outboard             = 0. 
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12 
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)
    
    # Wing Segments
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'winglet'
    segment.percent_span_location         = 0.98
    segment.twist                         = 1.  * Units.degrees 
    segment.root_chord_percent            = 0.630
    segment.dihedral_outboard             = 75. * Units.degrees 
    segment.sweeps.quarter_chord          = 30. * Units.degrees 
    segment.thickness_to_chord            = 0.12 
    segment.append_airfoil(airfoil)
    wing.append_segment(segment) 

    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0. * Units.degrees 
    segment.root_chord_percent            = 0.12
    segment.dihedral_outboard             = 0.
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)    
    
    # Fill out more segment properties automatically
    wing = wing_segmented_planform(wing)           
    
    # add to vehicle
    vehicle.append_component(wing)


    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Horizontal Tail
    #------------------------------------------------------------------------------------------------------------------------------------    
    wing                                  = RCAIDE.Library.Components.Wings.Wing()
    wing.tag                              = 'horizontal_stabilizer' 
    wing.sweeps.quarter_chord             = 0.0 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 2.540 
    wing.spans.projected                  = 3.3  * Units.meter 
    wing.sweeps.quarter_chord             = 0 * Units.deg 
    wing.chords.root                      = 0.769 * Units.meter 
    wing.chords.tip                       = 0.769 * Units.meter 
    wing.chords.mean_aerodynamic          = 0.769 * Units.meter  
    wing.taper                            = 1. 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 0.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[7.7, 0., 0.25]]
    wing.aerodynamic_center               = [7.8, 0., 0.25] 
    wing.vertical                         = False
    wing.winglet_fraction                 = 0.0  
    wing.symmetric                        = True
    wing.high_lift                        = False 
    wing.dynamic_pressure_ratio           = 0.9

    # add to vehicle
    vehicle.append_component(wing)


    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Vertical Stabilizer
    #------------------------------------------------------------------------------------------------------------------------------------ 
    wing                                  = RCAIDE.Library.Components.Wings.Wing()
    wing.tag                              = 'vertical_stabilizer'     
    wing.sweeps.quarter_chord             = 25. * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 2.258 * Units['meters**2']  
    wing.spans.projected                  = 1.854   * Units.meter  
    wing.chords.root                      = 1.6764 * Units.meter 
    wing.chords.tip                       = 0.6858 * Units.meter 
    wing.chords.mean_aerodynamic          = 1.21   * Units.meter 
    wing.taper                            = wing.chords.tip/wing.chords.root 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 0.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[6.75 ,0, 0.623]]
    wing.aerodynamic_center               = [0.508 ,0,0]  
    wing.vertical                         = True 
    wing.symmetric                        = False
    wing.t_tail                           = False
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0

    # add to vehicle
    vehicle.append_component(wing)

 
    # ##########################################################   Fuselage  ############################################################    
    fuselage = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.seats_abreast                      = 2.
    fuselage.fineness.nose                      = 1.6
    fuselage.fineness.tail                      = 2.
    fuselage.lengths.nose                       = 60.  * Units.inches
    fuselage.lengths.tail                       = 161. * Units.inches
    fuselage.lengths.cabin                      = 105. * Units.inches
    fuselage.lengths.total                      = 332.2* Units.inches
    fuselage.lengths.fore_space                 = 0.
    fuselage.lengths.aft_space                  = 0.
    fuselage.width                              = 42. * Units.inches
    fuselage.heights.maximum                    = 62. * Units.inches
    fuselage.heights.at_quarter_length          = 62. * Units.inches
    fuselage.heights.at_three_quarters_length   = 62. * Units.inches
    fuselage.heights.at_wing_root_quarter_chord = 23. * Units.inches
    fuselage.areas.side_projected               = 8000.  * Units.inches**2.
    fuselage.areas.wetted                       = 30000. * Units.inches**2.
    fuselage.areas.front_projected              = 42.* 62. * Units.inches**2.
    fuselage.effective_diameter                 = 50. * Units.inches 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0
    segment.percent_z_location                  = 0
    segment.height                              = 0.01
    segment.width                               = 0.01
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  = 0.007279116466
    segment.percent_z_location                  = 0.002502014453
    segment.height                              = 0.1669064748
    segment.width                               = 0.2780205877
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.01941097724
    segment.percent_z_location                  = 0.001216095397
    segment.height                              = 0.3129496403
    segment.width                               = 0.4365777215
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.06308567604
    segment.percent_z_location                  = 0.007395489231
    segment.height                              = 0.5841726619
    segment.width                               = 0.6735119903
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.1653761217
    segment.percent_z_location                  = 0.02891281352
    segment.height                              = 1.064028777
    segment.width                               = 1.067200529
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.2426372155
    segment.percent_z_location                  = 0.04214148761
    segment.height                              = 1.293766653
    segment.width                               = 1.183058255
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.2960174029
    segment.percent_z_location                  = 0.04705241831
    segment.height                              = 1.377026712
    segment.width                               = 1.181540054
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 0.3809404284
    segment.percent_z_location                  = 0.05313580461
    segment.height                              = 1.439568345
    segment.width                               = 1.178218989
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_8'
    segment.percent_x_location                  = 0.5046854083
    segment.percent_z_location                  = 0.04655492473
    segment.height                              = 1.29352518
    segment.width                               = 1.054390707
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_9'
    segment.percent_x_location                  = 0.6454149933
    segment.percent_z_location                  = 0.03741966266
    segment.height                              = 0.8971223022
    segment.width                               = 0.8501926505
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_10'
    segment.percent_x_location                  = 0.985107095
    segment.percent_z_location                  = 0.04540283436
    segment.height                              = 0.2920863309
    segment.width                               = 0.2012565415
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_11'
    segment.percent_x_location                  = 1
    segment.percent_z_location                  = 0.04787575562
    segment.height                              = 0.1251798561
    segment.width                               = 0.1206021048
    fuselage.Segments.append(segment)

    # add to vehicle
    vehicle.append_component(fuselage)
 
    # ########################################################  Energy Network  #########################################################  
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
    bat.pack.electrical_configuration.series               = 140   
    bat.pack.electrical_configuration.parallel             = 100
    initialize_from_circuit_configuration(bat)  
    bat.module.number_of_modules                           = 14  
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module.voltage                                     = bat.pack.maximum_voltage/bat.module.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
    bat.module.geometrtic_configuration.normal_count       = 24
    bat.module.geometrtic_configuration.parallel_count     = 40     
    bus.voltage                                            = bat.pack.maximum_voltage  
    bus.batteries.append(bat)            
    

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    starboard_propulsor                              = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    starboard_propulsor.tag                          = 'starboard_propulsor'
    starboard_propulsor.active_batteries             = ['li_ion_battery']   
  
    # Electronic Speed Controller       
    esc                                              = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc.tag                                          = 'esc_1'
    esc.efficiency                                   = 0.95 
    starboard_propulsor.electronic_speed_controller  = esc   
     
    # Propeller              
    propeller                                        = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    propeller.tag                                    = 'propeller_1'  
    propeller.tip_radius                             = 1.72/2   
    propeller.number_of_blades                       = 3
    propeller.hub_radius                             = 10.     * Units.inches 
    propeller.cruise.design_freestream_velocity      = 175.*Units['mph']   
    propeller.cruise.design_angular_velocity         = 2700. * Units.rpm 
    propeller.cruise.design_Cl                       = 0.7 
    propeller.cruise.design_altitude                 = 2500. * Units.feet 
    propeller.cruise.design_thrust                   = 2000   
    propeller.clockwise_rotation                     = False
    propeller.variable_pitch                         = True  
    propeller.origin                                 = [[2.,2.5,0.95]]   
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
    motor.origin                                     = [[2.,  2.5, 0.95]]
    motor.nominal_voltage                            = bat.pack.maximum_voltage*0.5
    motor.no_load_current                            = 1
    motor.rotor_radius                               = propeller.tip_radius
    motor.design_torque                              = propeller.cruise.design_torque
    motor.angular_velocity                           = propeller.cruise.design_angular_velocity 
    design_motor(motor)  
    motor.mass_properties.mass                       = nasa_motor(motor.design_torque) 
    starboard_propulsor.motor                        = motor 
 

 
    # ##########################################################   Nacelles  ############################################################    
    nacelle                        = RCAIDE.Library.Components.Nacelles.Stack_Nacelle()
    nacelle.tag                    = 'nacelle_1'
    nacelle.length                 = 2
    nacelle.diameter               = 42 * Units.inches
    nacelle.areas.wetted           = 0.01*(2*np.pi*0.01/2)
    nacelle.origin                 = [[2.5,2.5,1.0]]
    nacelle.flow_through           = False  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.1  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.3  
    nac_segment.height             = 0.52
    nac_segment.width              = 0.7
    nacelle.append_segment(nac_segment)  
     
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.5  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.7 
    nac_segment.height             = 0.4
    nac_segment.width              = 0.6
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.9 
    nac_segment.height             = 0.3
    nac_segment.width              = 0.5
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 1.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)    
    
    starboard_propulsor.nacelle = nacelle
    
    # append propulsor to distribution line 
    bus.propulsors.append(starboard_propulsor) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Port Propulsor
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
    # Payload 
    #------------------------------------------------------------------------------------------------------------------------------------  
    payload                      = RCAIDE.Library.Components.Payloads.Payload()
    payload.power_draw           = 10. # Watts
    payload.mass_properties.mass = 1.0 * Units.kg
    bus.payload                  = payload

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Avionics
    #------------------------------------------------------------------------------------------------------------------------------------  
    avionics                     = RCAIDE.Library.Components.Systems.Avionics()
    avionics.power_draw          = 20. # Watts
    bus.avionics                 = avionics   

    # append bus   
    net.busses.append(bus)
    
    vehicle.append_energy_network(net)

    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------
    
    return vehicle

# ---------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):

    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------  
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'  
    configs.append(base_config)   
    
    # done!
    return configs


# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ---------------------------------------------------------------------- 
def analyses_setup(configs,microphone_terrain_data,geospacial_data):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis      = base_analysis(config,microphone_terrain_data,geospacial_data) 
        analyses[tag] = analysis

    return analyses  


def base_analysis(vehicle,microphone_terrain_data,geospacial_data):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle() 
 
    # ------------------------------------------------------------------
    #  Weights
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics          = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)   
 
    #  Noise Analysis   
    noise = RCAIDE.Framework.Analyses.Noise.Frequency_Domain_Buildup()   
    noise.geometry = vehicle
    noise.settings.mean_sea_level_altitude          = False  
    noise.settings.aircraft_departure_location      = geospacial_data.departure_location   
    noise.settings.aircraft_destination_location    = geospacial_data.destination_location       
    noise.settings.aircraft_departure_coordinates   = geospacial_data.departure_coordinates
    noise.settings.aircraft_destination_coordinates = geospacial_data.destination_coordinates
    noise.settings.ground_microphone_x_resolution   = microphone_terrain_data.ground_microphone_x_resolution           
    noise.settings.ground_microphone_y_resolution   = microphone_terrain_data.ground_microphone_y_resolution          
    noise.settings.ground_microphone_x_stencil      = microphone_terrain_data.ground_microphone_x_stencil             
    noise.settings.ground_microphone_y_stencil      = microphone_terrain_data.ground_microphone_y_stencil             
    noise.settings.ground_microphone_min_y          = microphone_terrain_data.ground_microphone_min_x                 
    noise.settings.ground_microphone_max_y          = microphone_terrain_data.ground_microphone_max_x                 
    noise.settings.ground_microphone_min_x          = microphone_terrain_data.ground_microphone_min_y                 
    noise.settings.ground_microphone_max_x          = microphone_terrain_data.ground_microphone_max_y            
    noise.settings.topography_file                  = microphone_terrain_data.topography_file  
    noise.settings.ground_microphone_locations      = microphone_terrain_data.ground_microphone_locations 
    noise.settings.ground_microphone_min_lat        = microphone_terrain_data.ground_microphone_min_lat 
    noise.settings.ground_microphone_max_lat        = microphone_terrain_data.ground_microphone_max_lat 
    noise.settings.ground_microphone_min_long       = microphone_terrain_data.ground_microphone_min_long
    noise.settings.ground_microphone_min_long       = microphone_terrain_data.ground_microphone_min_long        
    noise.settings.ground_microphone_coordinates    = microphone_terrain_data.ground_microphone_coordinates   
    analyses.append(noise)

    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle  = vehicle 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Earth()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses    

# ----------------------------------------------------------------------
#  Set Up Mission 
# ---------------------------------------------------------------------- 
def mission_setup(analyses,geospacial_data):      
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission       = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag   = 'mission' 
    Segments      = RCAIDE.Framework.Mission.Segments  
    base_segment  = Segments.Segment()   
    base_segment.state.numerics.number_of_control_points  = 16  
    base_segment.state.numerics.discretization_method  = RCAIDE.Library.Methods.Utilities.Chebyshev.linear_data 
    
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 

    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = "climb"   
    segment.analyses.extend( analyses.base ) 
    segment.initial_battery_state_of_charge              = 0.89         
    segment.altitude_start                               = 10.0    * Units.feet  
    segment.altitude_end                                 = 100.0   * Units.feet 
    segment.air_speed_start                              = 100.    * Units['mph'] 
    segment.air_speed_end                                = 120.    * Units['mph'] 
    segment.climb_rate                                   = 50.     * Units['ft/min']         
    segment.true_course                                  = geospacial_data.true_course 
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                      = True  
    segment.flight_dynamics.force_z                      = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                
       
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Constant Altitude Cruises 
    # ------------------------------------------------------------------   
    segment                                               = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                           = "Cruise" 
    segment.analyses.extend( analyses.base)                 
    segment.initial_battery_state_of_charge               = 0.89         
    segment.altitude                                      = 100. * Units.ft 
    segment.air_speed                                     = 120.    * Units['mph'] 
    segment.distance                                      = 1000    
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                 
           
    mission.append_segment(segment)  
        
     
    return mission

# ----------------------------------------------------------------------
#  Set Up Missions 
# ---------------------------------------------------------------------- 
def missions_setup(mission): 
 
    missions     = RCAIDE.Framework.Mission.Missions() 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions  


# ----------------------------------------------------------------------
#  Plot Resuls 
# ---------------------------------------------------------------------- 
def plot_results(results): 
    
    noise_data   = post_process_noise_data(results)  
    
    # Plot noise hemisphere
    plot_noise_hemisphere(noise_data,
                          noise_level      = noise_data.SPL_dBA[1], 
                          min_noise_level  = 35,  
                          max_noise_level  = 90, 
                          noise_scale_label= 'SPL [dBA]')     
    

    # Plot noise hemisphere with vehicle 
    plot_noise_hemisphere(noise_data,
                          noise_level      = noise_data.SPL_dBA[1], 
                          min_noise_level  = 35,  
                          max_noise_level  = 90, 
                          noise_scale_label= 'SPL [dBA]',
                          save_filename    = "Noise_Hemisphere_With_Aircraft", 
                          vehicle          = results.segments.climb.analyses.aerodynamics.vehicle)      
    
    
    # Plot noise level
    flight_times = np.array(['06:00:00','07:00:00','08:00:00','09:00:00','10:00:00','11:00:00','12:00:00','13:00:00','14:00:00','15:00:00'])  
      
    noise_data      = post_process_noise_data(results)   
    noise_data      = DNL_noise_metric(noise_data, flight_times,time_period = 24*Units.hours)
    noise_data      = Equivalent_noise_metric(noise_data, flight_times,time_period = 15*Units.hours)
    noise_data      = SENEL_noise_metric(noise_data, flight_times,time_period = 24*Units.hours)
    
    plot_noise_level(noise_data,
                    noise_level  = noise_data.SPL_dBA[0], 
                    save_filename="Sideline_Noise_Levels")  
    
    # Maximum Sound Pressure Level   
    plot_3D_noise_contour(noise_data,
                          noise_level      = np.max(noise_data.SPL_dBA,axis=0), 
                          min_noise_level  = 35,  
                          max_noise_level  = 90, 
                          noise_scale_label= 'SPL [dBA]',
                          save_filename    = "SPL_max_Noise_3D_Contour")   
                        

    # Day Night Average Noise Level 
    plot_3D_noise_contour(noise_data,
                        noise_level      = noise_data.DNL,
                        min_noise_level  = 35,  
                        max_noise_level  = 90, 
                        noise_scale_label= 'DNL',
                        show_microphones = True, 
                        save_filename    = "DNL_Noise_3D_Contour") 
    

    # Equivalent Noise Level
    plot_3D_noise_contour(noise_data,
                        noise_level      = noise_data.L_AeqT,
                        min_noise_level  = 35,  
                        max_noise_level  = 90, 
                        noise_scale_label= 'LAeqT',
                        show_trajectory  = True,
                        save_filename    = "LAeqT_Noise_3D_Contour")    
    

    # 24-hr Equivalent Noise Level
    plot_3D_noise_contour(noise_data,
                       noise_level      = noise_data.L_AeqT,
                       min_noise_level  = 35,  
                       max_noise_level  = 90, 
                       noise_scale_label= '24hr-LAeqT',
                       save_filename    = "24hr_LAeqT_Noise_3D_Contour", 
                       use_lat_long_coordinates = False)      
    

    # Single Event Noise Exposure Level
    plot_3D_noise_contour(noise_data,
                       noise_level      = noise_data.SENEL,
                       min_noise_level  = 35,  
                       max_noise_level  = 90, 
                       noise_scale_label= 'SENEL',
                       save_filename    = "SENEL_Noise_3D_Contour")  


    noise_data      = post_process_noise_data(results)  
    noise_data      = Equivalent_noise_metric(noise_data, flight_times,time_period = 15*Units.hours)
    noise_data      = SENEL_noise_metric(noise_data, flight_times,time_period = 24*Units.hours)
    noise_data      = DNL_noise_metric(noise_data, flight_times,time_period = 24*Units.hours)
    
    # Maximum Sound Pressure Level   
    plot_2D_noise_contour(noise_data,
                        noise_level      = np.max(noise_data.SPL_dBA,axis=0), 
                        min_noise_level  = 35,  
                        max_noise_level  = 90, 
                        noise_scale_label= 'SPL [dBA]',
                        save_filename    = "SPL_max_Noise_2D_Contour",
                        show_elevation   = True,
                        use_lat_long_coordinates= False)   
                        

    # Day Night Average Noise Level 
    plot_2D_noise_contour(noise_data,
                        noise_level      = noise_data.DNL,
                        min_noise_level  = 35,  
                        max_noise_level  = 90, 
                        noise_scale_label= 'DNL',
                        save_filename    = "DNL_Noise_2D_Contour") 
    

    # Equivalent Noise Level
    plot_2D_noise_contour(noise_data,
                        noise_level      = noise_data.L_AeqT,
                        min_noise_level  = 35,  
                        max_noise_level  = 90, 
                        noise_scale_label= 'LAeqT',
                        save_filename    = "LAeqT_Noise_2D_Contour")    
    

    # 24-hr Equivalent Noise Level
    plot_2D_noise_contour(noise_data,
                       noise_level      = noise_data.L_AeqT,
                       min_noise_level  = 35,  
                       max_noise_level  = 90, 
                       noise_scale_label= '24hr-LAeqT',
                       save_filename    = "24hr_LAeqT_Noise_2D_Contour")      
    

    # Single Event Noise Exposure Level
    plot_2D_noise_contour(noise_data,
                       noise_level      = noise_data.SENEL,
                       min_noise_level  = 35,  
                       max_noise_level  = 90, 
                       noise_scale_label= 'SENEL',
                       save_filename    = "SENEL_Noise_2D_Contour")      
    return  

if __name__ == '__main__': 
    main()    
    plt.show()
