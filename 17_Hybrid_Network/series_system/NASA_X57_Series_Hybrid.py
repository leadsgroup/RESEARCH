# Regression/scripts/Vehicles/NASA_X57.py
# 
# 
# Created:  Jul 2023, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units  
from RCAIDE.Framework.Networks.Series_Hybrid_Network               import Series_Hybrid_Network
from RCAIDE.Library.Methods.Propulsors.Turboelectric_Propulsor   import design_turboelectric_turbine
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor          import design_propeller 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor       import design_motor 
from RCAIDE.Library.Methods.Propulsors.Converters.Generator      import design_generator
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion     import nasa_motor
from RCAIDE.Library.Methods.Energy.Sources.Battery.Common               import initialize_from_circuit_configuration
from RCAIDE.Library.Methods.Geometry.Two_Dimensional.Planform           import wing_segmented_planform 
from RCAIDE.Library.Plots import * 

# python imports 
import numpy as np 
from copy import deepcopy
import os
import matplotlib.pyplot        as plt 
 
def main():     
      
    # vehicle data
    vehicle  = vehicle_setup() 

    # Set up vehicle configs
    configs  = configs_setup(vehicle)

    # create analyses
    analyses = analyses_setup(configs)

    # mission analyses 
    mission = mission_setup(analyses)

    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 

    # mission analysis 
    results = missions.base_mission.evaluate()  

    # plot the results
    plot_mission(results)      
           
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
    vehicle.envelope.ultimate_load        = 5.7
    vehicle.envelope.limit_load           = 3.8 
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
    rel_path                              = os.path.dirname(ospath) + separator + '..' +  separator + '..' +  separator 
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
    segment.sweeps.quarter_chord          = 82. * Units.degrees 
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
    fuselage = RCAIDE.Library.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage'
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
 
    # ##########################################################   Nacelles  ############################################################    
    nacelle                    = RCAIDE.Library.Components.Nacelles.Nacelle()
    nacelle.tag                = 'nacelle_1'
    nacelle.length             = 2
    nacelle.diameter           = 42 * Units.inches
    nacelle.areas.wetted       = 0.01*(2*np.pi*0.01/2)
    nacelle.origin             = [[2.5,2.5,1.0]]
    nacelle.flow_through       = False  
    
    nac_segment                    = RCAIDE.Library.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.1  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.3  
    nac_segment.height             = 0.52
    nac_segment.width              = 0.7
    nacelle.append_segment(nac_segment)  
     
    nac_segment                    = RCAIDE.Library.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.5  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.7 
    nac_segment.height             = 0.4
    nac_segment.width              = 0.6
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.9 
    nac_segment.height             = 0.3
    nac_segment.width              = 0.5
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 1.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)    
    
    vehicle.append_component(nacelle)  

    nacelle_2          = deepcopy(nacelle)
    nacelle_2.tag      = 'nacelle_2'
    nacelle_2.origin   = [[2.5,-2.5,1.0]]
    vehicle.append_component(nacelle_2)    
 
    # ########################################################  Energy Network  #########################################################  
    net                              = Series_Hybrid_Network()
 
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                                                    = RCAIDE.Energy.Networks.Distribution.Electrical_Bus()  
    
    
    #------------------------------------------------------------------------------------------------------------------------------------           
    # Battery1
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat                                                    = RCAIDE.Energy.Sources.Batteries.Lithium_Ion_NMC()
    bat.tag                                                = 'li_ion_battery_1'
    bat.pack.electrical_configuration.series               = 140   
    bat.pack.electrical_configuration.parallel             = 100
    initialize_from_circuit_configuration(bat)  
    bat.pack.number_of_modules                           = 14  
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module.voltage                                     = bat.pack.maximum_voltage/bat.pack.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
    bat.module.geometrtic_configuration.normal_count       = 24
    bat.module.geometrtic_configuration.parallel_count     = 40
    bat.thermal_management_system.heat_acquisition_system  = RCAIDE.Energy.Thermal_Management.Batteries.Heat_Acquisition_Systems.Direct_Air()      
    bat.center_of_gravity = np.array([[0, 0, 0]])
    bus.voltage                                            = bat.pack.maximum_voltage
    bus.batteries.append(bat)            
    
    #------------------------------------------------------------------------------------------------------------------------------------           
    # Battery2
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat_2                                                  = deepcopy(bat)
    bat_2.tag                                              = 'li_ion_battery_2'
    bat.center_of_gravity = np.array([[0, 0, 1]])
    initialize_from_circuit_configuration(bat_2)        
    bus.batteries.append(bat_2) 
    
    #------------------------------------------------------------------------------------------------------------------------------------           
    # Battery3
    #------------------------------------------------------------------------------------------------------------------------------------  
    bat_3                                                  = deepcopy(bat)
    bat_3.tag                                              = 'li_ion_battery_3'
    bat.center_of_gravity = np.array([[0, 0, -1]])
    initialize_from_circuit_configuration(bat_3)        
    bus.batteries.append(bat_3)          

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    starboard_propulsor                              = RCAIDE.Energy.Propulsors.Electric_Rotor()  
    starboard_propulsor.tag                          = 'starboard_propulsor'
    starboard_propulsor.active_batteries             = ['li_ion_battery_1','li_ion_battery_2']   
  
    # Electronic Speed Controller       
    esc                                              = RCAIDE.Energy.Propulsors.Modulators.Electronic_Speed_Controller()
    esc.tag                                          = 'esc_1'
    esc.efficiency                                   = 0.95 
    starboard_propulsor.electronic_speed_controller  = esc   
     
    # Propeller              
    propeller                                        = RCAIDE.Energy.Propulsors.Converters.Propeller() 
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
    propeller                                        = design_propeller(propeller)    
    starboard_propulsor.rotor                        = propeller   
              
    # DC_Motor       
    motor                                            = RCAIDE.Energy.Propulsors.Converters.DC_Motor()
    motor.tag                                        = 'motor_1'
    motor.efficiency                                 = 0.98
    motor.origin                                     = [[2.,  2.5, 0.95]]
    motor.nominal_voltage                            = bat.pack.maximum_voltage
    motor.no_load_current                            = 1
    motor.rotor_radius                               = propeller.tip_radius
    motor.design_torque                              = propeller.cruise.design_torque
    motor.angular_velocity                           = propeller.cruise.design_angular_velocity 
    motor                                            = design_motor(motor)  
    motor.mass_properties.mass                       = nasa_motor(motor.design_torque) 
    starboard_propulsor.motor                        = motor 
 
    # append propulsor to distribution line 
    bus.propulsors.append(starboard_propulsor) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Port Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Energy.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_batteries            = ['li_ion_battery_2','li_ion_battery_3']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.tag                                  = 'esc_2'
    esc_2.origin                               = [[2., -2.5, 0.95]]      
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         = [[2.,-2.5,0.95]]
    propeller_2.clockwise_rotation             = False        
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.tag                                = 'motor_2'
    motor_2.origin                             = [[2., -2.5, 0.95]]      
    port_propulsor.motor                       = motor_2  
    
    # append propulsor to distribution line 
    bus.propulsors.append(port_propulsor) 


    #------------------------------------------------------------------------------------------------------------------------------------           
    # Payload 
    #------------------------------------------------------------------------------------------------------------------------------------  
    payload                      = RCAIDE.Energy.Peripherals.Payload()
    payload.power_draw           = 10. # Watts
    payload.mass_properties.mass = 1.0 * Units.kg
    bus.payload                  = payload

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Avionics
    #------------------------------------------------------------------------------------------------------------------------------------  
    avionics                     = RCAIDE.Energy.Peripherals.Avionics()
    avionics.power_draw          = 20. # Watts
    bus.avionics                 = avionics   

    # append bus   
    net.busses.append(bus)
    
    
    
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Fuel Distrubition Line 
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                     = RCAIDE.Energy.Networks.Distribution.Fuel_Line() 
    fuel_line.identical_propulsors                = False # only for regression
    
        
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Fuel Tank & Fuel
    #------------------------------------------------------------------------------------------------------------------------------------   
    fuel_tank                                      = RCAIDE.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_1'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[26.5,0,0]]) # Needs to be changed
    fuel_tank.mass_properties.fuel_mass_when_full  = 11096
    fuel_tank.fuel_selector_ratio                  = 1/2
    fuel_tank.fuel_type                            = RCAIDE.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_2'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[28.7,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 11943
    fuel_tank.fuel_selector_ratio                  = 1/2
    fuel_tank.fuel_type                            = RCAIDE.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
    
    fuel_tank                                      = RCAIDE.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                                  = 'tank_3'
    fuel_tank.mass_properties.center_of_gravity    = np.array([[31.0,0,0]])
    fuel_tank.mass_properties.fuel_mass_when_full  = 4198+4198
    fuel_tank.fuel_selector_ratio                  = 1/2
    fuel_tank.fuel_type                            = RCAIDE.Attributes.Propellants.Jet_A() 
    fuel_line.fuel_tanks.append(fuel_tank) 
            
    
     # Append fuel line to network      
    net.fuel_lines.append(fuel_line)      
    
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Right Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    right_turboelectric_generator                   = RCAIDE.Energy.Propulsors.Turboelectric_Generator()  
    right_turboelectric_generator.tag               = 'right_propulsor'   
    right_turboelectric_generator.active_fuel_tanks = ['tank_1','tank_2']
      
    right_gas_turbine                               = RCAIDE.Energy.Propulsors.Converters.Turbine() 
    right_gas_turbine.tag                           = 'right_turbine'
    right_gas_turbine.engine_length                 = 4.039
    right_gas_turbine.mechanical_efficiency         = 0.99
    right_gas_turbine.nacelle_diameter              = 1.3
    right_gas_turbine.inlet_diameter                = 1.212 
    right_gas_turbine.areas.wetted                  = 30
    right_gas_turbine.mass_properties.mass          = 0.3
    right_gas_turbine.design_altitude               = 60000.0*Units.ft
    right_gas_turbine.design_mach_number            = 2.02
    #right_gas_turbine.design_thrust                 = 40000. * Units.lbf 
    right_gas_turbine.origin                        = [[37.,6.,-1.3]] 
    right_gas_turbine.working_fluid                 = RCAIDE.Attributes.Gases.Air() 
    right_gas_turbine.shaft_angular_velocity        = 2700. * Units.rpm
    
    
    # Ram  
    ram     = RCAIDE.Energy.Propulsors.Converters.Ram()
    ram.tag = 'ram' 
    right_gas_turbine.append(ram) 
 
    # Inlet Nozzle 
    inlet_nozzle                          = RCAIDE.Energy.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                      = 'inlet_nozzle' 
    inlet_nozzle.polytropic_efficiency    = 1.0
    inlet_nozzle.pressure_ratio           = 1.0
    inlet_nozzle.pressure_recovery        = 0.94
    right_gas_turbine.append(inlet_nozzle)     
          
    #  Low Pressure Compressor      
    lp_compressor                            = RCAIDE.Energy.Propulsors.Converters.Compressor()    
    lp_compressor.tag                        = 'low_pressure_compressor' 
    lp_compressor.polytropic_efficiency      = 0.88
    lp_compressor.pressure_ratio             = 3.1     
    right_gas_turbine.append(lp_compressor)
    
    # High Pressure Compressor        
    hp_compressor                                 = RCAIDE.Energy.Propulsors.Converters.Compressor()    
    hp_compressor.tag                             = 'high_pressure_compressor' 
    hp_compressor.polytropic_efficiency           = 0.88
    hp_compressor.pressure_ratio                  = 5.0      
    right_gas_turbine.append(hp_compressor)
    
    # Low Pressure Turbine 
    lp_turbine                               = RCAIDE.Energy.Propulsors.Converters.Turbine()   
    lp_turbine.tag                           ='low_pressure_turbine' 
    lp_turbine.mechanical_efficiency         = 0.99
    lp_turbine.polytropic_efficiency         = 0.89
    right_gas_turbine.append(lp_turbine)
                 
    # High Pressure Turbine         
    hp_turbine                               = RCAIDE.Energy.Propulsors.Converters.Turbine()   
    hp_turbine.tag                           ='high_pressure_turbine' 
    hp_turbine.mechanical_efficiency         = 0.99
    hp_turbine.polytropic_efficiency         = 0.87 
    right_gas_turbine.append(hp_turbine)
          
    # Combustor   
    combustor                             = RCAIDE.Energy.Propulsors.Converters.Combustor()   
    combustor.tag                         = 'combustor' 
    combustor.efficiency                  = 0.94
    combustor.alphac                      = 1.0     
    combustor.turbine_inlet_temperature   = 1440.
    combustor.pressure_ratio              = 0.92
    combustor.fuel_data                   = RCAIDE.Attributes.Propellants.Jet_A()     
    right_gas_turbine.append(combustor)
 
    # Core Nozzle 
    nozzle                                = RCAIDE.Energy.Propulsors.Converters.Supersonic_Nozzle()   
    nozzle.tag                            = 'core_nozzle' 
    nozzle.pressure_recovery              = 0.95
    nozzle.pressure_ratio                 = 1.    
    right_gas_turbine.append(nozzle) 
    

    # design turboelectric 
    right_gas_turbine = design_turboelectric_turbine(right_gas_turbine)   
    right_turboelectric_generator.turbines.append(right_gas_turbine)
    
    right_gas_turbine.shaft_torque              = right_gas_turbine.design_power[0][0] * 0.737 / right_gas_turbine.shaft_angular_velocity
    generator                                  = RCAIDE.Energy.Propulsors.Converters.Generator()
    generator.tag                              = 'right_generator'
    generator.origin                           = [[35.,6.,-1.3]] 
    generator.connected_engine                 = ['right_turbine'] 
    generator.mass_properties.mass             = motor.mass_properties.mass
     
    net.hybrid_electrical_power_split_ratio 
    generator.no_load_torque                   = 20
    generator.efficiency                       = 0.95    
    generator.nominal_voltage                  = bat.pack.maximum_voltage
    
    generator.shaft_radius = 0.1
    right_gas_turbine.shaft_radius = 0.08
    
    
    
    generator.shaft_input_power                = right_gas_turbine.low_pressure_turbine.inputs.shaft_power_off_take * right_gas_turbine.mechanical_efficiency
    generator.design_power = 60000 #net.generator_power  # generator.design_power = generator.shaft_input_power * generator.efficiency
    
    generator.no_load_current = 1
    generator.current = generator.design_power / generator.nominal_voltage
    
    generator.design_torque = 400    #right_gas_turbine.shaft_torque * generator.shaft_radius / right_gas_turbine.shaft_radius
    generator.design_omega = generator.design_power/generator.design_torque
    
    #generator.current = generator.design_power / generator.nominal_voltag
    
    generator                                  = design_generator(generator)
    
    right_turboelectric_generator.generators.append(generator)
    
    # append turboelectric    
    fuel_line.propulsors.append(right_turboelectric_generator)
    
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Left Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    left_turboelectric_generator                    = deepcopy(right_turboelectric_generator) 
    left_turboelectric_generator.tag                 = 'left_turbine'       
    left_turboelectric_generator.origin              = [[37.,5.3,-1.3]]   
    left_turboelectric_generator.active_fuel_tanks   = ['tank_2','tank_3']    
    fuel_line.propulsors.append(left_turboelectric_generator)      
    
    vehicle.append_energy_network(net)

    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------
    
    return vehicle
# ---------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):

    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'  
    configs.append(base_config) 

    # done!
    return configs
 
 
# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle()
 
    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Framework.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle = vehicle 
    analyses.append(aerodynamics)  
 
    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Framework.Analyses.Energy.Energy()
    energy.networks= vehicle.networks 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   
    

    # done!
    return analyses    

# ----------------------------------------------------------------------
#   Define the Mission
# ---------------------------------------------------------------------- 

def mission_setup(analyses):   
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission' 

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment() 
    
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'DER'       
    segment.analyses.extend( analyses.base )
    segment.initial_battery_state_of_charge                  = 0.89   
    segment.altitude_start                                   = 0.0 * Units.feet
    segment.altitude_end                                     = 50.0 * Units.feet
    segment.air_speed_start                                  = 45  * Units['m/s'] 
    segment.air_speed_end                                    = 45        
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Initial Climb Area Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'ICA' 
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                   = 50.0 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = 45  * Units['m/s']   
    segment.air_speed_end                                    = 50 * Units['m/s']   
    segment.climb_rate                                       = 600 * Units['ft/min']    
    mission.append_segment(segment) 
             
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Climb'        
    segment.analyses.extend( analyses.base )      
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 2500 * Units.feet
    segment.air_speed                                        = 120 * Units['mph']
    segment.climb_rate                                       = 500* Units['ft/min']   
    mission.append_segment(segment) 
    
    # ------------------------------------------------------------------
    #   Climb 1 : constant Speed, constant rate segment 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "Ascent"
    segment.analyses.extend( analyses.base )
    segment.initial_battery_state_of_charge  = 0.89 
    segment.altitude_start                   = 2500.0  * Units.feet
    segment.altitude_end                     = 8012    * Units.feet 
    segment.air_speed                        = 96.4260 * Units['mph'] 
    segment.climb_rate                       = 700.034 * Units['ft/min']      
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Cruise Segment: constant Speed, constant altitude
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "Cruise" 
    segment.analyses.extend(analyses.base) 
    segment.altitude                  = 8012   * Units.feet
    segment.air_speed                 = 120.91 * Units['mph'] 
    segment.distance                  = 50.   * Units.nautical_mile        
    mission.append_segment(segment)    


    # ------------------------------------------------------------------
    #   Descent Segment Flight 1   
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = "Decent"  
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                   = 8012 * Units.feet  
    segment.altitude_end                                     = 1000.0 * Units.feet
    segment.air_speed_start                                  = 130.* Units['mph']  
    segment.air_speed_end                                    = 120 * Units['mph']   
    segment.climb_rate                                       = -200 * Units['ft/min']     
    mission.append_segment(segment)   

    # ------------------------------------------------------------------
    #  Downleg_Altitude Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag = 'Downleg'
    segment.analyses.extend(analyses.base) 
    segment.air_speed                                        = 120 * Units['mph']           
    segment.distance                                         = 6000 * Units.feet
    segment.acceleration                                     = -0.025  * Units['m/s/s']   
    segment.descent_rate                                     = 300 * Units['ft/min']           
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #  Reserve Climb 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Reserve_Climb'        
    segment.analyses.extend( analyses.base )      
    segment.altitude_end                                     = 1500 * Units.feet
    segment.air_speed                                        = 120 * Units['mph']
    segment.climb_rate                                       = 500* Units['ft/min']   
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #  Researve Cruise Segment 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag = 'Reserve_Cruise'  
    segment.analyses.extend(analyses.base) 
    segment.air_speed                                        = 145* Units['mph']
    segment.distance                                         = 60 * Units.miles * 0.1  
    mission.append_segment(segment)     
    
    # ------------------------------------------------------------------
    #  Researve Descent
    # ------------------------------------------------------------------ 
    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Reserve_Descent'
    segment.analyses.extend( analyses.base )    
    segment.altitude_end                                     = 1000 * Units.feet 
    segment.air_speed                                        = 110 * Units['mph']
    segment.descent_rate                                     = 300 * Units['ft/min']   
    mission.append_segment(segment)  

    
    # ------------------------------------------------------------------
    #  Baseleg Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = 'Baseleg'
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                                   = 1000 * Units.feet
    segment.altitude_end                                     = 500.0 * Units.feet
    segment.air_speed_start                                  = 45  
    segment.air_speed_end                                    = 40    
    segment.climb_rate                                       = -350 * Units['ft/min'] 
    mission.append_segment(segment) 

    # ------------------------------------------------------------------
    #  Final Approach Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment_name = 'Final_Approach'
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)      
    segment.altitude_start                                   = 500.0 * Units.feet
    segment.altitude_end                                     = 00.0 * Units.feet
    segment.air_speed_start                                  = 40 
    segment.air_speed_end                                    = 35   
    segment.climb_rate                                       = -300 * Units['ft/min']        
    mission.append_segment(segment)  
 
    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------ 
    
    return mission 
 

def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions 


def plot_mission(results):  
    
    plot_flight_conditions(results) 

    plot_aerodynamic_coefficients(results)  
    
    plot_aircraft_velocities(results)
    
    plot_battery_pack_conditions(results)
    
    plot_battery_cell_conditions(results)
    
    plot_battery_degradation(results)

    plot_rotor_conditions(results) 
     
    return
 
if __name__ == '__main__': 
    main()    
    plt.show()


