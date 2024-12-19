 

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
# RCAIDE imports 
import RCAIDE 
from RCAIDE.Framework.Core import Units   
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor             import design_propeller 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor          import design_motor 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion import compute_motor_weight
from RCAIDE.Library.Methods.Geometry.Planform                       import wing_segmented_planform 
from RCAIDE.Library.Plots                                           import *     

# python imports 
import os
import numpy as np 
from copy import deepcopy
import pickle
import  pandas as pd
import matplotlib.pyplot as plt  

 
def main():
    
    # vehicle data
    new_geometry = True
    if new_geometry :
        vehicle  = vehicle_setup()
        save_aircraft_geometry(vehicle , 'NASA_X57')
    else: 
        vehicle = load_aircraft_geometry('NASA_X57')
        
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
     
    save_filename = 'NASA_X57'
    save_csv      =  True 
    write_results_to_csv(results, save_filename, save_csv) 
           
    return 



# ----------------------------------------------------------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------------------------------------------------------
def vehicle_setup():
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------

    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'X57_Maxwell_Mod2' 
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

         
    #------------------------------------------------------------------------------------------------------------------------------------
    # ######################################################## Wings ####################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------
    
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
    rel_path                              = os.path.dirname(ospath) + separator + '..'    + separator   + '..'    + separator 
    airfoil                               = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                           = 'NACA_63_412.txt' 
    airfoil.coordinate_file               = rel_path + 'Airfoils' + separator + 'NACA_63_412.txt'   # absolute path     
    cg_x                                  = wing.origin[0][0] + 0.25*wing.chords.mean_aerodynamic
    cg_z                                  = wing.origin[0][2] - 0.2*wing.chords.mean_aerodynamic
    vehicle.mass_properties.center_of_gravity = [[cg_x,   0.  ,  cg_z ]]  

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
 
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Electric Network
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
    bat                                                    = RCAIDE.Library.Components.Energy.Sources.Battery_Modules.Lithium_Ion_NMC()
    number_of_modules                                      = 8 
    bat.tag                                                = 'li_ion_battery'
    bat.electrical_configuration.series                    = 16   
    bat.electrical_configuration.parallel                  = 40
    bat.cell.maximum_voltage                               = 4.2                                                                          
    bat.cell.nominal_capacity                              = 3.0                                                                          
    bat.cell.nominal_voltage                               = 3.6   
    bat.geometrtic_configuration.total                      = bat.electrical_configuration.total
    bat.voltage                                             = bat.maximum_voltage 
    bat.geometrtic_configuration.normal_count               = 20
    bat.geometrtic_configuration.parallel_count             = 32
     
    for _ in range(number_of_modules):
        bus.battery_modules.append(deepcopy(bat))    
    
    bus.battery_module_electric_configuration = 'Series' 
    bus.initialize_bus_properties()      
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    starboard_propulsor                              = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    starboard_propulsor.tag                          = 'starboard_propulsor' 
  
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
    motor.nominal_voltage                            = bus.voltage*0.5
    motor.no_load_current                            = 1
    motor.rotor_radius                               = propeller.tip_radius
    motor.design_torque                              = propeller.cruise.design_torque
    motor.angular_velocity                           = propeller.cruise.design_angular_velocity 
    design_motor(motor)  
    motor.mass_properties.mass                       = compute_motor_weight(motor.design_torque) 
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

    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'  
    configs.append(base_config) 
 
    # ------------------------------------------------------------------
    #   Hover Climb Configuration
    # ------------------------------------------------------------------
    config                                      = RCAIDE.Library.Components.Configs.Config(vehicle)
    config.tag                                  = 'takeoff' 
    for network in  config.networks: 
        for bus in network.busses: 
            for propulsor in  bus.propulsors: 
                propulsor.rotor.pitch_command   = -5 *Units.degrees 
    configs.append(config)
    
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
    weights = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
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
    energy.vehicle =  vehicle  
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
    
    ## ------------------------------------------------------------------------------------------------------------------------------------ 
    ##   Takeoff Roll
    ## ------------------------------------------------------------------------------------------------------------------------------------ 

    #segment = Segments.Ground.Takeoff(base_segment)
    #segment.tag = "Takeoff" 
    #segment.analyses.extend( analyses.takeoff )
    #segment.velocity_start                                                = 50 * Units['mph']
    #segment.velocity_end                                                  = 105 * Units['mph']
    #segment.friction_coefficient                                          = 0.04
    #segment.altitude                                                      = 0.0   
    #segment.initial_battery_state_of_charge                               = 1.0 
    #segment.assigned_control_variables.elapsed_time.active                = True  
    #segment.assigned_control_variables.elapsed_time.initial_guess_values  = [[20.]]  
    #mission.append_segment(segment) 
    
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'DER'       
    segment.analyses.extend( analyses.base )  
    segment.altitude_start                                           = 0.0 * Units.feet
    segment.altitude_end                                             = 50.0 * Units.feet 
    segment.air_speed_start                                          = 100 * Units['mph'] 
    segment.air_speed_end                                            = 105 * Units['mph']  
    segment.initial_battery_state_of_charge                               = 1.0 
                       
    # define flight dynamics to model             
    segment.flight_dynamics.force_x                                  = True  
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                  
      
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Initial Climb Area Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'ICA' 
    segment.analyses.extend( analyses.base )   
    segment.altitude_start                                           = 50.0 * Units.feet
    segment.altitude_end                                             = 500.0 * Units.feet 
    segment.air_speed_end                                            = 110 * Units['mph']
    segment.climb_rate                                               = 600 * Units['ft/min']   
               
    # define flight dynamics to model            
    segment.flight_dynamics.force_x                                  = True  
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                  
          
    mission.append_segment(segment)  
   
             
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
    segment.tag = 'climb_1'        
    segment.analyses.extend( analyses.base )      
    segment.altitude_start                                           = 500.0 * Units.feet
    segment.altitude_end                                             = 2500 * Units.feet
    segment.air_speed                                                = 120 * Units['mph']
    segment.climb_rate                                               = 500* Units['ft/min']  
               
    # define flight dynamics to model            
    segment.flight_dynamics.force_x                                  = True  
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                 
           
    mission.append_segment(segment)
    
        
    # ------------------------------------------------------------------
    #   Climb 1 : constant Speed, constant rate segment 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2"
    segment.analyses.extend( analyses.base )
    segment.altitude_start                                = 2500.0  * Units.feet
    segment.altitude_end                                  = 5000    * Units.feet 
    segment.air_speed                                     = 96.4260 * Units['mph'] 
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
    segment.tag = "cruise" 
    segment.analyses.extend(analyses.base) 
    segment.altitude                                      = 5000   * Units.feet
    segment.air_speed                                     = 120.91 * Units['mph'] 
    segment.distance                                      = 30.   * Units.nautical_mile  
    
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
    segment.tag = "decent"  
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                = 5000  * Units.feet  
    segment.altitude_end                                  = 1000 * Units.feet  
    segment.air_speed_end                                 = 110 * Units['mph']   
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
    #  Downleg_Altitude Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment) 
    segment.tag = 'Downleg'
    segment.analyses.extend(analyses.base)   
    segment.air_speed_end                                 = 45.0 * Units['m/s']            
    segment.distance                                      = 6000 * Units.feet
    segment.acceleration                                  = -0.025  * Units['m/s/s']   
    segment.descent_rate                                  = 300 * Units['ft/min']   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                   
            
    mission.append_segment(segment)     
    
    ## ------------------------------------------------------------------
    ##  Reserve Climb 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'Reserve_Climb'        
    #segment.analyses.extend( analyses.base )      
    #segment.altitude_end                                  = 1500 * Units.feet
    #segment.air_speed                                     = 120 * Units['mph']
    #segment.climb_rate                                    = 500* Units['ft/min']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
        
    #mission.append_segment(segment)
    
    ## ------------------------------------------------------------------
    ##  Researve Cruise Segment 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    #segment.tag = 'Reserve_Cruise'  
    #segment.analyses.extend(analyses.base) 
    #segment.air_speed                                     = 145* Units['mph']
    #segment.distance                                      = 60 * Units.miles * 0.1  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                  
       
    #mission.append_segment(segment)     
    
    ## ------------------------------------------------------------------
    ##  Researve Descent
    ## ------------------------------------------------------------------ 
    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'Reserve_Descent'
    #segment.analyses.extend( analyses.base )    
    #segment.altitude_end                                  = 1000 * Units.feet 
    #segment.air_speed                                     = 110 * Units['mph']
    #segment.descent_rate                                  = 300 * Units['ft/min']   
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    #mission.append_segment(segment)  

    
    # ------------------------------------------------------------------
    #  Baseleg Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    segment.tag = 'Baseleg'
    segment.analyses.extend( analyses.base)   
    segment.altitude_start                                = 1000 * Units.feet
    segment.altitude_end                                  = 500.0 * Units.feet
    segment.air_speed_start                               = 45 
    segment.air_speed_end                                 = 40    
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
    segment_name = 'Final_Approach'
    segment.tag = segment_name          
    segment.analyses.extend( analyses.base)      
    segment.altitude_start                                           = 500.0 * Units.feet
    segment.altitude_end                                             = 0.0 * Units.feet
    segment.air_speed_start                                          = 40 
    segment.air_speed_end                                            = 35   
    segment.climb_rate                                               = -300 * Units['ft/min']   
                
    # define flight dynamics to model             
    segment.flight_dynamics.force_x                                  = True  
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                      
    mission.append_segment(segment)   
 
 
    ## ------------------------------------------------------------------------------------------------------------------------------------ 
    ##   Landing Roll
    ## ------------------------------------------------------------------------------------------------------------------------------------ 

    #segment = Segments.Ground.Landing(base_segment)
    #segment.tag = "Landing" 
    #segment.analyses.extend( analyses.base ) 
    #segment.velocity_start                                                = 35 * Units['m/s']
    #segment.velocity_end                                                  = 10 * Units.knots 
    #segment.friction_coefficient                                          = 0.4
    #segment.altitude                                                      = 0.0   
    #segment.assigned_control_variables.elapsed_time.active                = True  
    #segment.assigned_control_variables.elapsed_time.initial_guess_values  = [[30.]]  
    #mission.append_segment(segment)     
 
    
    return mission 
 

def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
  
    return missions 



# ----------------------------------------------------------------------
#   Plot Results
# ----------------------------------------------------------------------

def plot_results(results):
    # Plots fligh conditions 
    #plot_flight_conditions(results) 
    
    # Plot arcraft trajectory
    #plot_flight_trajectory(results)   

    plot_propulsor_throttles(results)
    
    # Plot Aircraft Electronics
    #plot_battery_module_conditions(results) 
    #plot_battery_temperature(results)
    plot_battery_cell_conditions(results) 
    #plot_battery_module_C_rates(results)
    #plot_battery_degradation(results) 
    
    # Plot Propeller Conditions 
    #plot_rotor_conditions(results) 
    #plot_disc_and_power_loading(results)
    
    # Plot Electric Motor and Propeller Efficiencies 
    #plot_electric_propulsor_efficiencies(results)
    
      
    return 


def write_results_to_csv(results, save_filename, save_csv):
    time        = []
    current     = []
    power       = []
    voltage     = [] 
    temperature = [] 
    
    for network in results.segments[0].analyses.energy.vehicle.networks: 
        busses  = network.busses
        for bus in busses: 
            for b_i, battery in enumerate(bus.battery_modules):
                if b_i == 0 or bus.identical_batteries == False: 
                    for i in range(len(results.segments)):  
                        battery_conditions   = results.segments[i].conditions.energy[bus.tag].battery_modules[battery.tag]  
                        time.append(results.segments[i].conditions.frames.inertial.time[:,0])   
                        power.append(battery_conditions.cell.power[:,0])  
                        voltage.append(battery_conditions.cell.voltage_under_load[:,0])  
                        current.append(battery_conditions.cell.current[:,0])
                        temperature.append(battery_conditions.cell.temperature[:,0]) 
                    
    voltage       = np.hstack(voltage)
    current       = np.hstack(current)
    time          = np.hstack(time)
    power         = np.hstack(power)
    temperature   = np.hstack(temperature)
    df      = pd.DataFrame({'Time': time, 'Current': current, 'Voltage': voltage, 'Power': power, 'Temperature': temperature}) 
    if save_csv:
        df.to_csv(save_filename + '.csv', index=False)                       
    return
 
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
 
 
if __name__ == '__main__': 
    main()    
    plt.show()


