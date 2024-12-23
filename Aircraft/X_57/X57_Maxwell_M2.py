# X57_Maxwell.py
#
# Created: Feb 2020, M. Clarke
#          Sep 2020, M. Clarke 

""" setup file for the X57-Maxwell Electric Aircraft 
"""
 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import RCAIDE
from RCAIDE.Core import Units 
import numpy as np   
import matplotlib.pyplot        as plt
import os
import pickle
from RCAIDE.Core                                         import Data 
from RCAIDE.Visualization.Performance.Aerodynamics.Vehicle                 import *  
from RCAIDE.Visualization.Performance.Mission                              import *  
from RCAIDE.Visualization.Performance.Aerodynamics.Rotor import *  
from RCAIDE.Visualization.Performance.Energy.Battery                       import *   
from RCAIDE.Visualization.Performance.Noise                                import *  
from RCAIDE.Visualization.Geometry.Three_Dimensional.plot_3d_vehicle       import plot_3d_vehicle 
from RCAIDE.Visualization.Geometry                                         import *
from RCAIDE.Components.Energy.Networks.Battery_Electric_Rotor                       import Battery_Electric_Rotor
from RCAIDE.Methods.Propulsion                                             import propeller_design 
from RCAIDE.Methods.Power.Battery.Sizing                                   import initialize_from_mass
#from RCAIDE.Methods.Power.Battery.Sizing                                   import initialize_from_circuit_configuration 
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform                      import segment_properties
from RCAIDE.Methods.Propulsion.electric_motor_sizing                       import size_optimal_motor
from RCAIDE.Methods.Weights.Correlations.Propulsion                        import nasa_motor
from RCAIDE.Methods.Weights.Buildups.eVTOL.empty                           import empty

import vsp 
from RCAIDE.Input_Output.OpenVSP.vsp_write import write 

#try:
    #import vsp 
    #from RCAIDE.Input_Output.OpenVSP.vsp_write import write 
#except ImportError:
    ## This allows RCAIDE to build without OpenVSP
    #pass  
from copy import deepcopy

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():   
    configs, analyses = full_setup() 

    configs.finalize()
    analyses.finalize()  

    # weight analysis
    #weight  = General_Aviation.empty(configs.base)	 
    
    # mission analysis
    mission = analyses.missions.base
    results = mission.evaluate() 
    
    # save results  
    #save_results(results)
    
    # plot the results
    plot_results(results)  
    return


# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup():

    # vehicle data
    vehicle  = vehicle_setup() 
    write(vehicle, "X57_Maxwell_M2")  
    
    # Set up configs
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(configs_analyses,vehicle)
    missions_analyses = missions_setup(mission)

    analyses = RCAIDE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    analyses.missions = missions_analyses

    return configs, analyses

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

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
    weights = RCAIDE.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.settings.plot_vortex_distribution  = True     
    aerodynamics.vehicle = vehicle 
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)

    #aerodynamics = RCAIDE.Analyses.Aerodynamics.AVL()
    #aerodynamics.process.compute.lift.inviscid.settings.filenames.avl_bin_name =  '/Users/matthewclarke/Documents/AVL/avl3.35'
    #aerodynamics.settings.print_output = True 
    #aerodynamics.settings.keep_files = True  
    #aerodynamics.vehicle = vehicle
    #analyses.append(aerodynamics)


    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = RCAIDE.Analyses.Stability.Fidelity_Zero()
    stability.geometry = vehicle
    analyses.append(stability)

    #stability = RCAIDE.Analyses.Stability.AVL() 
    #stability.settings.filenames.avl_bin_name = '/Users/matthewclarke/Documents/AVL/avl3.35' 
    #stability.settings.print_output = True 
    #stability.settings.keep_files = True 
    #stability.geometry = vehicle
    #analyses.append(stability)
    
    # ------------------------------------------------------------------
    #  Noise Analysis
    noise = RCAIDE.Analyses.Noise.Fidelity_Zero()   
    noise.geometry = vehicle
    analyses.append(noise)

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Analyses.Energy.Energy()
    energy.network = vehicle.networks 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Analyses.Planets.Earth()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses    


def vehicle_setup():

    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------

    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'X57_Maxwell_Mod2'


    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------

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
    atmo                                  = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    freestream                            = atmo.compute_values (0.)
    freestream0                           = atmo.compute_values (altitude)
    mach_number                           = (cruise_speed/freestream.speed_of_sound)[0][0] 
    vehicle.design_dynamic_pressure       = ( .5 *freestream0.density*(cruise_speed*cruise_speed))[0][0]
    vehicle.design_mach_number            =  mach_number
    
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------    
    wing                                  = RCAIDE.Components.Wings.Main_Wing()
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
    airfoil                               = RCAIDE.Components.Airfoils.Airfoil()
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    airfoil.coordinate_file               = '..' + separator + 'Airfoils' + separator + 'NACA_63_412.txt'
    
    cg_x = wing.origin[0][0] + 0.25*wing.chords.mean_aerodynamic
    cg_z = wing.origin[0][2] - 0.2*wing.chords.mean_aerodynamic
    vehicle.mass_properties.center_of_gravity = [[cg_x,   0.  ,  cg_z ]]  # SOURCE: Design and aerodynamic analysis of a twin-engine commuter aircraft

    # Wing Segments
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'inboard'
    segment.percent_span_location         = 0.0 
    segment.twist                         = 3. * Units.degrees   
    segment.root_chord_percent            = 1. 
    segment.dihedral_outboard             = 0.  
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)

    segment                               = RCAIDE.Components.Wings.Segment()
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
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'winglet'
    segment.percent_span_location         = 0.98
    segment.twist                         = 1.  * Units.degrees 
    segment.root_chord_percent            = 0.630
    segment.dihedral_outboard             = 75. * Units.degrees 
    segment.sweeps.quarter_chord          = 55. * Units.degrees 
    segment.thickness_to_chord            = 0.12 
    segment.append_airfoil(airfoil)
    wing.append_segment(segment) 

    segment                               = RCAIDE.Components.Wings.Segment()
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
    wing = segment_properties(wing)           
    
    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------       
    wing                                  = RCAIDE.Components.Wings.Wing()
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


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------ 
    wing                                  = RCAIDE.Components.Wings.Wing()
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
    wing.origin                           = [[6.75 ,0,0.435]]
    wing.aerodynamic_center               = [0.508 ,0,0]  
    wing.vertical                         = True 
    wing.symmetric                        = False
    wing.t_tail                           = False
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    fuselage = RCAIDE.Components.Fuselages.Fuselage()
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
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0
    segment.percent_z_location                  = 0
    segment.height                              = 0.01
    segment.width                               = 0.01
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  = 0.007279116466
    segment.percent_z_location                  = 0.002502014453
    segment.height                              = 0.1669064748
    segment.width                               = 0.2780205877
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.01941097724
    segment.percent_z_location                  = 0.001216095397
    segment.height                              = 0.3129496403
    segment.width                               = 0.4365777215
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.06308567604
    segment.percent_z_location                  = 0.007395489231
    segment.height                              = 0.5841726619
    segment.width                               = 0.6735119903
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.1653761217
    segment.percent_z_location                  = 0.02891281352
    segment.height                              = 1.064028777
    segment.width                               = 1.067200529
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.2426372155
    segment.percent_z_location                  = 0.04214148761
    segment.height                              = 1.293766653
    segment.width                               = 1.183058255
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.2960174029
    segment.percent_z_location                  = 0.04705241831
    segment.height                              = 1.377026712
    segment.width                               = 1.181540054
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 0.3809404284
    segment.percent_z_location                  = 0.05313580461
    segment.height                              = 1.439568345
    segment.width                               = 1.178218989
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_8'
    segment.percent_x_location                  = 0.5046854083
    segment.percent_z_location                  = 0.04655492473
    segment.height                              = 1.29352518
    segment.width                               = 1.054390707
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_9'
    segment.percent_x_location                  = 0.6454149933
    segment.percent_z_location                  = 0.03741966266
    segment.height                              = 0.8971223022
    segment.width                               = 0.8501926505
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_10'
    segment.percent_x_location                  = 0.985107095
    segment.percent_z_location                  = 0.04540283436
    segment.height                              = 0.2920863309
    segment.width                               = 0.2012565415
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_11'
    segment.percent_x_location                  = 1
    segment.percent_z_location                  = 0.04787575562
    segment.height                              = 0.1251798561
    segment.width                               = 0.1206021048
    fuselage.Segments.append(segment)

    # add to vehicle
    vehicle.append_component(fuselage)

    # ------------------------------------------------------------------
    #   Nacelles
    # ------------------------------------------------------------------ 
    nacelle                = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.tag            = 'nacelle_1'
    nacelle.length         = 2
    nacelle.diameter       = 42 * Units.inches
    nacelle.areas.wetted   = 0.01*(2*np.pi*0.01/2)
    nacelle.origin         = [[2.5,2.5,1.0]]
    nacelle.flow_through   = False  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.1  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.3  
    nac_segment.height             = 0.52
    nac_segment.width              = 0.7
    nacelle.append_segment(nac_segment)  
     
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.5  
    nac_segment.height             = 0.5
    nac_segment.width              = 0.65
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.7 
    nac_segment.height             = 0.4
    nac_segment.width              = 0.6
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.9 
    nac_segment.height             = 0.3
    nac_segment.width              = 0.5
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
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
    
    #---------------------------------------------------------------------------------------------
    # DEFINE PROPELLER
    #---------------------------------------------------------------------------------------------
    net = Battery_Electric_Rotor()  
    net.rotor_group_indexes      = [0,0]
    net.motor_group_indexes      = [0,0]
    net.esc_group_indexes        = [0,0]

    # Component 1 the ESC
    esc_1            = RCAIDE.Components.Energy.Distributors.Electronic_Speed_Controller()
    esc_1.tag        = 'esc_1'
    esc_1.efficiency = 0.95 
    net.electronic_speed_controllers.append(esc_1)  

    esc_2            = RCAIDE.Components.Energy.Distributors.Electronic_Speed_Controller()
    esc_2.tag        = 'esc_2'
    esc_2.efficiency = 0.95 
    net.electronic_speed_controllers.append(esc_2)     
    
    # Component 2 the Propeller 
    prop                                  = RCAIDE.Components.Energy.Converters.Propeller()
    prop.tag                              = 'propeller_1'
    prop.number_of_blades                 = 3.0
    prop.tip_radius                       = 1.72/2  
    prop.hub_radius                       = 10.     * Units.inches
    prop.cruise.design_freestream_velocity= 150.   * Units.knots
    prop.cruise.design_angular_velocity   = 2400. * Units.rpm
    prop.cruise.design_Cl                 = 0.8
    prop.cruise.design_altitude           = 9000. * Units.feet  
    prop.cruise.design_power              = 98 * 0.65  * Units.hp # assume 65 BHP at cruise
    prop.origin                           = [[2.,2.5,0.784]]
    prop.rotation                         = -1
    prop.symmetry                         = True
    prop.variable_pitch                   = True 
    airfoil                               = RCAIDE.Components.Airfoils.Airfoil()   
    airfoil.number_of_points              = 102 
    airfoil.coordinate_file               = '..' + separator + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                   = ['..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_50000.txt' ,
                                                  '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_100000.txt' ,
                                                  '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_200000.txt' ,
                                                  '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_500000.txt' ,
                                                  '..' + separator + 'Airfoils' + separator + 'Polars' + separator +'NACA_4412_polar_Re_1000000.txt' ]
    prop.append_airfoil(airfoil)        
    prop.airfoil_polar_stations           = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    prop                                  = propeller_design(prop)  

    prop_left = deepcopy(prop)
    prop_left.tag = 'propeller_2' 
    prop_left.origin   = [[2.,-2.5,0.784]]
    prop_left.rotation = 1
    
    net.rotors.append(prop)
    net.rotors.append(prop_left)


    # Component 3 the Battery
    bat = RCAIDE.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion_LiNiMnCoO2_18650()
    bat.mass_properties.mass = 500. * Units.kg  
    bat.pack.max_voltage     = 400. 
    initialize_from_mass(bat)
    net.battery              = bat
    net.voltage              = bat.pack.max_voltage

    # Component 4 Miscellaneous Systems
    sys = RCAIDE.Components.Systems.System()
    sys.mass_properties.mass = 5 # kg
 
    # Component 5 the Motor  
    motor                         = RCAIDE.Components.Energy.Converters.Motor()
    motor.efficiency              = 0.95
    motor.gearbox_efficiency      = 1.
    motor.origin                  = [[2.,  2.5, 0.95],[2.,  -2.5, 0.95]]
    motor.nominal_voltage         = bat.pack.max_voltage*0.8 
    motor.no_load_current         = 0.1 
    motor.rotor_radius            = prop.tip_radius
    motor.design_torque           = prop.cruise.design_torque
    motor.angular_velocity        = prop.cruise.design_angular_velocity/motor.gear_ratio
    motor                         = size_optimal_motor(motor)  
    motor.mass_properties.mass    = 10. * Units.kg 
    
    # append right motor
    net.motors.append(motor)
    
    # append left motor 
    motor_left = deepcopy(motor)
    motor_left.origin = [[2., -2.5, 0.95]] 
    net.motors.append(motor_left) 

    # Component 6 the Payload
    payload = RCAIDE.Components.Energy.Peripherals.Payload()
    payload.power_draw           = 10. # Watts
    payload.mass_properties.mass = 1.0 * Units.kg
    net.payload                  = payload

    # Component 7 the Avionics
    avionics = RCAIDE.Components.Energy.Peripherals.Avionics()
    avionics.power_draw = 20. # Watts
    net.avionics        = avionics

    # add the solar network to the vehicle
    vehicle.append_component(net)

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

    configs = RCAIDE.Components.Configs.Config.Container()

    base_config = RCAIDE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    base_config.networks.battery_electric_rotor.pitch_command = 0
    configs.append(base_config) 

 
    # ------------------------------------------------------------------
    #   Hover Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag = 'landing' 
    config.networks.battery_electric_rotor.pitch_command                = 0.0 * Units.degrees 
    config.networks.battery_electric_rotor.motors.motor.gear_ratio   = 1.1
    config.networks.battery_electric_rotor.motors.motor2.gear_ratio   = 1.1
    configs.append(config)                  


    # done!
    return configs 

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses 


 

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses,vehicle):   
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # airport
    airport = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   =  0. * Units.ft
    airport.delta_isa  =  0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976()

    mission.airport = airport    

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments 
    
    # base segment
    base_segment = Segments.Segment()
    ones_row     = base_segment.state.ones_row
    base_segment.process.initialize.initialize_battery       = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery  
    base_segment.state.numerics.number_of_control_points        = 4 
    base_segment.battery_age_in_days                         = 1 # optional but added for regression
    base_segment.temperature_deviation                       = 1 # Kelvin #  optional but added for regression
    
    # ------------------------------------------------------------------
    #   Climb 1 : constant Speed, constant rate segment 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1"
    segment.analyses.extend( analyses.base )
    segment.battery_energy                   = vehicle.networks.battery_electric_rotor.battery.pack.max_energy * 0.89
    segment.altitude_start                   = 2500.0  * Units.feet
    segment.altitude_end                     = 8012    * Units.feet 
    segment.air_speed                        = 96.4260 * Units['mph'] 
    segment.climb_rate                       = 700.034 * Units['ft/min']     
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)  
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Cruise Segment: constant Speed, constant altitude
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend(analyses.base) 
    segment.altitude                  = 8012   * Units.feet
    segment.air_speed                 = 120.91 * Units['mph'] 
    segment.distance                  =  20.   * Units.nautical_mile   
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)    
    mission.append_segment(segment)    


    # ------------------------------------------------------------------
    #   Descent Segment Flight 1   
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = "decent"  
    segment.analyses.extend( analyses.base )       
    segment.altitude_start                                   = 8012 * Units.feet  
    segment.altitude_end                                     = 2500.0 * Units.feet
    segment.air_speed_start                                  = 130.* Units['mph']  
    segment.air_speed_end                                    = 110 * Units['mph']   
    segment.climb_rate                                       = -200 * Units['ft/min']  
    segment.state.unknowns.throttle                          = 0.8 * ones_row(1)  
    segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)   
    mission.append_segment(segment)  
    
     
    ## ------------------------------------------------------------------
    ##   Departure End of Runway Segment Flight 1 : 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'DER'       
    #segment.analyses.extend( analyses.base ) 
    #segment.battery_energy                                   = vehicle.networks.battery_electric_rotor.battery.pack.max_energy  
    #segment.state.unknowns.throttle                          = 1.0 * ones_row(1) 
    #segment.altitude_start                                   = 0.0 * Units.feet
    #segment.altitude_end                                     = 50.0 * Units.feet
    #segment.air_speed_start                                  = Vstall  
    #segment.air_speed_end                                    = 45      
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)   
    #mission.append_segment(segment)
    
    ## ------------------------------------------------------------------
    ##   Initial Climb Area Segment Flight 1  
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'ICA' 
    #segment.analyses.extend( analyses.base )  
    #segment.state.unknowns.throttle                          = 0.85  * ones_row(1)  
    #segment.altitude_start                                   = 50.0 * Units.feet
    #segment.altitude_end                                     = 500.0 * Units.feet
    #segment.air_speed_start                                  = 45  * Units['m/s']   
    #segment.air_speed_end                                    = 50 * Units['m/s']   
    #segment.climb_rate                                       = 600 * Units['ft/min']   
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)   
    #mission.append_segment(segment) 
             
    ## ------------------------------------------------------------------
    ##   Climb Segment Flight 1 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'Climb'        
    #segment.analyses.extend( analyses.base )     
    #segment.state.unknowns.throttle                          = 0.85 * ones_row(1)
    #segment.state.unknowns.propeller_power_coefficient       = 0.3  * ones_row(1)   
    #segment.altitude_start                                   = 500.0 * Units.feet
    #segment.altitude_end                                     = 2500 * Units.feet
    #segment.air_speed                                        = 120 * Units['mph']
    #segment.climb_rate                                       = 500* Units['ft/min']  
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment)   
    #mission.append_segment(segment)
    
    ## ------------------------------------------------------------------
    ##   Cruise Segment Flight 1 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    #segment.tag = 'Cruise'  
    #segment.analyses.extend(analyses.base) 
    #segment.air_speed                                        = 145* Units['mph']
    #segment.distance                                         = 60 * Units.miles
    #segment.state.unknowns.throttle                          = 0.6 * ones_row(1)    
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient= 0.12)   
    #mission.append_segment(segment)     
    
    ## ------------------------------------------------------------------
    ##   Descent Segment Flight 1   
    ## ------------------------------------------------------------------ 
    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'Descent'
    #segment.analyses.extend( analyses.base )    
    #segment.altitude_end                                     = 1000 * Units.feet 
    #segment.air_speed                                        = 110 * Units['mph']
    #segment.descent_rate                                     = 300 * Units['ft/min'] 
    #segment.state.unknowns.throttle                          = 0.5  * ones_row(1)   
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient=0.3)  
    #mission.append_segment(segment)  

    ## ------------------------------------------------------------------
    ##  Downleg_Altitude Segment Flight 1 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Cruise.Constant_Acceleration_Constant_Altitude(base_segment) 
    #segment.tag = 'Downleg'
    #segment.analyses.extend(analyses.base) 
    #segment.air_speed_start                                  = 110 * Units['mph'] * Units['m/s']
    #segment.air_speed_end                                    = 45.0 * Units['m/s']            
    #segment.distance                                         = 6000 * Units.feet
    #segment.acceleration                                     = -0.05307 * Units['m/s/s']  
    #segment.air_speed                                        = 49.174
    #segment.descent_rate                                     = 300 * Units['ft/min']
    #segment.state.unknowns.throttle                          = 0.5   * ones_row(1)        
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient=0.3)      
    #mission.append_segment(segment)     
    
    ## ------------------------------------------------------------------
    ##  Reserve Climb 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'Reserve_Climb'        
    #segment.analyses.extend( analyses.base )     
    #segment.state.unknowns.throttle                          = 0.85 * ones_row(1) 
    #segment.altitude_end                                     = 1500 * Units.feet
    #segment.air_speed                                        = 120 * Units['mph']
    #segment.climb_rate                                       = 500* Units['ft/min']  
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient=0.3 )  
    #mission.append_segment(segment)
    
    ## ------------------------------------------------------------------
    ##  Researve Cruise Segment 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    #segment.tag = 'Reserve_Cruise'  
    #segment.analyses.extend(analyses.base) 
    #segment.air_speed                                        = 145* Units['mph']
    #segment.distance                                         = 60 * Units.miles * 0.1
    #segment.state.unknowns.throttle                          = 0.6 * ones_row(1)    
    #segment.state.unknowns.propeller_power_coefficient       = 0.4 * ones_row(1)   
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient=0.4 ) 
    #mission.append_segment(segment)     
    
    ## ------------------------------------------------------------------
    ##  Researve Descent
    ## ------------------------------------------------------------------ 
    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'Reserve_Descent'
    #segment.analyses.extend( analyses.base )    
    #segment.altitude_end                                     = 1000 * Units.feet 
    #segment.air_speed                                        = 110 * Units['mph']
    #segment.descent_rate                                     = 300 * Units['ft/min'] 
    #segment.state.unknowns.throttle                          = 0.5  * ones_row(1)  
    #segment.state.unknowns.propeller_power_coefficient       = 0.3     * ones_row(1)   
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient=0.3) 
    #mission.append_segment(segment)  

    
    ## ------------------------------------------------------------------
    ##  Baseleg Segment Flight 1  
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    #segment.tag = 'Baseleg'
    #segment.analyses.extend( analyses.base)  
    #segment.state.unknowns.throttle                          = 0.8 * ones_row(1) 
    #segment.altitude_start                                   = 1000 * Units.feet
    #segment.altitude_end                                     = 500.0 * Units.feet
    #segment.air_speed_start                                  = 45 
    #segment.air_speed_end                                    = 40    
    #segment.climb_rate                                       = -350 * Units['ft/min']
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient=0.3) 
    #mission.append_segment(segment) 

    ## ------------------------------------------------------------------
    ##  Final Approach Segment Flight 1  
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    #segment_name = 'Final_Approach'
    #segment.tag = segment_name          
    #segment.analyses.extend( analyses.landing)     
    #segment.state.unknowns.throttle                          = 0.8 * ones_row(1) 
    #segment.altitude_start                                   = 500.0 * Units.feet
    #segment.altitude_end                                     = 00.0 * Units.feet
    #segment.air_speed_start                                  = 40 
    #segment.air_speed_end                                    = 35   
    #segment.climb_rate                                       = -300 * Units['ft/min']       
    #segment = vehicle.networks.battery_electric_rotor.add_unknowns_and_residuals_to_segment(segment,initial_power_coefficient=0.3)  
    #mission.append_segment(segment)  
 
    ## ------------------------------------------------------------------
    ##   Mission definition complete    
    ## ------------------------------------------------------------------ 
    
    return mission

def missions_setup(base_mission):

    # the mission container
    missions = RCAIDE.Analyses.Mission.Mission.Container()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------
    missions.base = base_mission

    # done!
    return missions  


def plot_results(results,line_style = 'bo-'):  
    
    
    plot_flight_conditions(results, line_style) 
    
    # Plot Flight Conditions 
    plot_flight_conditions(results, line_style) 
    
    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results, line_style)  
    
    # Plot Aircraft Flight Speed
    plot_aircraft_velocities(results, line_style)
    
    # Plot Aircraft Electronics
    plot_battery_pack_conditions(results, line_style)
    
    # Plot Propeller Conditions 
    plot_rotor_conditions(results, line_style) 
    
    # Plot Electric Motor and Propeller Efficiencies 
    plot_electric_motor_and_rotor_efficiencies(results, line_style)
    
    # Plot propeller Disc and Power Loading
    plot_disc_power_loading(results, line_style)   
    
    plot_stability_coefficients(results)
                        
    return

 
def save_results(results):

    # Store data (serialize)
    with open('Maxwell_M2_Results.pkl', 'wb') as file:
        pickle.dump(results, file)
        
    return   

if __name__ == '__main__': 
    main()    
    plt.show()