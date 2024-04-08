# Navion.py
# 
# Created: Dec 2021, M. Clarke

""" setup file for a mission with a Navion
"""
# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
import RCAIDE 
from RCAIDE.Core import Units   
from RCAIDE.Methods.Energy.Propulsors.Converters.Rotor import design_propeller
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform  import segment_properties
from RCAIDE.Visualization       import *  

# python imports 
import os 
import numpy as np
import pylab as plt
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    run_stability     = True 
    stability_method  = 'vlm' # either "avl" , "vlm" or "analytical"
    
    # vehicle data
    vehicle  = vehicle_setup() 

    # Set up vehicle configs
    configs  = configs_setup(vehicle)

    # create analyses
    analyses = analyses_setup(configs,stability_method,run_stability)

    # mission analyses
    mission  = mission_setup(analyses) 

    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 

    # mission analysis 
    results = missions.base_mission.evaluate()   

    # plt results
    plot_mission(results,run_stability)


    return  
# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs,stability_method,run_stability):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config, stability_method,run_stability, configs)
        analyses[tag] = analysis

    return analyses


def base_analysis(vehicle,stability_method, run_stability,configs):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Analyses.Vehicle() 

    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Analyses.Aerodynamics.Subsonic_VLM() 
    aerodynamics.geometry                            = vehicle
    aerodynamics.settings.number_spanwise_vortices   = 30
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics) 
    
    if run_stability:
        # ------------------------------------------------------------------
        #  Stability Analysis
        if stability_method == "avl":
            stability                                   = RCAIDE.Analyses.Stability.AVL()     
            stability.settings.filenames.avl_bin_name   = 'C:\\Users\\Matteo\\Documents\\UIUC\\avl.exe' #/Users/matthewclarke/Documents/AVL/avl3.35'    
            #stability.settings.trim_aircraft            = False 
            stability.settings.print_output             = False 
        elif stability_method == "vlm":
            stability                                   = RCAIDE.Analyses.Stability.VLM_Perturbation_Method()
            #stability.settings.trim_aircraft            = False 
        elif stability_method == "analytical":
            stability                                   = RCAIDE.Analyses.Stability.Analytical_Approximation()
            #stability.settings.trim_aircraft            = False 
    
        stability.configuration                         = configs
        stability.geometry                              = vehicle
        analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Analyses.Energy.Energy()
    energy.networks = vehicle.networks  
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

    # done!
    return analyses 
  

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup(): 
       # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------ 
    vehicle     = RCAIDE.Vehicle()
    vehicle.tag = 'Navion' 

    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------

    # mass properties
    vehicle.mass_properties.max_takeoff               = 2948 * Units.pounds
    vehicle.mass_properties.takeoff                   = 2948 * Units.pounds
    vehicle.mass_properties.moments_of_inertia.tensor = np.array([[164627.7,0.0,0.0],[0.0,471262.4,0.0],[0.0,0.0,554518.7]])
    vehicle.mass_properties.center_of_gravity         = [[2.239696797,0,-0.131189711 ]]
    vehicle.envelope.ultimate_load                    = 5.7
    vehicle.envelope.limit_load                       = 3.8
    vehicle.reference_area                            = 17.112 
    vehicle.passengers                                = 2 
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------   

    wing                                  = RCAIDE.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.sweeps.quarter_chord             = 0.165 * Units.degrees 
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 17.112
    wing.chords.mean_aerodynamic          = 1.74 
    wing.taper                            = 0.54 
    wing.aspect_ratio                     = 6.04  
    wing.spans.projected                  = 10.166
    wing.chords.root                      = 2.1944 
    wing.chords.tip                       = 1.1850
    wing.twists.root                      = 2 * Units.degrees  
    wing.twists.tip                       = -1 * Units.degrees   
    wing.dihedral                         = 7.5 * Units.degrees   
    wing.origin                           = [[1.652555594, 0.,-0.6006666]]
    wing.aerodynamic_center               = [1.852555594, 0., 6006666 ] # INCORRECT 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0   
    
    tip_airfoil                           = RCAIDE.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = 'Airfoils/NACA_6410.txt' 
 
    root_airfoil                          = RCAIDE.Components.Airfoils.Airfoil()
    root_airfoil.coordinate_file          = 'Airfoils/NACA_4415.txt' 
    
    # Wing Segments 
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'root_segment'
    segment.percent_span_location         = 0.0
    segment.twist                         = 2 * Units.degrees  
    segment.root_chord_percent            = 1.0
    segment.dihedral_outboard             = 7.5 * Units.degrees  
    segment.sweeps.quarter_chord          = 0.165 * Units.degrees  
    segment.thickness_to_chord            = .15 
    wing.append_segment(segment)  
         
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'tip'
    segment.percent_span_location         = 1.0
    segment.twist                         = -1.0 * Units.degrees
    segment.root_chord_percent            = 0.54  
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = 0 * Units.degrees  
    segment.thickness_to_chord            = .12
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)     
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)    
    
                                          
    # control surfaces ------------------------------------------- 
    flap                          = RCAIDE.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.5
    flap.deflection               = 0.0 * Units.degrees 
    flap.chord_fraction           = 0.20
    wing.append_control_surface(flap)  
    

    aileron                       = RCAIDE.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.9 
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.2
    wing.append_control_surface(aileron)      

    # add to vehicle
    vehicle.append_component(wing) 
    

    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------       
    wing                                  = RCAIDE.Components.Wings.Wing()
    wing.tag                              = 'horizontal_stabilizer'  
    wing.sweeps.leading_edge              = 6 * Units.degrees 
    wing.thickness_to_chord               = 12
    wing.areas.reference                  = 4   
    wing.spans.projected                  = 4 
    wing.chords.root                      = 1.2394
    wing.chords.mean_aerodynamic          = 1.0484
    wing.chords.tip                       = 0.8304 
    wing.taper                            = wing.chords.tip/wing.chords.root
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference
    wing.twists.root                      = 0 * Units.degrees  
    wing.twists.tip                       = 0 * Units.degrees   
    wing.origin                           = [[ 6.54518625 , 0., 0.203859697]]
    wing.aerodynamic_center               = [[ 6.545186254 + 0.25*wing.spans.projected, 0., 0.203859697]] 
    wing.vertical                         = False 
    wing.symmetric                        = True
    wing.high_lift                        = False 
    wing.dynamic_pressure_ratio           = 0.9  
    
    elevator                              = RCAIDE.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                          = 'elevator'
    elevator.span_fraction_start          = 0.1
    elevator.span_fraction_end            = 0.9
    elevator.deflection                   = 0.0  * Units.deg
    elevator.chord_fraction               = 0.3
    wing.append_control_surface(elevator)       

    RCAIDE.Methods.Geometry.Two_Dimensional.Planform.wing_planform(wing)     

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------ 
    wing                                  = RCAIDE.Components.Wings.Wing()
    wing.tag                              = 'vertical_stabilizer'   
    wing.sweeps.leading_edge              = 20 * Units.degrees 
    wing.thickness_to_chord               = 12.5
    wing.areas.reference                  = 1.163 
    wing.spans.projected                  = 1.4816
    wing.chords.root                      = 1.2176
    wing.chords.tip                       = 1.2176
    wing.chords.tip                       = 0.5870 
    wing.aspect_ratio                     = 1.8874 
    wing.taper                            = 0.4820 
    wing.chords.mean_aerodynamic          = 0.9390 
    wing.twists.root                      = 0 * Units.degrees  
    wing.twists.tip                       = 0 * Units.degrees   
    wing.origin                           = [[ 7.127369987, 0., 0.303750948]]
    wing.aerodynamic_center               = [ 7.49778005775, 0., 0.67416101875] 
    wing.vertical                         = True 
    wing.symmetric                        = False
    wing.t_tail                           = False
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0  
    
    rudder                                = RCAIDE.Components.Wings.Control_Surfaces.Rudder()
    rudder.tag                            = 'rudder'
    rudder.span_fraction_start            = 0.2
    rudder.span_fraction_end              = 0.8
    rudder.deflection                     = 0.0  * Units.deg
    rudder.chord_fraction                 = 0.2
    wing.append_control_surface(rudder) 
    
    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    fuselage = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage'
    fuselage.seats_abreast                      = 2
    fuselage.lengths.total                      = 8.349950916 
    fuselage.width                              = 1.22028016 
    fuselage.heights.maximum                    = 1.634415138  
    fuselage.areas.wetted                       = 12. # ESTIMATED 
    fuselage.areas.front_projected              = fuselage.width*fuselage.heights.maximum
    fuselage.effective_diameter                 = 1.22028016 

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0
    segment.percent_z_location                  = 0
    segment.height                              = 0.529255748
    segment.width                               = 0.575603849
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  =  0.028527593
    segment.percent_z_location                  =  0
    segment.height                              =  0.737072721
    segment.width                               =  0.921265952 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.187342754 
    segment.percent_z_location                  = 0 
    segment.height                              = 1.174231852 
    segment.width                               = 1.196956212
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.242034847 
    segment.percent_z_location                  = 0.011503528 
    segment.height                              = 1.450221906 
    segment.width                               = 1.173932059 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.296715183 
    segment.percent_z_location                  = 0.015984303 
    segment.height                              = 1.634415138 
    segment.width                               = 1.22028016 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.510275342 
    segment.percent_z_location                  = -0.005
    segment.height                              = 1.082135236 
    segment.width                               = 1.013062774 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.833284347 
    segment.percent_z_location                  = 0.014138855 
    segment.height                              = 0.621652157 
    segment.width                               = 0.414134978
    fuselage.Segments.append(segment)
 
    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 1
    segment.percent_z_location                  = 0.018978667 
    segment.height                              = 0.092096616 
    segment.width                               = 0.046048308 
    fuselage.Segments.append(segment)
    
    # add to vehicle
    vehicle.append_component(fuselage)
 
    # ------------------------------------------------------------------
    #   Fuel
    # ------------------------------------------------------------------    
    # define fuel weight needed to size fuel system
    fuel                                        = RCAIDE.Attributes.Propellants.Aviation_Gasoline()
    fuel.mass_properties                        = RCAIDE.Components.Mass_Properties() 
    fuel.number_of_tanks                        = 1.
    fuel.origin                                 = wing.origin
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density #all of the fuel volume is internal
    fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel.mass_properties.mass                   = 319 *Units.lbs
    vehicle.fuel                                = fuel




    # ########################################################  Energy Network  #########################################################  
    net                                         = RCAIDE.Energy.Networks.Internal_Combustion_Engine_Network()   

    # add the network to the vehicle
    vehicle.append_energy_network(net) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                   = RCAIDE.Energy.Networks.Distribution.Fuel_Line()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Fuel Tank & Fuel
    #------------------------------------------------------------------------------------------------------------------------------------       
    fuel_tank                                   = RCAIDE.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                            = wing.origin  
    fuel                                        = RCAIDE.Attributes.Propellants.Aviation_Gasoline() 
    fuel.mass_properties.mass                   = 319 *Units.lbs 
    fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    fuel_tank.fuel                              = fuel     
    fuel_line.fuel_tanks.append(fuel_tank)  
    net.fuel_lines.append(fuel_line)    

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    ice_prop    = RCAIDE.Energy.Propulsors.ICE_Propeller()     
    ice_prop.active_fuel_tanks                 = ['fuel_tank']   
                                                     
    # Engine                     
    engine                                     = RCAIDE.Energy.Propulsors.Converters.Engine()

    engine.sea_level_power                     = 185. * Units.horsepower 
    engine.rated_speed                         = 2300. * Units.rpm 
    engine.power_specific_fuel_consumption     = 0.01
    ice_prop.engine                            = engine 
     
    # Propeller 
    prop = RCAIDE.Energy.Propulsors.Converters.Propeller()
    prop.tag                                = 'propeller'
    prop.number_of_blades                   = 2.0
    prop.tip_radius                         = 76./2. * Units.inches
    prop.hub_radius                         = 8.     * Units.inches
    prop.cruise.design_freestream_velocity  = 119.   * Units.knots
    prop.cruise.design_angular_velocity     = 2650.  * Units.rpm
    prop.cruise.design_Cl                   = 0.8
    prop.cruise.design_altitude             = 12000. * Units.feet
    prop.cruise.design_power                = .64 * 180. * Units.horsepower
    prop.variable_pitch                     = True  
    ospath                                  = os.path.abspath(__file__)
    separator                               = os.path.sep
    rel_path                                = os.path.dirname(ospath) + separator 
    airfoil                                 = RCAIDE.Components.Airfoils.Airfoil()
    airfoil.tag                             = 'NACA_4412'   
    airfoil.coordinate_file                 =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                     =[ rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt']  
    prop.append_airfoil(airfoil)           
    prop.airfoil_polar_stations             = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  
    prop                                    = design_propeller(prop)    
    ice_prop.propeller                      = prop 
    
    fuel_line.propulsors.append(ice_prop)

    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Avionics
    #------------------------------------------------------------------------------------------------------------------------------------ 
    Wuav                                        = 2. * Units.lbs
    avionics                                    = RCAIDE.Energy.Peripherals.Avionics()
    avionics.mass_properties.uninstalled        = Wuav
    vehicle.avionics                            = avionics     

    #------------------------------------------------------------------------------------------------------------------------------------ 
    #   Vehicle Definition Complete
    #------------------------------------------------------------------------------------------------------------------------------------ 
     
    return vehicle


# ----------------------------------------------------------------------
#   Define the Configurations
# --------------------------------------------------------------------- 

def configs_setup(vehicle):
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs                                                    = RCAIDE.Components.Configs.Config.Container() 
    base_config                                                = RCAIDE.Components.Configs.Config(vehicle) 
    base_config.tag                                            = 'base'
    configs.append(base_config)
    
    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------ 
    config                                                     = RCAIDE.Components.Configs.Config(base_config)
    config.tag                                                 = 'cruise' 
    configs.append(config)
    
    
    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------ 
    config                                                     = RCAIDE.Components.Configs.Config(base_config)
    config.tag                                                 = 'takeoff' 
    config.wings['main_wing'].control_surfaces.flap.deflection = 20. * Units.deg
    config.V2_VS_ratio                                         = 1.21
    config.maximum_lift_coefficient                            = 2.
    
    configs.append(config)
    
    
    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config                                                     = RCAIDE.Components.Configs.Config(base_config)
    config.tag                                                 = 'landing' 
    config.wings['main_wing'].control_surfaces.flap.deflection = 20. * Units.deg
    config.Vref_VS_ratio                                       = 1.23
    config.maximum_lift_coefficient                            = 2.
                                                               
    configs.append(config) 
     
    return configs
# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results,run_stability):
    plot_flight_conditions(results) 
        
    # Plot Flight Conditions 
    plot_flight_conditions(results) 
    
    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results)  
    
    # Plot Aircraft Flight Speed
    plot_aircraft_velocities(results) 

    # Plot Aircraft Stability
    if run_stability:
        plot_stability_coefficients(results)  
    
    # Plot Propeller Conditions 
    plot_rotor_conditions(results) 
     
    # Plot Throttle
    plot_propulsor_throttles(results)   
      
    return
 
# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
 

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment
    base_segment = Segments.Segment() 

    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment     = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.base )   
    segment.altitude                                      = 12000. * Units.feet
    segment.air_speed                                     = 119.   * Units.knots
    segment.distance                                      = 10 * Units.nautical_mile
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    segment.flight_dynamics.moment_y                      = True 
    
    # define flight controls  
    segment.flight_controls.RPM.active                    = True           
    segment.flight_controls.RPM.assigned_propulsors       = [['ice_propeller']]
    segment.flight_controls.RPM.initial_guess             = True 
    segment.flight_controls.RPM.initial_guess_values      = [[2500]] 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['ice_propeller']]  
    segment.flight_controls.body_angle.active             = True   
    mission.append_segment(segment)


    return mission 

def missions_setup(mission): 
 
    missions         = RCAIDE.Analyses.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions  



if __name__ == '__main__': 
    main()    
    plt.show()