# Cessna_172.py
# 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import RCAIDE
from RCAIDE.Framework.Core                                     import Units  
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor        import design_propeller 
from RCAIDE.Library.Plots                                      import *    
 
import matplotlib.pyplot as plt
import numpy as np 
import os 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():   
    
    # vehicle data
    vehicle  = vehicle_setup() 
    
    # Set up vehicle configs
    configs  = configs_setup(vehicle)

    # create analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(analyses) 
    
    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 
     
    # mission analysis 
    results = missions.base_mission.evaluate()  
    
    plot_mission(results)


def vehicle_setup(): 
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties ########################################################  
    #------------------------------------------------------------------------------------------------------------------------------------     
    vehicle                                     = RCAIDE.Vehicle()
    vehicle.tag                                 = 'Cessna_172' 
    vehicle.mass_properties.max_takeoff         = 2550. * Units.pounds
    vehicle.mass_properties.takeoff             = 2550. * Units.pounds
    vehicle.mass_properties.max_zero_fuel       = 2550. * Units.pounds
    vehicle.mass_properties.cargo               = 0. 
                                               
    # envelope properties                       
    vehicle.flight_envelope.ultimate_load              = 5.7
    vehicle.flight_envelope.positive_limit_load                 = 3.8
                                                
    cruise_speed                                = 124. * Units.kts
    altitude                                    = 8500. * Units.ft
    atmo                                        = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    freestream                                  = atmo.compute_values (0.)
    freestream0                                 = atmo.compute_values (altitude)
    mach_number                                 = (cruise_speed/freestream.speed_of_sound)[0][0] 
    vehicle.design_dynamic_pressure             = ( .5 *freestream0.density*(cruise_speed*cruise_speed))[0][0]
    vehicle.flight_envelope.design_mach_number                  =  mach_number
                                                
    # basic parameters                          
    vehicle.reference_area                      = 174. * Units.feet**2       
    vehicle.passengers                          = 4


    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################################### Landing Gear ################################################################    
    #------------------------------------------------------------------------------------------------------------------------------------
    landing_gear                                = RCAIDE.Library.Components.Landing_Gear.Landing_Gear()
    main_gear                                   = RCAIDE.Library.Components.Landing_Gear.Main_Landing_Gear()
    nose_gear                                   = RCAIDE.Library.Components.Landing_Gear.Nose_Landing_Gear()
    main_gear.strut_length                      = 12. * Units.inches  
    nose_gear.strut_length                      = 6. * Units.inches 
                                                
    landing_gear.main                           = main_gear
    landing_gear.nose                           = nose_gear
                                                
    #add to vehicle                             
    vehicle.landing_gear                        = landing_gear 


    #------------------------------------------------------------------------------------------------------------------------------------
    # ######################################################## Wings ####################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------   

    wing                                        = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                                    = 'main_wing'    
    wing.sweeps.quarter_chord                   = 0.0 * Units.deg
    wing.thickness_to_chord                     = 0.12
    wing.areas.reference                        = 174. * Units.feet**2
    wing.spans.projected                        = 36.  * Units.feet + 1. * Units.inches
    wing.chords.root                            = 66. * Units.inches
    wing.chords.tip                             = 45. * Units.inches
    wing.chords.mean_aerodynamic                = 58. * Units.inches
    wing.taper                                  = wing.chords.tip/wing.chords.root
    wing.aspect_ratio                           = wing.spans.projected**2. / wing.areas.reference
    wing.twists.root                            = 3.0 * Units.degrees
    wing.twists.tip                             = 1.5 * Units.degrees
    wing.origin                                 = [[80.* Units.inches,0,0]]
    wing.aerodynamic_center                     = [22.* Units.inches,0,0]
    wing.vertical                               = False
    wing.symmetric                              = True
    wing.high_lift                              = True 
    wing.dynamic_pressure_ratio                 = 1.0 
                                          
    # control surfaces -------------------------------------------
    flap                                        = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap() 
    flap.tag                                    = 'flap' 
    flap.span_fraction_start                    = 0.15 
    flap.span_fraction_end                      = 0.324    
    flap.deflection                             = 1.0 * Units.deg
    flap.chord_fraction                         = 0.19    
    wing.append_control_surface(flap)           
                                                
    slat                                        = RCAIDE.Library.Components.Wings.Control_Surfaces.Slat() 
    slat.tag                                    = 'slat' 
    slat.span_fraction_start                    = 0.324 
    slat.span_fraction_end                      = 0.963     
    slat.deflection                             = 1.0 * Units.deg
    slat.chord_fraction                         = 0.1      
    wing.append_control_surface(slat)
    
    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.963
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.16
    wing.append_control_surface(aileron)
    
    RCAIDE.Library.Methods.Geometry.Planform.wing_planform(wing) 

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------        
                                                
    wing                                        = RCAIDE.Library.Components.Wings.Horizontal_Tail()
    wing.tag                                    = 'horizontal_stabilizer' 
    wing.sweeps.quarter_chord                   = 0.0 * Units.deg
    wing.thickness_to_chord                     = 0.12
    wing.areas.reference                        = 5800. * Units.inches**2
    wing.spans.projected                        = 136.  * Units.inches
    wing.chords.root                            = 55. * Units.inches
    wing.chords.tip                             = 30. * Units.inches
    wing.chords.mean_aerodynamic                = 43. * Units.inches 
    wing.taper                                  = wing.chords.tip/wing.chords.root
    wing.aspect_ratio                           = wing.spans.projected**2. / wing.areas.reference
    wing.twists.root                            = 0.0 * Units.degrees
    wing.twists.tip                             = 0.0 * Units.degrees
    wing.origin                                 = [[246.* Units.inches,0,0]]
    wing.aerodynamic_center                     = [20.* Units.inches,0,0]
    wing.vertical                               = False
    wing.symmetric                              = True
    wing.high_lift                              = False 
    wing.dynamic_pressure_ratio                 = 0.9
    vehicle.append_component(wing)
    
    # control surfaces -------------------------------------------
    elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.09
    elevator.span_fraction_end     = 0.92
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    wing.append_control_surface(elevator)    


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing                                        = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag                                    = 'vertical_stabilizer' 
    wing.sweeps.quarter_chord                   = 25. * Units.deg
    wing.thickness_to_chord                     = 0.12
    wing.areas.reference                        = 3500. * Units.inches**2
    wing.spans.projected                        = 73.   * Units.inches
    wing.chords.root                            = 66. * Units.inches
    wing.chords.tip                             = 27. * Units.inches
    wing.chords.mean_aerodynamic                = 48. * Units.inches 
    wing.taper                                  = wing.chords.tip/wing.chords.root
    wing.aspect_ratio                           = wing.spans.projected**2. / wing.areas.reference
    wing.twists.root                            = 0.0 * Units.degrees
    wing.twists.tip                             = 0.0 * Units.degrees
    wing.origin                                 = [[237.* Units.inches,0,0]]
    wing.aerodynamic_center                     = [20.* Units.inches,0,0] 
    wing.vertical                               = True 
    wing.symmetric                              = False
    wing.t_tail                                 = False 
    wing.dynamic_pressure_ratio                 = 1.0
    
    rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
    rudder.tag                   = 'rudder'
    rudder.span_fraction_start   = 0.1
    rudder.span_fraction_end     = 0.963
    rudder.deflection            = 0.0 * Units.degrees
    rudder.chord_fraction        = 0.3
    wing.append_control_surface(rudder)
    

    # add to vehicle
    vehicle.append_component(wing)


    #------------------------------------------------------------------------------------------------------------------------------------
    # ########################################################## Fuselage ############################################################### 
    #------------------------------------------------------------------------------------------------------------------------------------
    
    fuselage                                    = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.number_coach_seats                 = 4.        
    fuselage.differential_pressure              = 8*Units.psi                    # Maximum differential pressure
    fuselage.width                              = 42.         * Units.inches     # Width of the fuselage
    fuselage.heights.maximum                    = 62. * Units.inches    # Height of the fuselage
    fuselage.lengths.total                      = 326.         * Units.inches            # Length of the fuselage
    fuselage.lengths.empennage                  = 161. * Units.inches  
    fuselage.lengths.cabin                      = 105. * Units.inches
    fuselage.lengths.structure                  = fuselage.lengths.total-fuselage.lengths.empennage 
    fuselage.mass_properties.volume             = .4*fuselage.lengths.total*(np.pi/4.)*(fuselage.heights.maximum**2.) #try this as approximation
    fuselage.mass_properties.internal_volume    = .3*fuselage.lengths.total*(np.pi/4.)*(fuselage.heights.maximum**2.)
    fuselage.areas.wetted                       = 30000. * Units.inches**2.
    fuselage.seats_abreast                      = 2.
    fuselage.fineness.nose                      = 1.6
    fuselage.fineness.tail                      = 2.
    fuselage.lengths.nose                       = 60.  * Units.inches
    fuselage.heights.at_quarter_length          = 62. * Units.inches
    fuselage.heights.at_three_quarters_length   = 62. * Units.inches
    fuselage.heights.at_wing_root_quarter_chord = 23. * Units.inches
    fuselage.areas.front_projected              = fuselage.width* fuselage.heights.maximum
    fuselage.effective_diameter                 = 50. * Units.inches

    # add to vehicle
    vehicle.append_component(fuselage)
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ########################################################## Energy Network ######################################################### 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    #initialize the fuel network
    net                                         = RCAIDE.Framework.Networks.Fuel()   

    # add the network to the vehicle
    vehicle.append_energy_network(net) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                   = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Fuel Tank & Fuel
    #------------------------------------------------------------------------------------------------------------------------------------       
    fuel_tank                                   = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                            = wing.origin  
    fuel                                        = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline() 
    fuel.mass_properties.mass                   = 319 *Units.lbs 
    fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    fuel_tank.fuel                              = fuel     
    fuel_line.fuel_tanks.append(fuel_tank)  
    net.fuel_lines.append(fuel_line)    

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    ice_prop    = RCAIDE.Library.Components.Propulsors.ICE_Propeller()     
    ice_prop.active_fuel_tanks                 = ['fuel_tank']   
                                                     
    # Engine                     
    engine                                     = RCAIDE.Library.Components.Propulsors.Converters.Engine()
    engine.sea_level_power                     = 180. * Units.horsepower
    engine.flat_rate_altitude                  = 0.0
    engine.rated_speed                         = 2700. * Units.rpm
    engine.power_specific_fuel_consumption     = 0.52 
    ice_prop.engine                            = engine 
     
    # Propeller 
    prop = RCAIDE.Library.Components.Propulsors.Converters.Propeller()
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
    airfoil                                 = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                             = 'NACA_4412' 
    airfoil.coordinate_file                 =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                     =[ rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                               rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt']  
    prop.append_airfoil(airfoil)      
    prop.airfoil_polar_stations             = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  
    design_propeller(prop)    
    ice_prop.propeller                      = prop 
    
    fuel_line.propulsors.append(ice_prop)

    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Avionics
    #------------------------------------------------------------------------------------------------------------------------------------ 
    Wuav                                        = 2. * Units.lbs
    avionics                                    = RCAIDE.Library.Components.Systems.Avionics()
    avionics.mass_properties.uninstalled        = Wuav
    vehicle.avionics                            = avionics     

    #------------------------------------------------------------------------------------------------------------------------------------ 
    #   Vehicle Definition Complete
    #------------------------------------------------------------------------------------------------------------------------------------ 

    return vehicle
  
  
def configs_setup(vehicle):
     # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs                                                    = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config                                                = RCAIDE.Library.Components.Configs.Config(vehicle) 
    base_config.tag                                            = 'base'
    configs.append(base_config)
    
    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------ 
    config                                                     = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag                                                 = 'cruise' 
    configs.append(config)
    
    
    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------ 
    config                                                     = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag                                                 = 'takeoff' 
    config.wings['main_wing'].control_surfaces.flap.deflection = 20. * Units.deg
    config.V2_VS_ratio                                         = 1.21
    config.maximum_lift_coefficient                            = 2.
    
    configs.append(config)
    
    
    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config                                                     = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag                                                 = 'landing' 
    config.wings['main_wing'].control_surfaces.flap.deflection = 20. * Units.deg
    config.Vref_VS_ratio                                       = 1.23
    config.maximum_lift_coefficient                            = 2.
                                                               
    configs.append(config) 
    
    # done!
    return configs
 
# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
 

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    # base segment
    base_segment = Segments.Segment()
    


    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment     = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.base ) 
    segment.altitude                                = 12000. * Units.feet
    segment.air_speed                               = 119.   * Units.knots
    segment.distance                                = 10 * Units.nautical_mile  
    segment.state.numerics.number_control_points    = 4  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls  
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['ice_propeller']] 
    segment.assigned_control_variables.body_angle.active             = True                   
    
    
    mission.append_segment(segment)


    return mission

def analyses_setup(configs):
    """Set up analyses for each of the different configurations."""

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # Build a base analysis for each configuration. Here the base analysis is always used, but
    # this can be modified if desired for other cases.
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):
    """This is the baseline set of analyses to be used with this vehicle. Of these, the most
    commonly changed are the weights and aerodynamics methods."""

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Athena_Vortex_Lattice()
    aerodynamics.settings.filenames.avl_bin_name   =  '/Users/aidanmolloy/Documents/LEADS/Codes/AVL/avl3.35'
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.number_of_spanwise_vortices   = 60
    aerodynamics.settings.number_of_chordwise_vortices  = 5   
    analyses.append(aerodynamics)
 
    # ------------------------------------------------------------------
    #  Energy
    energy = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle
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

    return analyses    
    
    

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):
    """This function defines the baseline mission that will be flown by the aircraft in order
    to compute performance."""

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'the_mission'
  
    Segments = RCAIDE.Framework.Mission.Segments 
    base_segment = Segments.Segment()

    # ------------------------------------------------------------------------------------------------------------------------------------ 
    #   Takeoff Roll
    # ------------------------------------------------------------------------------------------------------------------------------------ 

    #segment = Segments.Ground.Takeoff(base_segment)
    #segment.tag = "Takeoff" 
    #segment.analyses.extend( analyses.takeoff )
    #segment.velocity_start           = 10.* Units.knots
    #segment.velocity_end             = 125.0 * Units['m/s']
    #segment.friction_coefficient     = 0.04
    #segment.altitude                 = 0.0   
    #mission.append_segment(segment)

    ## ------------------------------------------------------------------
    ##   First Climb Segment: Constant Speed Constant Rate  
    ## ------------------------------------------------------------------

    #segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag = "climb_1" 
    #segment.analyses.extend( analyses.takeoff ) 
    #segment.altitude_start = 0.0   * Units.km
    #segment.altitude_end   = 3.0   * Units.km
    #segment.air_speed      = 125.0 * Units['m/s']
    #segment.climb_rate     = 6.0   * Units['m/s']  
     
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                      = True  
    #segment.flight_dynamics.force_z                      = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                 
    
    #mission.append_segment(segment)


    ## ------------------------------------------------------------------
    ##   Second Climb Segment: Constant Speed Constant Rate  
    ## ------------------------------------------------------------------    

    #segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag = "climb_2" 
    #segment.analyses.extend( analyses.cruise ) 
    #segment.altitude_end   = 8.0   * Units.km
    #segment.air_speed      = 190.0 * Units['m/s']
    #segment.climb_rate     = 6.0   * Units['m/s']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                      = True  
    #segment.flight_dynamics.force_z                      = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                  
    
    #mission.append_segment(segment)


    ## ------------------------------------------------------------------
    ##   Third Climb Segment: Constant Speed Constant Rate  
    ## ------------------------------------------------------------------    

    #segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag = "climb_3" 
    #segment.analyses.extend( analyses.cruise ) 
    #segment.altitude_end = 10.5   * Units.km
    #segment.air_speed    = 226.0  * Units['m/s']
    #segment.climb_rate   = 3.0    * Units['m/s']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                      = True  
    #segment.flight_dynamics.force_z                      = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)


    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise" 
    segment.analyses.extend( analyses.cruise ) 
    segment.altitude                                      = 1. * Units.km  
    segment.air_speed                                     = 60 * Units['m/s']
    segment.distance                                      = 20 * Units.nmi   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['ice_propeller']] 
    segment.assigned_control_variables.body_angle.active             = True                
    
    mission.append_segment(segment)


    ## ------------------------------------------------------------------
    ##   First Descent Segment: Constant Speed Constant Rate  
    ## ------------------------------------------------------------------

    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag = "descent_1" 
    #segment.analyses.extend( analyses.cruise ) 
    #segment.altitude_start                                = 10.5 * Units.km 
    #segment.altitude_end                                  = 8.0   * Units.km
    #segment.air_speed                                     = 220.0 * Units['m/s']
    #segment.descent_rate                                  = 4.5   * Units['m/s']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)


    ## ------------------------------------------------------------------
    ##   Second Descent Segment: Constant Speed Constant Rate  
    ## ------------------------------------------------------------------

    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag  = "descent_2" 
    #segment.analyses.extend( analyses.landing ) 
    #segment.altitude_end                                  = 6.0   * Units.km
    #segment.air_speed                                     = 195.0 * Units['m/s']
    #segment.descent_rate                                  = 5.0   * Units['m/s']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)


    ## ------------------------------------------------------------------
    ##   Third Descent Segment: Constant Speed Constant Rate  
    ## ------------------------------------------------------------------

    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag = "descent_3"  
    #segment.analyses.extend( analyses.landing ) 
    #segment.altitude_end                                  = 4.0   * Units.km
    #segment.air_speed                                     = 170.0 * Units['m/s']
    #segment.descent_rate                                  = 5.0   * Units['m/s']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)


    ## ------------------------------------------------------------------
    ##   Fourth Descent Segment: Constant Speed Constant Rate  
    ## ------------------------------------------------------------------

    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag = "descent_4" 
    #segment.analyses.extend( analyses.landing ) 
    #segment.altitude_end                                  = 2.0   * Units.km
    #segment.air_speed                                     = 150.0 * Units['m/s']
    #segment.descent_rate                                  = 5.0   * Units['m/s']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)



    ## ------------------------------------------------------------------
    ##   Fifth Descent Segment:Constant Speed Constant Rate  
    ## ------------------------------------------------------------------

    #segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    #segment.tag = "descent_5" 
    #segment.analyses.extend( analyses.landing ) 
    #segment.altitude_end                                  = 0.0   * Units.km
    #segment.air_speed                                     = 145.0 * Units['m/s']
    #segment.descent_rate                                  = 3.0   * Units['m/s']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['outer_starboard_propulsor','inner_starboard_propulsor','outer_port_propulsor','inner_port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)
    ## ------------------------------------------------------------------------------------------------------------------------------------ 
    ##   Landing Roll
    ## ------------------------------------------------------------------------------------------------------------------------------------ 

    #segment = Segments.Ground.Landing(base_segment)
    #segment.tag = "Landing"

    #segment.analyses.extend( analyses.landing ) 
    #segment.velocity_end                                  = 10 * Units.knots 
    #segment.friction_coefficient                          = 0.4
    #segment.altitude                                      = 0.0   
    #segment.assigned_control_variables.elapsed_time.active           = True  
    #segment.assigned_control_variables.elapsed_time.initial_guess_values   = [[30.]]  
    #mission.append_segment(segment)     


    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------

    return mission

def missions_setup(mission):
    """This allows multiple missions to be incorporated if desired, but only one is used here."""

    missions     = RCAIDE.Framework.Mission.Missions() 
    mission.tag  = 'base_mission'
    missions.append(mission)

    return missions

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------
def plot_mission(results):
    """This function plots the results of the mission analysis and saves those results to 
    png files."""

    # Plot Flight Conditions 
    plot_flight_conditions(results)
    
    # Plot Aerodynamic Forces 
    plot_aerodynamic_forces(results)
    
    # Plot Aerodynamic Coefficients 
    plot_aerodynamic_coefficients(results)     
    
    # Drag Components
    plot_drag_components(results)
    
    # Plot Altitude, sfc, vehicle weight 
    plot_altitude_sfc_weight(results)
    
    # Plot Velocities 
    plot_aircraft_velocities(results)  
        
    return

if __name__ == '__main__': 
    main()
    plt.show()
