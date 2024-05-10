# Navion.py
# 
# Created: Dec 2021, M. Clarke

""" setup file for a mission with a Navion
"""
# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
import RCAIDE 
from RCAIDE.Framework.Core import Units   
from RCAIDE.Library.Methods.Energy.Propulsors.Converters.Rotor import design_propeller
from RCAIDE.Library.Methods.Geometry.Planform  import segment_properties
from RCAIDE.Library.Plots       import *  

# python imports 
import os 
import numpy as np
import pylab as plt
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    run_stability     =  True 
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

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config, stability_method,run_stability, configs)
        analyses[tag] = analysis

    return analyses


def base_analysis(vehicle,stability_method, run_stability,configs):

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
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Subsonic_VLM() 
    aerodynamics.geometry                            = vehicle
    aerodynamics.settings.number_spanwise_vortices   = 30
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics) 
    
    if run_stability:
        # ------------------------------------------------------------------
        #  Stability Analysis
        if stability_method == "avl":
            stability                                   = RCAIDE.Framework.Analyses.Stability.AVL()     
            stability.settings.filenames.avl_bin_name   = 'C:\\Users\\Matteo\\Documents\\UIUC\\avl.exe' #/Users/matthewclarke/Documents/AVL/avl3.35'     
            stability.settings.print_output             = False 
        elif stability_method == "vlm":
            stability                                   = RCAIDE.Framework.Analyses.Stability.VLM_Perturbation_Method() 
        #elif stability_method == "analytical":
            #stability                                   = RCAIDE.Framework.Analyses.Stability.Analytical_Approximation() 
    
        stability.configuration                         = configs
        stability.geometry                              = vehicle
        analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Framework.Analyses.Energy.Energy()
    energy.networks = vehicle.networks  
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

    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.sweeps.quarter_chord             = 0 * Units.degrees 
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 17.112
    wing.chords.mean_aerodynamic          = 1.74 
    wing.taper                            = 1 
    wing.aspect_ratio                     = 6.04  
    wing.spans.projected                  = 10.166
    wing.chords.root                      = 2.1944 
    wing.chords.tip                       = 2.1944
    wing.twists.root                      = 0 * Units.degrees  
    wing.twists.tip                       = 0 * Units.degrees   
    wing.dihedral                         = 0 * Units.degrees   
    wing.origin                           = [[1.652555594, 0.,-0.6006666]]
    wing.aerodynamic_center               = [1.852555594, 0., 6006666 ] # INCORRECT 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0   
    
    tip_airfoil                           = RCAIDE.Library.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = 'Airfoils/NACA_6410.txt' 
 
    root_airfoil                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    root_airfoil.coordinate_file          = 'Airfoils/NACA_6410.txt' 
    
    # add to vehicle
    vehicle.append_component(wing) 
 
    # ------------------------------------------------------------------
    #   Fuel
    # ------------------------------------------------------------------    
    # define fuel weight needed to size fuel system
    fuel                                        = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline()
    fuel.mass_properties                        = RCAIDE.Library.Components.Mass_Properties() 
    fuel.number_of_tanks                        = 1.
    fuel.origin                                 = wing.origin
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density #all of the fuel volume is internal
    fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel.mass_properties.mass                   = 319 *Units.lbs
    vehicle.fuel                                = fuel

    # ########################################################  Energy Network  #########################################################  
    net                                         = RCAIDE.Framework.Networks.Internal_Combustion_Engine_Network()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                   = RCAIDE.Library.Components.Energy.Distribution.Fuel_Line()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Fuel Tank & Fuel
    #------------------------------------------------------------------------------------------------------------------------------------       
    fuel_tank                                   = RCAIDE.Library.Components.Energy.Fuel_Tanks.Fuel_Tank()
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
    ice_prop                                   = RCAIDE.Library.Components.Propulsors.ICE_Propeller()     
    ice_prop.active_fuel_tanks                 = ['fuel_tank']   
                                                     
    # Engine                     
    engine                                     = RCAIDE.Library.Components.Propulsors.Converters.Engine()

    engine.sea_level_power                     = 185. * Units.horsepower 
    engine.rated_speed                         = 2300. * Units.rpm 
    engine.power_specific_fuel_consumption     = 0.01  * Units['lb/hp/hr']
    ice_prop.engine                            = engine 
     
    # Propeller 
    prop                                    = RCAIDE.Library.Components.Propulsors.Converters.Propeller()
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
    airfoil                                 = RCAIDE.Library.Components.Airfoils.NACA_4_Series_Airfoil()
    airfoil.NACA_4_Series_code              = '4412'   
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

    # add the network to the vehicle
    vehicle.append_energy_network(net) 

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


# ----------------------------------------------------------------------
#   Define the Configurations
# --------------------------------------------------------------------- 

def configs_setup(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs                                                    = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config                                                = RCAIDE.Library.Components.Configs.Config(vehicle) 
    base_config.tag                                            = 'base'
    configs.append(base_config)

    return configs
# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results,run_stability):
    
    # Plot Flight Conditions 
    plot_flight_conditions(results) 
    
    # Plot Aerodynamic Coefficients
    plot_aerodynamic_coefficients(results)  
    
    # Plot Aircraft Flight Speed
    plot_aircraft_velocities(results) 

    # Plot Aircraft Stability
    if run_stability:
        plot_longitudinal_stability(results)  
        
        plot_lateral_stability(results) 
        
        plot_flight_forces_and_moments(results)
    
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
    segment.tag = "base" 
    segment.analyses.extend( analyses.base )   
    segment.altitude                                                 = 12000. * Units.feet
    segment.air_speed                                                = 120 * Units['mph']
    segment.distance                                                 = 10 * Units.nautical_mile 
    segment.sideslip_angle                                           = 0 * Units.degrees 
                
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True    
    segment.flight_dynamics.force_z                       = True    
                
    # define flight controls              
    #segment.flight_controls.RPM.active                               = True           
    #segment.flight_controls.RPM.assigned_propulsors                  = [['ice_propeller']]
    #segment.flight_controls.RPM.initial_guess                        = True 
    #segment.flight_controls.RPM.initial_guess_values                 = [[2500]] 
    #segment.flight_controls.throttle.active                          = True           
    #segment.flight_controls.throttle.assigned_propulsors             = [['ice_propeller']]    
    #segment.flight_controls.body_angle.active                        = True   
    
    ## Longidinal Flight Mechanics
    segment.flight_dynamics.moment_y                                 = True 
    #segment.flight_controls.elevator_deflection.active               = True    
    #segment.flight_controls.elevator_deflection.assigned_surfaces    = [['elevator']]
    #segment.flight_controls.elevator_deflection.initial_guess_values = [[0]]     
   
    ## Lateral Flight Mechanics 
    segment.flight_dynamics.force_y                       = True     
    segment.flight_dynamics.moment_x                                 = True
    segment.flight_dynamics.moment_z                                 = True 
    #segment.flight_controls.aileron_deflection.active               = True    
    #segment.flight_controls.aileron_deflection.assigned_surfaces    = [['aileron']]
    #segment.flight_controls.aileron_deflection.initial_guess_values = [[0]] 
    #segment.flight_controls.rudder_deflection.active               = True    
    #segment.flight_controls.rudder_deflection.assigned_surfaces    = [['rudder']]
    #segment.flight_controls.rudder_deflection.initial_guess_values = [[0]]
    #segment.flight_controls.bank_angle.active                      = True    
    #segment.flight_controls.bank_angle.initial_guess_values        = [[0]]     
    
    mission.append_segment(segment)

    return mission 

def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions  



if __name__ == '__main__': 
    main()    
    plt.show()