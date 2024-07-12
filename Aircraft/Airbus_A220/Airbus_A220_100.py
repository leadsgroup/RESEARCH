# Airbus A220-100
# 
# 
# Created:  Mar 2024, M. Guidotti, S. Shekar 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 

# RCAIDE imports 
import RCAIDE
from RCAIDE.Core                                           import Units       
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform      import segment_properties  , segment_properties 
from RCAIDE.Methods.Energy.Propulsors.Turbofan_Propulsor   import design_turbofan , compute_nacelle_geometry
from RCAIDE.Methods.Stability.Center_of_Gravity            import compute_component_centers_of_gravity
from RCAIDE.Visualization                                  import *     

# Python imports 
import numpy                                               as np
import matplotlib.pyplot                                   as plt
from copy                                                  import deepcopy 
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
    mission = mission_setup(analyses)

    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 

    # mission analysis 
    results = missions.base_mission.evaluate()  

    # plot the results
    plot_mission(results) 

    # plot vehicle 
    plot_3d_vehicle(configs.base,
                    show_wing_control_points               = False,
                    show_rotor_wake_vortex_core            = False,
                    min_x_axis_limit                       = 0,
                    max_x_axis_limit                       = 40,
                    min_y_axis_limit                       = -20,
                    max_y_axis_limit                       = 20,
                    min_z_axis_limit                       = -20,
                    max_z_axis_limit                       = 20)         

    return


def analyses_setup(configs):
    """Set up analyses for each of the different configurations."""

    analyses = RCAIDE.Analyses.Analysis.Container()

    # Build a base analysis for each configuration. Here the base analysis is always used, but
    # this can be modified if desired for other cases.
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):
    """This is the baseline set of analyses to be used with this vehicle. Of these, the most
    commonly changed are the weights and aerodynamics methods."""


    return 
def vehicle_setup(): 

    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    

    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Airbus_220-100'   


    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------

    # mass properties
    vehicle.mass_properties.max_takeoff               = 63100 
    vehicle.mass_properties.takeoff                   = 63100  
    vehicle.mass_properties.operating_empty           = 56100 # Not sure  
    vehicle.mass_properties.max_zero_fuel             = 52200 
    vehicle.mass_properties.cargo                     = 15100 
    vehicle.mass_properties.center_of_gravity         = [[0,0,0]] # Unknown 
    vehicle.mass_properties.moments_of_inertia.tensor = [[0,0,0]] # Unknown 
    vehicle.mass_properties.max_fuel                  = 5000
    vehicle.design_mach_number                        = 0.82 
    vehicle.design_range                              = 63960000*Units.meter  
    vehicle.design_cruise_alt                         = 35000 *Units.feet

    # envelope properties
    vehicle.envelope.ultimate_load        = 3.75
    vehicle.envelope.limit_load           = 1.5

    # basic parameters       
    vehicle.reference_area                = 112.3* Units['meters**2']
    vehicle.passengers                    = 135
    vehicle.systems.control               = "fully powered"
    vehicle.systems.accessories           = "medium range"  




    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------


    # 737
    wing                                  = RCAIDE.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.aspect_ratio                     = 10.97
    wing.spand.projected                  = 35.1
    wing.sweeps.quarter_chord             = 25 * Units.deg
    wing.thickness_to_chord               = 0.1
    wing.taper                            = 0.1 
    wing.spans.projected                  = 34.32 
    wing.chords.root                      = 7.760 * Units.meter
    wing.chords.tip                       = 0.782 * Units.meter
    wing.chords.mean_aerodynamic          = 4.235 * Units.meter 
    wing.areas.reference                  = 112.3
    wing.areas.exposed                    = 2 * wing.areas.reference
    wing.areas.wetted                     = 2 * wing.areas.reference 
    wing.twists.root                      = 4.0 * Units.degrees
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[13.61,0,-0.5]]
    wing.aerodynamic_center               = [0,0,0] 
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.dynamic_pressure_ratio           = 1.0

    # atr 
    wing                                  = RCAIDE.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing'
    wing.areas.reference                  = 112.3  
    wing.spans.projected                  = 27.05 
    wing.aspect_ratio                     = (wing.spans.projected**2) /  wing.areas.reference
    wing.sweeps.quarter_chord             = 0.0 
    wing.thickness_to_chord               = 0.1 
    wing.chords.root                      = 2.7 
    wing.chords.tip                       = 1.35 
    wing.total_length                     = wing.chords.root  
    wing.taper                            = wing.chords.tip/wing.chords.root 
    wing.chords.mean_aerodynamic          = wing.chords.root * 2/3 * (( 1 + wing.taper + wing.taper**2 )/( 1 + wing.taper )) 
    wing.areas.exposed                    = 2 * wing.areas.reference
    wing.areas.wetted                     = 2 * wing.areas.reference 
    wing.twists.root                      = 0 * Units.degrees  
    wing.twists.tip                       = 0 * Units.degrees   
    wing.origin                           = [[11.52756129,0,2.009316366]]  
    wing.aerodynamic_center               = [[11.52756129 + 0.25*wing.chords.root ,0,2.009316366]]   
    wing.vertical                         = False   
    wing.symmetric                        = True  
    wing.dynamic_pressure_ratio           = 1.0 


    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

    wing     = RCAIDE.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'

    wing.aspect_ratio            = 
    wing.sweeps.quarter_chord    = 
    wing.thickness_to_chord      = 
    wing.taper                   = 
    wing.spans.projected         = 
    wing.chords.root             = 
    wing.chords.tip              = 
    wing.chords.mean_aerodynamic = 
    wing.areas.reference         = 
    wing.areas.exposed           = 
    wing.areas.wetted            = 
    wing.twists.root             = 
    wing.twists.tip              = 
    wing.origin                  = 
    wing.aerodynamic_center      = 
    wing.vertical                = 
    wing.symmetric               = 
    wing.dynamic_pressure_ratio  = 





    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing = RCAIDE.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'

    wing.aspect_ratio            = 
    wing.sweeps.quarter_chord    = 
    wing.thickness_to_chord      = 
    wing.taper                   = 

    wing.spans.projected         = 
    wing.total_length            = 

    wing.chords.root             = 
    wing.chords.tip              = 
    wing.chords.mean_aerodynamic = 

    wing.areas.reference         = 
    wing.areas.wetted            = 

    wing.twists.root             = 
    wing.twists.tip              = 

    wing.origin                  = 
    wing.aerodynamic_center      = 

    wing.vertical                = 
    wing.symmetric               = 
    wing.t_tail                  = 

    wing.dynamic_pressure_ratio  = 



    # add to vehicle
    vehicle.append_component(wing)

    # ################################################# Fuselage ################################################################ 

    fuselage                                    = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage' 
    fuselage.number_coach_seats                 = vehicle.passengers 
    fuselage.seats_abreast                      = 4
    fuselage.seat_pitch                         = 32.0    *Units.inches
    fuselage.fineness.nose                      = 1.1
    fuselage.fineness.tail                      = 2.1
    fuselage.lengths.nose                       = 3.8     *Units.meter
    fuselage.lengths.tail                       = 7.4     *Units.meter
    fuselage.lengths.total                      = 34.9    *Units.meter
    # fuselage.lengths.fore_space                 = #?
    # fuselage.lengths.aft_space                  = #?
    fuselage.lengths.cabin                      = fuselage.lengths.total- (fuselage.lengths.nose + fuselage.lengths.tail  )
    fuselage.width                              = 3.5
    fuselage.heights.maximum                    = 3.7
    fuselage.effective_diameter                 = 3.7
    fuselage.areas.side_projected               = 1.0 # incorrect 
    fuselage.areas.wetted                       = 1.0 # incorrect 
    fuselage.areas.front_projected              = 1.0 # incorrect 
    fuselage.differential_pressure              = 58.6e3 *Units.pascal
    #fuselage.heights.at_quarter_length          = 
    #fuselage.heights.at_three_quarters_length   = 
    #fuselage.heights.at_wing_root_quarter_chord = 

    # Segment  
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment() 
    segment.tag                                 = 'segment_1'    
    segment.percent_x_location                  = 0.0000
    segment.percent_z_location                  = 0.0000
    segment.height                              = 1E-3
    segment.width                               = 1E-3  
    fuselage.Segments.append(segment)       

    # add to vehicle
    vehicle.append_component(fuselage)


    # ################################################# Energy Network #######################################################         
    # Step 1: Define network
    # Step 2: Define Distribution Type
    # Step 3: Define Propulsors 
    # Step 4: Define Enegy Source 

    #------------------------------------------------------------------------------------------------------------------------- 
    #  Turbofan Network
    #-------------------------------------------------------------------------------------------------------------------------   
    net                                         = RCAIDE.Energy.Networks.Turbofan_Engine_Network() 

    # Append energy network to aircraft 
    vehicle.append_energy_network(net)   

    #------------------------------------------------------------------------------------------------------------------------- 
    # Fuel Distrubition Line 
    #------------------------------------------------------------------------------------------------------------------------- 
    fuel_line                                   = RCAIDE.Energy.Networks.Distribution.Fuel_Line()  

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------         
    turbofan                                    = RCAIDE.Energy.Propulsors.Turbofan() 
    turbofan.tag                                = 'starboard_propulsor'
    turbofan.active_fuel_tanks                  = ['fuel_tank']   
    turbofan.origin                             = [[13.72, 4.86,-1.1]] 
    turbofan.engine_length                      = 2.71     
    turbofan.bypass_ratio                       = 5.4    
    turbofan.design_altitude                    = 35000.0*Units.ft
    turbofan.design_mach_number                 = 0.78   
    turbofan.design_thrust                      = 35000.0* Units.N 

    # fan                
    fan                                         = RCAIDE.Energy.Propulsors.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.7   
    turbofan.fan                                = fan        

    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Attributes.Gases.Air() 
    ram                                         = RCAIDE.Energy.Propulsors.Converters.Ram()
    ram.tag                                     = 'ram' 
    turbofan.ram                                = ram 

    # inlet nozzle          
    inlet_nozzle                                = RCAIDE.Energy.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                            = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency          = 0.98
    inlet_nozzle.pressure_ratio                 = 0.98 
    turbofan.inlet_nozzle                       = inlet_nozzle 

    # low pressure compressor    
    low_pressure_compressor                       = RCAIDE.Energy.Propulsors.Converters.Compressor()    
    low_pressure_compressor.tag                   = 'lpc'
    low_pressure_compressor.polytropic_efficiency = 0.91
    low_pressure_compressor.pressure_ratio        = 1.9   
    turbofan.low_pressure_compressor              = low_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Energy.Propulsors.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 10.0    
    turbofan.high_pressure_compressor              = high_pressure_compressor

    # low pressure turbine  
    low_pressure_turbine                           = RCAIDE.Energy.Propulsors.Converters.Turbine()   
    low_pressure_turbine.tag                       ='lpt'
    low_pressure_turbine.mechanical_efficiency     = 0.99
    low_pressure_turbine.polytropic_efficiency     = 0.93 
    turbofan.low_pressure_turbine                  = low_pressure_turbine

    # high pressure turbine     
    high_pressure_turbine                          = RCAIDE.Energy.Propulsors.Converters.Turbine()   
    high_pressure_turbine.tag                      ='hpt'
    high_pressure_turbine.mechanical_efficiency    = 0.99
    high_pressure_turbine.polytropic_efficiency    = 0.93 
    turbofan.high_pressure_turbine                 = high_pressure_turbine 

    # combustor  
    combustor                                      = RCAIDE.Energy.Propulsors.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1500
    combustor.pressure_ratio                       = 0.95
    combustor.fuel_data                            = RCAIDE.Attributes.Propellants.Jet_A()  
    turbofan.combustor                             = combustor

    # core nozzle
    core_nozzle                                    = RCAIDE.Energy.Propulsors.Converters.Expansion_Nozzle()   
    core_nozzle.tag                                = 'core nozzle'
    core_nozzle.polytropic_efficiency              = 0.95
    core_nozzle.pressure_ratio                     = 0.99  
    turbofan.core_nozzle                           = core_nozzle

    # fan nozzle             
    fan_nozzle                                     = RCAIDE.Energy.Propulsors.Converters.Expansion_Nozzle()   
    fan_nozzle.tag                                 = 'fan nozzle'
    fan_nozzle.polytropic_efficiency               = 0.95
    fan_nozzle.pressure_ratio                      = 0.99 
    turbofan.fan_nozzle                            = fan_nozzle 

    # design turbofan
    design_turbofan(turbofan)  
    # append propulsor to distribution line 



    # Nacelle 
    nacelle                                     = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.diameter                            = 2.05
    nacelle.length                              = 2.71
    nacelle.tag                                 = 'nacelle_1'
    nacelle.inlet_diameter                      = 2.0
    nacelle.origin                              = [[13.5,4.38,-1.5]] 
    nacelle.areas.wetted                        = 1.1*np.pi*nacelle.diameter*nacelle.length
    nacelle.Airfoil.NACA_4_series_flag          = True 
    nacelle.Airfoil.coordinate_file             = '2410' 


    nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                             = 'segment_1'
    nac_segment.percent_x_location              = 0.0  
    nac_segment.height                          = 2.05
    nac_segment.width                           = 2.05
    nacelle.append_segment(nac_segment)         

    nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                             = 'segment_2'
    nac_segment.percent_x_location              = 0.3
    nac_segment.height                          = 2.1  
    nac_segment.width                           = 2.1 
    nacelle.append_segment(nac_segment)         

    nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                             = 'segment_3'
    nac_segment.percent_x_location              = 0.4  
    nac_segment.height                          = 2.05
    nac_segment.width                           = 2.05 
    nacelle.append_segment(nac_segment)         

    nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                             = 'segment_4'
    nac_segment.percent_x_location              = 0.75  
    nac_segment.height                          = 1.9
    nac_segment.width                           = 1.9
    nacelle.append_segment(nac_segment)         

    nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                             = 'segment_5'
    nac_segment.percent_x_location              = 1.0
    nac_segment.height                          = 1.7 
    nac_segment.width                           = 1.7
    nacelle.append_segment(nac_segment)           
    compute_nacelle_geometry(turbofan, nacelle)
    turbofan.nacelle                            = nacelle

    fuel_line.propulsors.append(turbofan)  

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Port Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------      
    # copy turbofan
    turbofan_2                                  = deepcopy(turbofan)
    turbofan_2.active_fuel_tanks                = ['fuel_tank'] 
    turbofan_2.tag                              = 'port_propulsor' 
    turbofan_2.origin                           = [[13.72,-4.38,-1.1]]  # change origin 
    turbofan_2.nacelle.origin                   = [[13.5,-4.38,-1.5]]

    # append propulsor to distribution line 
    fuel_line.propulsors.append(turbofan_2)

    #------------------------------------------------------------------------------------------------------------------------- 
    #  Energy Source: Fuel Tank
    #------------------------------------------------------------------------------------------------------------------------- 
    # fuel tank
    fuel_tank                                   = RCAIDE.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                            = wing.origin 

    # append fuel 
    fuel                                        = RCAIDE.Attributes.Propellants.Aviation_Gasoline()   
    fuel.mass_properties.mass                   = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    fuel.origin                                 = vehicle.wings.main_wing.mass_properties.center_of_gravity      
    fuel.mass_properties.center_of_gravity      = vehicle.wings.main_wing.aerodynamic_center
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    fuel_tank.fuel                              = fuel            

    # apend fuel tank to dataclass of fuel tanks on fuel line 
    fuel_line.fuel_tanks.append(fuel_tank) 

    # Append fuel line to Network      
    net.fuel_lines.append(fuel_line)  

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


    return   

#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):
    """This function defines the baseline mission that will be flown by the aircraft in order
    to compute performance."""



    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------

    return mission

def missions_setup(mission):
    """This allows multiple missions to be incorporated if desired, but only one is used here."""


    return   

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

    # Plot Static Stability Coefficients 
    plot_stability_coefficients(results)    

    # Drag Components
    plot_drag_components(results)

    # Plot Altitude, sfc, vehicle weight 
    plot_altitude_sfc_weight(results)

    # Plot Velocities 
    plot_aircraft_velocities(results)  

    return

# This section is needed to actually run the various functions in the file
if __name__ == '__main__': 
    main()
    plt.show()