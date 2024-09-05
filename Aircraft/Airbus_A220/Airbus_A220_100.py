# Airbus A220-100
# 
# 
# Created:  Mar 2024, S. Shekar 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 

# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core                                           import Units
from RCAIDE.Library.Plots                                           import *  
from RCAIDE.Library.Methods.Geometry.Planform.wing_segmented_planform      import segment_properties   
from RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor   import design_turbofan 
#from RCAIDE.Methods.Stability.Center_of_Gravity            import compute_component_centers_of_gravity


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
    # plot vehicle 
    plot_3d_vehicle(vehicle)       
    

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

    # ################################################# Vehicle-level Properties ########################################################  

    # mass properties
    vehicle.mass_properties.max_takeoff   = 63100  # kg 
    vehicle.mass_properties.takeoff       = 63100  # kg 
    vehicle.mass_properties.max_zero_fuel = 52200  # kg 
    vehicle.envelope.ultimate_load        = 3.75
    vehicle.envelope.limit_load           = 1.5
    vehicle.reference_area                = 112.3* Units['meters**2']
    vehicle.passengers                    = 135
    vehicle.systems.control               = "fully powered"
    vehicle.systems.accessories           = "medium range"   
    
    cruise_speed                          = 470 * Units.kts
    altitude                              = 40000 * Units.feet
    atmo                                  = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    freestream                            = atmo.compute_values (0.)
    freestream0                           = atmo.compute_values (altitude)
    mach_number                           = (cruise_speed/freestream.speed_of_sound)[0][0] 
    vehicle.design_dynamic_pressure       = ( .5 *freestream0.density*(cruise_speed*cruise_speed))[0][0]
    vehicle.design_mach_number            =  mach_number

   

    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------

    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing'
    wing.aspect_ratio                     = 9.24167
    wing.sweeps.quarter_chord             = 25 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.spans.projected                  = 33.70371 
    wing.chords.root                      = 7.3 * Units.meter
    wing.chords.tip                       = 0.4 * Units.meter
    wing.taper                            = wing.chords.tip / wing.chords.root
    wing.chords.mean_aerodynamic          = 4.88* Units.meter 
    wing.areas.reference                  = 122.915
    wing.areas.wetted                     = 390#? 
    wing.twists.root                      = 4.0 * Units.degrees 
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.origin                           = [[ 10.543,0,  -0.652]]
    wing.aerodynamic_center               = [0,0,0] 
    wing.vertical                         = False
    wing.dihedral                         = 7.0 * Units.degrees 
    wing.symmetric                        = True 
    wing.high_lift                        = True 
    wing.dynamic_pressure_ratio           = 1.0
        
    # Wing Segments
    root_airfoil                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator  + '..'  + separator 
    root_airfoil.coordinate_file          = rel_path  + 'Airfoils' + separator + 'B737a.txt'
    
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 0.0 * Units.degrees
    segment.root_chord_percent            = 1.
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 7.82609 * Units.degrees
    segment.sweeps.quarter_chord          = 24.48701 * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(root_airfoil)
    wing.append_segment(segment)

    
    #yehudi_airfoil                        = RCAIDE.Library.Components.Airfoils.Airfoil()
    #yehudi_airfoil.coordinate_file        = rel_path+ 'Airfoils' + separator + 'B737b.txt'
    #segment                               = RCAIDE.Library.Components.Wings.Segment()
    #segment.tag                           = 'Yehudi'
    #segment.percent_span_location         = 0.324
    #segment.twist                         = 0.047193 * Units.deg
    #segment.root_chord_percent            = 0.5
    #segment.thickness_to_chord            = 0.1
    #segment.dihedral_outboard             = 5.5 * Units.degrees
    #segment.sweeps.quarter_chord          = 25. * Units.degrees
    #segment.thickness_to_chord            = .1
    #segment.append_airfoil(yehudi_airfoil)
    #wing.append_segment(segment)

    mid_airfoil                           = RCAIDE.Library.Components.Airfoils.Airfoil()
    mid_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737c.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Section_2'
    segment.percent_span_location         = 0.368 
    segment.twist                         = 0.00258 * Units.deg
    segment.root_chord_percent            = 0.5136
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 6.41304 * Units.degrees
    segment.sweeps.quarter_chord          = 26.54545 * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(mid_airfoil)
    wing.append_segment(segment)

    tip_airfoil                           =  RCAIDE.Library.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737d.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Tip'
    segment.percent_span_location         = 0.96
    segment.twist                         = 5 * Units.degrees
    segment.root_chord_percent            = 0.1986
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 53.55* Units.degrees
    segment.sweeps.quarter_chord          = 48.0
    segment.thickness_to_chord            = .1
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)
    
    tip_airfoil                           =  RCAIDE.Library.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737d.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Tip'
    segment.percent_span_location         = 1.0
    segment.twist                         = 5 * Units.degrees
    segment.root_chord_percent            = 0.1066
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 53.55* Units.degrees
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = .1
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)    
    

    ## control surfaces -------------------------------------------
    #slat                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Slat()
    #slat.tag                      = 'slat'
    #slat.span_fraction_start      = 0.2
    #slat.span_fraction_end        = 0.963
    #slat.deflection               = 0.0 * Units.degrees
    #slat.chord_fraction           = 0.075
    #wing.append_control_surface(slat)

    #flap                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    #flap.tag                      = 'flap'
    #flap.span_fraction_start      = 0.2
    #flap.span_fraction_end        = 0.7
    #flap.deflection               = 0.0 * Units.degrees
    #flap.configuration_type       = 'double_slotted'
    #flap.chord_fraction           = 0.30
    #wing.append_control_surface(flap)

    #aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    #aileron.tag                   = 'aileron'
    #aileron.span_fraction_start   = 0.7
    #aileron.span_fraction_end     = 0.963
    #aileron.deflection            = 0.0 * Units.degrees
    #aileron.chord_fraction        = 0.16
    #wing.append_control_surface(aileron) 

   
    
    # add to vehicle
    vehicle.append_component(wing)





    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

    wing     = RCAIDE.Library.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'

    wing.aspect_ratio            = 2.55111
    wing.sweeps.quarter_chord    = 30
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.355 
    wing.spans.projected         = 12.30000
    wing.chords.root             = 3.55556
    wing.chords.tip              = 1.26587
    wing.chords.mean_aerodynamic = 2.27
    wing.areas.reference         = 16.47
    wing.areas.exposed           = 25.8 #?????
    wing.areas.wetted            = 29.65179#?
    wing.twists.root             = 3.0 * Units.degrees
    wing.twists.tip              = 3.0 * Units.degrees
    wing.origin                  = [[ 28.852,0,  1.230]]
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
    segment.dihedral_outboard      = 6.35714 * Units.degrees
    segment.sweeps.quarter_chord   = 30.00000 * Units.degrees 
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)

    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'tip_segment'
    segment.percent_span_location  = 1.
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 0.355             
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 0 * Units.degrees  
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)
    
    ## Fill out more segment properties automatically
    #wing = segment_properties(wing)        

    ## control surfaces -------------------------------------------
    #elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    #elevator.tag                   = 'elevator'
    #elevator.span_fraction_start   = 0.09
    #elevator.span_fraction_end     = 0.92
    #elevator.deflection            = 0.0  * Units.deg
    #elevator.chord_fraction        = 0.3
    #wing.append_control_surface(elevator)

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'

    wing.aspect_ratio            = 3.60335
    wing.sweeps.quarter_chord    = 38.92857
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.25

    wing.spans.projected         = 7.0709
    wing.total_length            = 7.0787

    wing.chords.root             = 6.27778
    wing.chords.tip              = 1.5714
    wing.chords.mean_aerodynamic = 4.0

    wing.areas.reference         = 30.83
    wing.areas.wetted            = 55.5

    wing.twists.root             = 3.0 * Units.degrees
    wing.twists.tip              = 3.0 * Units.degrees

    wing.origin                  = [[25.328,0, 1.311]]
    wing.aerodynamic_center      = [0,0,0]

    wing.vertical                = True
    wing.symmetric               = False
    wing.t_tail                  = False

    wing.dynamic_pressure_ratio  = 1.0


    # Wing Segments
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 1.
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = 38.92857 * Units.degrees  
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'segment_1'
    segment.percent_span_location         = 1.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.25
    segment.dihedral_outboard             = 0.0 * Units.degrees
    segment.sweeps.quarter_chord          = 0.0    
    segment.thickness_to_chord            = .1  
    wing.append_segment(segment)
    
 
    ## control surfaces -------------------------------------------
    #rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
    #rudder.tag                   = 'rudder'
    #rudder.span_fraction_start   = 0.1 
    #rudder.span_fraction_end     = 0.95 
    #rudder.deflection            = 0 
    #rudder.chord_fraction        = 0.33  
    #wing.append_control_surface(rudder)    
    
    
    ## Fill out more segment properties automatically
    #wing = segment_properties(wing)        


    # add to vehicle
    vehicle.append_component(wing)

    # ##########################################################   Fuselage  ############################################################    
    fuselage = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.seats_abreast                      = 5
    fuselage.fineness.nose                      = 0.58
    fuselage.fineness.tail                      = 1.75
    fuselage.lengths.nose                       = 0.921 
    fuselage.lengths.tail                       = 4.181
    fuselage.lengths.cabin                      = 29.798 #m
    fuselage.lengths.total                      = 34.9
    fuselage.width                              = 3.95 
    fuselage.heights.maximum                    = 4.19   
    fuselage.heights.at_quarter_length          = 3.71  
    fuselage.heights.at_three_quarters_length   = 3.83
    fuselage.heights.at_wing_root_quarter_chord = 4.19 
    fuselage.areas.side_projected               = fuselage.lengths.total *fuselage.heights.maximum  # estimate    
    fuselage.areas.wetted                       = 2 * np.pi * fuselage.width *  fuselage.lengths.total +  2 * np.pi * fuselage.width ** 2
    fuselage.areas.front_projected              =  np.pi * fuselage.width ** 2 
    fuselage.effective_diameter                 = 3.6 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0
    segment.percent_z_location                  = 0
    segment.height                              = 0
    segment.width                               = 0
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  = 0.01
    segment.percent_z_location                  = 0
    segment.height                              = 0.90000
    segment.width                               = 0.68182
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.02000
    segment.percent_z_location                  = 0.00100	 
    segment.height                              = 1.30000
    segment.width                               = 1.25000
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.03000
    segment.percent_z_location                  = 0.00300 
    segment.height                              = 1.65000
    segment.width                               = 1.64773
    fuselage.Segments.append(segment)  

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.04000 
    segment.percent_z_location                  = 0.00500	 
    segment.height                              = 1.96000
    segment.width                               = 1.93182
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.06000
    segment.percent_z_location                  = 0.01080 
    segment.height                              = 2.60000
    segment.width                               = 2.50000
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.07500 
    segment.percent_z_location                  = 0.01350	 
    segment.height                              = 2.92000
    segment.width                               = 2.67045 
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 0.10000
    segment.percent_z_location                  = 0.01750
    segment.height                              = 3.30000  
    segment.width                               = 3.01136
    fuselage.Segments.append(segment)
  

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_8'
    segment.percent_x_location                  = 0.16000
    segment.percent_z_location                  = 0.02174	 	 
    segment.height                              = 3.70000	 
    segment.width                               = 3.50000
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_9'
    segment.percent_x_location                  = 0.67000
    segment.percent_z_location                  = 0.02170
    segment.height                              = 3.70000
    segment.width                               = 3.50000
    fuselage.Segments.append(segment)
  

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_10'
    segment.percent_x_location                  = 0.73000
    segment.percent_z_location                  = 0.02600	 	 
    segment.height                              = 3.30000
    segment.width                               = 3.50000
    fuselage.Segments.append(segment)
    

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_11'
    segment.percent_x_location                  = 0.89582
    segment.percent_z_location                  = 0.04200	 	 
    segment.height                              = 1.95000 
    segment.width                               = 2.10000
    fuselage.Segments.append(segment)
    

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_12'
    segment.percent_x_location                  = 0.93715
    segment.percent_z_location                  = 0.04348 
    segment.height                              = 1.54000	 
    segment.width                               = 1.80000 
    fuselage.Segments.append(segment)
    
    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_13'
    segment.percent_x_location                  = 0.98463
    segment.percent_z_location                  = 0.04800
    segment.height                              = 0.85000	 
    segment.width                               = 0.80000 
    fuselage.Segments.append(segment)
    
    
    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_14'
    segment.percent_x_location                  = 1.0
    segment.percent_z_location                  = 0.04800
    segment.height                              = 0.0	 
    segment.width                               = 0.0000 
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
    net                                         = RCAIDE.Framework.Networks.Fuel() 
    
    #------------------------------------------------------------------------------------------------------------------------- 
    # Fuel Distrubition Line 
    #------------------------------------------------------------------------------------------------------------------------- 
    fuel_line                                   = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()  
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------         
    turbofan                                    = RCAIDE.Library.Components.Propulsors.Turbofan() 
    turbofan.tag                                = 'starboard_propulsor'
    turbofan.active_fuel_tanks                  = ['fuel_tank'] 
    turbofan.origin                             = [[ 10.150,  5.435, -1.087]] 
    turbofan.engine_length                      = 2.71     
    turbofan.bypass_ratio                       = 5.4    
    turbofan.design_altitude                    = 35000.0*Units.ft
    turbofan.design_mach_number                 = 0.78   
    turbofan.design_thrust                      = 35000.0* Units.N 

    # fan                
    fan                                         = RCAIDE.Library.Components.Propulsors.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.7   
    turbofan.fan                                = fan        

    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Library.Attributes.Gases.Air() 
    ram                                         = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    ram.tag                                     = 'ram' 
    turbofan.ram                                = ram 

    # inlet nozzle          
    inlet_nozzle                                = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                            = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency          = 0.98
    inlet_nozzle.pressure_ratio                 = 0.98 
    turbofan.inlet_nozzle                       = inlet_nozzle 

    # low pressure compressor    
    low_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    low_pressure_compressor.tag                   = 'lpc'
    low_pressure_compressor.polytropic_efficiency = 0.91
    low_pressure_compressor.pressure_ratio        = 1.9   
    turbofan.low_pressure_compressor              = low_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 10.0    
    turbofan.high_pressure_compressor              = high_pressure_compressor

    # low pressure turbine  
    low_pressure_turbine                           = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    low_pressure_turbine.tag                       ='lpt'
    low_pressure_turbine.mechanical_efficiency     = 0.99
    low_pressure_turbine.polytropic_efficiency     = 0.93 
    turbofan.low_pressure_turbine                  = low_pressure_turbine

    # high pressure turbine     
    high_pressure_turbine                          = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    high_pressure_turbine.tag                      ='hpt'
    high_pressure_turbine.mechanical_efficiency    = 0.99
    high_pressure_turbine.polytropic_efficiency    = 0.93 
    turbofan.high_pressure_turbine                 = high_pressure_turbine 

    # combustor  
    combustor                                      = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1500
    combustor.pressure_ratio                       = 0.95
    combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A1()  
    turbofan.combustor                             = combustor

    # core nozzle
    core_nozzle                                    = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    core_nozzle.tag                                = 'core nozzle'
    core_nozzle.polytropic_efficiency              = 0.95
    core_nozzle.pressure_ratio                     = 0.99  
    turbofan.core_nozzle                           = core_nozzle

    # fan nozzle             
    fan_nozzle                                     = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    fan_nozzle.tag                                 = 'fan nozzle'
    fan_nozzle.polytropic_efficiency               = 0.95
    fan_nozzle.pressure_ratio                      = 0.99 
    turbofan.fan_nozzle                            = fan_nozzle 

    # design turbofan
    design_turbofan(turbofan)  




    # Nacelle 
    nacelle                                     = RCAIDE.Library.Components.Nacelles.Body_of_Revolution_Nacelle()
    nacelle.diameter                            = 1.918
    nacelle.length                              = 3.258
    nacelle.tag                                 = 'nacelle_1'
    nacelle.inlet_diameter                      = 1.5 #?????
    nacelle.origin                              = [[ 10.150, 5.435, -1.087]] 
    nacelle.areas.wetted                        = 1.1*np.pi*nacelle.diameter*nacelle.length
    nacelle_airfoil                             = RCAIDE.Library.Components.Airfoils.NACA_4_Series_Airfoil()
    nacelle_airfoil.NACA_4_Series_code          = '2410'
    nacelle.Airfoil.append(nacelle_airfoil) 
    turbofan.nacelle                            = nacelle
    
    fuel_line.propulsors.append(turbofan)

        #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Port Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------      
    # copy turbofan
    turbofan_2                                  = deepcopy(turbofan)
    turbofan_2.active_fuel_tanks                = ['fuel_tank'] 
    turbofan_2.tag                              = 'port_propulsor' 
    turbofan_2.origin                           = [[10.150, -5.435, -1.087]]  # change origin 
    turbofan_2.nacelle.origin                   = [[10.150, -5.435, -1.087]]
         
    # append propulsor to distribution line 
    fuel_line.propulsors.append(turbofan_2)
  
    #------------------------------------------------------------------------------------------------------------------------- 
    #  Energy Source: Fuel Tank
    #------------------------------------------------------------------------------------------------------------------------- 
    # fuel tank
    fuel_tank                                   = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                            = wing.origin 
    
    # append fuel 
    fuel                                        = RCAIDE.Library.Attributes.Propellants.Jet_A1()   
    fuel.mass_properties.mass                   = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    fuel.origin                                 = vehicle.wings.main_wing.mass_properties.center_of_gravity      
    fuel.mass_properties.center_of_gravity      = vehicle.wings.main_wing.aerodynamic_center
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    fuel_tank.fuel                              = fuel            
    
    # apend fuel tank to dataclass of fuel tanks on fuel line 
    fuel_line.fuel_tanks.append(fuel_tank) 

    # Append fuel line to Network      
    net.fuel_lines.append(fuel_line)   

    # Append energy network to aircraft 
    vehicle.append_energy_network(net)    
    
    #------------------------------------------------------------------------------------------------------------------------- 
    # Compute Center of Gravity of aircraft (Optional)
    #------------------------------------------------------------------------------------------------------------------------- 
   
    vehicle.center_of_gravity()    
    #compute_component_centers_of_gravity(vehicle)
    
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