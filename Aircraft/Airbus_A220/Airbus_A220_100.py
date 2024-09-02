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
#from RCAIDE.Library.Methods.Geometry.Two_Dimensional.Planform      import segment_properties   
#from RCAIDE.Methods.Energy.Propulsors.Turbofan_Propulsor   import design_turbofan , compute_nacelle_geometry
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

   

    ## ------------------------------------------------------------------
    ##   Main Wing
    ## ------------------------------------------------------------------

    #wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    #wing.tag                              = 'main_wing'
    #wing.aspect_ratio                     = 10.97
    #wing.sweeps.quarter_chord             = 25 * Units.deg
    #wing.thickness_to_chord               = 0.12
    #wing.spans.projected                  = 35.1 
    #wing.chords.root                      = 5.848 * Units.meter
    #wing.chords.tip                       = 0.554 * Units.meter
    #wing.taper                            = wing.chords.tip / wing.chords.root
    #wing.chords.mean_aerodynamic          = * Units.meter 
    #wing.areas.reference                  = 112.3
    #wing.areas.wetted                     = 390 
    #wing.twists.root                      = 4.0 * Units.degrees 
    #wing.twists.tip                       = 0.0 * Units.degrees 
    #wing.origin                           = [[12.96441762,0, -0.170]]
    #wing.aerodynamic_center               = [0,0,0] 
    #wing.vertical                         = False
    #wing.dihedral                         = 3.5 * Units.degrees 
    #wing.symmetric                        = True 
    #wing.high_lift                        = True 
    #wing.dynamic_pressure_ratio           = 1.0
        
    ## Wing Segments
    #root_airfoil                          = RCAIDE.Library.Components.Airfoils.Airfoil()
    #ospath                                = os.path.abspath(__file__)
    #separator                             = os.path.sep
    #rel_path                              = os.path.dirname(ospath) + separator  + '..'  + separator 
    #root_airfoil.coordinate_file          = rel_path  + 'Airfoils' + separator + 'B737a.txt' # what kind of airfoil for a220?
    #segment                               = RCAIDE.Library.Components.Wings.Segment()
    #segment.tag                           = 'Root'
    #segment.percent_span_location         = 0.0
    #segment.twist                         = 4. * Units.deg #...?
    #segment.root_chord_percent            = 1.
    #segment.thickness_to_chord            = 0.1
    #segment.dihedral_outboard             = 2.5 * Units.degrees
    #segment.sweeps.quarter_chord          = 0 * Units.degrees
    #segment.thickness_to_chord            = .1
    #segment.append_airfoil(root_airfoil)
    #wing.append_segment(segment)

    
    ##yehudi_airfoil                        = RCAIDE.Library.Components.Airfoils.Airfoil()
    ##yehudi_airfoil.coordinate_file        = rel_path+ 'Airfoils' + separator + 'B737b.txt'
    ##segment                               = RCAIDE.Library.Components.Wings.Segment()
    ##segment.tag                           = 'Yehudi'
    ##segment.percent_span_location         = 0.324
    ##segment.twist                         = 0.047193 * Units.deg
    ##segment.root_chord_percent            = 0.5
    ##segment.thickness_to_chord            = 0.1
    ##segment.dihedral_outboard             = 5.5 * Units.degrees
    ##segment.sweeps.quarter_chord          = 25. * Units.degrees
    ##segment.thickness_to_chord            = .1
    ##segment.append_airfoil(yehudi_airfoil)
    ##wing.append_segment(segment)

    #mid_airfoil                           = RCAIDE.Library.Components.Airfoils.Airfoil()
    #mid_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737c.txt'
    #segment                               = RCAIDE.Library.Components.Wings.Segment()
    #segment.tag                           = 'Section_2'
    #segment.percent_span_location         = 0.963
    #segment.twist                         = 0.00258 * Units.deg
    #segment.root_chord_percent            = 0.220
    #segment.thickness_to_chord            = 0.1
    #segment.dihedral_outboard             = 5.5 * Units.degrees
    #segment.sweeps.quarter_chord          = 56.75 * Units.degrees
    #segment.thickness_to_chord            = .1
    #segment.append_airfoil(mid_airfoil)
    #wing.append_segment(segment)

    #tip_airfoil                           =  RCAIDE.Library.Components.Airfoils.Airfoil()
    #tip_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737d.txt'
    #segment                               = RCAIDE.Library.Components.Wings.Segment()
    #segment.tag                           = 'Tip'
    #segment.percent_span_location         = 1.
    #segment.twist                         = 0. * Units.degrees
    #segment.root_chord_percent            = 0.10077
    #segment.thickness_to_chord            = 0.1
    #segment.dihedral_outboard             = 0.
    #segment.sweeps.quarter_chord          = 0.
    #segment.thickness_to_chord            = .1
    #segment.append_airfoil(tip_airfoil)
    #wing.append_segment(segment)
    

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

    ## Fill out more segment properties automatically
    #wing = segment_properties(wing)   
    
    ## add to vehicle
    #vehicle.append_component(wing)





    ## ------------------------------------------------------------------
    ##  Horizontal Stabilizer
    ## ------------------------------------------------------------------

    #wing     = RCAIDE.Components.Wings.Horizontal_Tail()
    #wing.tag = 'horizontal_stabilizer'

    #wing.aspect_ratio            = 
    #wing.sweeps.quarter_chord    = 
    #wing.thickness_to_chord      = 
    #wing.taper                   = 
    #wing.spans.projected         = 
    #wing.chords.root             = 
    #wing.chords.tip              = 
    #wing.chords.mean_aerodynamic = 
    #wing.areas.reference         = 
    #wing.areas.exposed           = 
    #wing.areas.wetted            = 
    #wing.twists.root             = 
    #wing.twists.tip              = 
    #wing.origin                  = 
    #wing.aerodynamic_center      = 
    #wing.vertical                = 
    #wing.symmetric               = 
    #wing.dynamic_pressure_ratio  = 





    ## add to vehicle
    #vehicle.append_component(wing)


    ## ------------------------------------------------------------------
    ##   Vertical Stabilizer
    ## ------------------------------------------------------------------

    #wing = RCAIDE.Components.Wings.Vertical_Tail()
    #wing.tag = 'vertical_stabilizer'

    #wing.aspect_ratio            = 
    #wing.sweeps.quarter_chord    = 
    #wing.thickness_to_chord      = 
    #wing.taper                   = 

    #wing.spans.projected         = 
    #wing.total_length            = 

    #wing.chords.root             = 
    #wing.chords.tip              = 
    #wing.chords.mean_aerodynamic = 

    #wing.areas.reference         = 
    #wing.areas.wetted            = 

    #wing.twists.root             = 
    #wing.twists.tip              = 

    #wing.origin                  = 
    #wing.aerodynamic_center      = 

    #wing.vertical                = 
    #wing.symmetric               = 
    #wing.t_tail                  = 

    #wing.dynamic_pressure_ratio  = 



    ## add to vehicle
    #vehicle.append_component(wing)

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
    segment.height                              = 0.0134
    segment.width                               = 0.017
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  = 0.016231
    segment.percent_z_location                  = -0.027433333/ fuselage.lengths.total	 
    segment.height                              = 1.556338028
    segment.width                               = 1.436619718
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.048379142
    segment.percent_z_location                  = 6.66667E-05/ fuselage.lengths.total	 
    segment.height                              = 2.633802817
    segment.width                               = 2.274647887
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.099627853
    segment.percent_z_location                  = 0.009233333/ fuselage.lengths.total	 
    segment.height                              = 3.471830986
    segment.width                               = 3.471830986
    fuselage.Segments.append(segment)  

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.15379097 
    segment.percent_z_location                  = 0.0459 / fuselage.lengths.total	 
    segment.height                              = 3.711267606
    segment.width                               = 3.711267606 
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.548580908
    segment.percent_z_location                  = 0.055066667/ fuselage.lengths.total	 
    segment.height                              = 4.190140845
    segment.width                               = 3.950704225
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.602699188 
    segment.percent_z_location                  = 0.091733333 / fuselage.lengths.total	 
    segment.height                              = 3.950704225 
    segment.width                               = 3.830985915  
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 0.774604313
    segment.percent_z_location                  = 0.251733333 / fuselage.lengths.total	 
    segment.height                              = 3.830985915  
    segment.width                               = 3.711267606
    fuselage.Segments.append(segment)
  

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_8'
    segment.percent_x_location                  = 0.870017486 
    segment.percent_z_location                  = 0.311733333/ fuselage.lengths.total	 	 
    segment.height                              = 2.992957746	 
    segment.width                               = 2.633802817
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_9'
    segment.percent_x_location                  = 0.952696947
    segment.percent_z_location                  = 0.311733333 / fuselage.lengths.total	  
    segment.height                              = 1.915492958	 
    segment.width                               = 1.676056338 
    fuselage.Segments.append(segment)
    
    

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_10'
    segment.percent_x_location                  = 0.984441555
    segment.percent_z_location                  = 0.311733333/ fuselage.lengths.total	 	 
    segment.height                              = 1.197183099 
    segment.width                               = 0.838028169 
    fuselage.Segments.append(segment)
    

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_11'
    segment.percent_x_location                  = 1.0 
    segment.percent_z_location                  = 0.394233333/ fuselage.lengths.total	 	 
    segment.height                              = 0.838028169	 
    segment.width                               = 0.478873239
    fuselage.Segments.append(segment)
    

    ## Segment
    #segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    #segment.tag                                 = 'segment_12'
    #segment.percent_x_location                  = 1.0
    #segment.percent_z_location                  = 0.440066667 / fuselage.lengths.total	 
    #segment.height                              = 0.11	 
    #segment.width                               = 0.05 
    #fuselage.Segments.append(segment) 
          
 

    # add to vehicle
    vehicle.append_component(fuselage) 


    ## ################################################# Energy Network #######################################################         
    ## Step 1: Define network
    ## Step 2: Define Distribution Type
    ## Step 3: Define Propulsors 
    ## Step 4: Define Enegy Source 

    ##------------------------------------------------------------------------------------------------------------------------- 
    ##  Turbofan Network
    ##-------------------------------------------------------------------------------------------------------------------------   
    #net                                         = RCAIDE.Energy.Networks.Turbofan_Engine_Network() 

    ## Append energy network to aircraft 
    #vehicle.append_energy_network(net)   

    ##------------------------------------------------------------------------------------------------------------------------- 
    ## Fuel Distrubition Line 
    ##------------------------------------------------------------------------------------------------------------------------- 
    #fuel_line                                   = RCAIDE.Energy.Networks.Distribution.Fuel_Line()  

    ##------------------------------------------------------------------------------------------------------------------------------------  
    ## Propulsor: Starboard Propulsor
    ##------------------------------------------------------------------------------------------------------------------------------------         
    #turbofan                                    = RCAIDE.Energy.Propulsors.Turbofan() 
    #turbofan.tag                                = 'starboard_propulsor'
    #turbofan.active_fuel_tanks                  = ['fuel_tank']   
    #turbofan.origin                             = [[13.72, 4.86,-1.1]] 
    #turbofan.engine_length                      = 2.71     
    #turbofan.bypass_ratio                       = 5.4    
    #turbofan.design_altitude                    = 35000.0*Units.ft
    #turbofan.design_mach_number                 = 0.78   
    #turbofan.design_thrust                      = 35000.0* Units.N 

    ## fan                
    #fan                                         = RCAIDE.Energy.Propulsors.Converters.Fan()   
    #fan.tag                                     = 'fan'
    #fan.polytropic_efficiency                   = 0.93
    #fan.pressure_ratio                          = 1.7   
    #turbofan.fan                                = fan        

    ## working fluid                   
    #turbofan.working_fluid                      = RCAIDE.Attributes.Gases.Air() 
    #ram                                         = RCAIDE.Energy.Propulsors.Converters.Ram()
    #ram.tag                                     = 'ram' 
    #turbofan.ram                                = ram 

    ## inlet nozzle          
    #inlet_nozzle                                = RCAIDE.Energy.Propulsors.Converters.Compression_Nozzle()
    #inlet_nozzle.tag                            = 'inlet nozzle'
    #inlet_nozzle.polytropic_efficiency          = 0.98
    #inlet_nozzle.pressure_ratio                 = 0.98 
    #turbofan.inlet_nozzle                       = inlet_nozzle 

    ## low pressure compressor    
    #low_pressure_compressor                       = RCAIDE.Energy.Propulsors.Converters.Compressor()    
    #low_pressure_compressor.tag                   = 'lpc'
    #low_pressure_compressor.polytropic_efficiency = 0.91
    #low_pressure_compressor.pressure_ratio        = 1.9   
    #turbofan.low_pressure_compressor              = low_pressure_compressor

    ## high pressure compressor  
    #high_pressure_compressor                       = RCAIDE.Energy.Propulsors.Converters.Compressor()    
    #high_pressure_compressor.tag                   = 'hpc'
    #high_pressure_compressor.polytropic_efficiency = 0.91
    #high_pressure_compressor.pressure_ratio        = 10.0    
    #turbofan.high_pressure_compressor              = high_pressure_compressor

    ## low pressure turbine  
    #low_pressure_turbine                           = RCAIDE.Energy.Propulsors.Converters.Turbine()   
    #low_pressure_turbine.tag                       ='lpt'
    #low_pressure_turbine.mechanical_efficiency     = 0.99
    #low_pressure_turbine.polytropic_efficiency     = 0.93 
    #turbofan.low_pressure_turbine                  = low_pressure_turbine

    ## high pressure turbine     
    #high_pressure_turbine                          = RCAIDE.Energy.Propulsors.Converters.Turbine()   
    #high_pressure_turbine.tag                      ='hpt'
    #high_pressure_turbine.mechanical_efficiency    = 0.99
    #high_pressure_turbine.polytropic_efficiency    = 0.93 
    #turbofan.high_pressure_turbine                 = high_pressure_turbine 

    ## combustor  
    #combustor                                      = RCAIDE.Energy.Propulsors.Converters.Combustor()   
    #combustor.tag                                  = 'Comb'
    #combustor.efficiency                           = 0.99 
    #combustor.alphac                               = 1.0     
    #combustor.turbine_inlet_temperature            = 1500
    #combustor.pressure_ratio                       = 0.95
    #combustor.fuel_data                            = RCAIDE.Attributes.Propellants.Jet_A()  
    #turbofan.combustor                             = combustor

    ## core nozzle
    #core_nozzle                                    = RCAIDE.Energy.Propulsors.Converters.Expansion_Nozzle()   
    #core_nozzle.tag                                = 'core nozzle'
    #core_nozzle.polytropic_efficiency              = 0.95
    #core_nozzle.pressure_ratio                     = 0.99  
    #turbofan.core_nozzle                           = core_nozzle

    ## fan nozzle             
    #fan_nozzle                                     = RCAIDE.Energy.Propulsors.Converters.Expansion_Nozzle()   
    #fan_nozzle.tag                                 = 'fan nozzle'
    #fan_nozzle.polytropic_efficiency               = 0.95
    #fan_nozzle.pressure_ratio                      = 0.99 
    #turbofan.fan_nozzle                            = fan_nozzle 

    ## design turbofan
    #design_turbofan(turbofan)  
    ## append propulsor to distribution line 



    ## Nacelle 
    #nacelle                                     = RCAIDE.Components.Nacelles.Nacelle()
    #nacelle.diameter                            = 2.05
    #nacelle.length                              = 2.71
    #nacelle.tag                                 = 'nacelle_1'
    #nacelle.inlet_diameter                      = 2.0
    #nacelle.origin                              = [[13.5,4.38,-1.5]] 
    #nacelle.areas.wetted                        = 1.1*np.pi*nacelle.diameter*nacelle.length
    #nacelle.Airfoil.NACA_4_series_flag          = True 
    #nacelle.Airfoil.coordinate_file             = '2410' 


    #nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    #nac_segment.tag                             = 'segment_1'
    #nac_segment.percent_x_location              = 0.0  
    #nac_segment.height                          = 2.05
    #nac_segment.width                           = 2.05
    #nacelle.append_segment(nac_segment)         

    #nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    #nac_segment.tag                             = 'segment_2'
    #nac_segment.percent_x_location              = 0.3
    #nac_segment.height                          = 2.1  
    #nac_segment.width                           = 2.1 
    #nacelle.append_segment(nac_segment)         

    #nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    #nac_segment.tag                             = 'segment_3'
    #nac_segment.percent_x_location              = 0.4  
    #nac_segment.height                          = 2.05
    #nac_segment.width                           = 2.05 
    #nacelle.append_segment(nac_segment)         

    #nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    #nac_segment.tag                             = 'segment_4'
    #nac_segment.percent_x_location              = 0.75  
    #nac_segment.height                          = 1.9
    #nac_segment.width                           = 1.9
    #nacelle.append_segment(nac_segment)         

    #nac_segment                                 = RCAIDE.Components.Lofted_Body_Segment.Segment()
    #nac_segment.tag                             = 'segment_5'
    #nac_segment.percent_x_location              = 1.0
    #nac_segment.height                          = 1.7 
    #nac_segment.width                           = 1.7
    #nacelle.append_segment(nac_segment)           
    #compute_nacelle_geometry(turbofan, nacelle)
    #turbofan.nacelle                            = nacelle

    #fuel_line.propulsors.append(turbofan)  

    ##------------------------------------------------------------------------------------------------------------------------------------  
    ## Propulsor: Port Propulsor
    ##------------------------------------------------------------------------------------------------------------------------------------      
    ## copy turbofan
    #turbofan_2                                  = deepcopy(turbofan)
    #turbofan_2.active_fuel_tanks                = ['fuel_tank'] 
    #turbofan_2.tag                              = 'port_propulsor' 
    #turbofan_2.origin                           = [[13.72,-4.38,-1.1]]  # change origin 
    #turbofan_2.nacelle.origin                   = [[13.5,-4.38,-1.5]]

    ## append propulsor to distribution line 
    #fuel_line.propulsors.append(turbofan_2)

    ##------------------------------------------------------------------------------------------------------------------------- 
    ##  Energy Source: Fuel Tank
    ##------------------------------------------------------------------------------------------------------------------------- 
    ## fuel tank
    #fuel_tank                                   = RCAIDE.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    #fuel_tank.origin                            = wing.origin 

    ## append fuel 
    #fuel                                        = RCAIDE.Attributes.Propellants.Aviation_Gasoline()   
    #fuel.mass_properties.mass                   = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    #fuel.origin                                 = vehicle.wings.main_wing.mass_properties.center_of_gravity      
    #fuel.mass_properties.center_of_gravity      = vehicle.wings.main_wing.aerodynamic_center
    #fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    #fuel_tank.fuel                              = fuel            

    ## apend fuel tank to dataclass of fuel tanks on fuel line 
    #fuel_line.fuel_tanks.append(fuel_tank) 

    ## Append fuel line to Network      
    #net.fuel_lines.append(fuel_line)  

    ##------------------------------------------------------------------------------------------------------------------------- 
    ## Compute Center of Gravity of aircraft (Optional)
    ##------------------------------------------------------------------------------------------------------------------------- 

    #vehicle.center_of_gravity()    
    #compute_component_centers_of_gravity(vehicle)

    ##------------------------------------------------------------------------------------------------------------------------- 
    ## Done ! 
    ##------------------------------------------------------------------------------------------------------------------------- 

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