# Twin_Otter.py
# 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
# RCAIDE imports 
import RCAIDE      
from RCAIDE.Framework.Core import Units  
from   RCAIDE.Library.Methods.Propulsors.Turboprop_Propulsor        import design_turboprop  
from RCAIDE.Library.Methods.Performance.estimate_stall_speed        import estimate_stall_speed  
from RCAIDE.Library.Methods.Geometry.Planform                       import wing_segmented_planform  
from RCAIDE.Library.Plots                                           import *   

# python imports 
import numpy as np 
from copy import deepcopy
import os 
import matplotlib.pyplot        as plt 

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
    
    # plot vehicle 
    plot_3d_vehicle(vehicle,
                    min_x_axis_limit            = -5,
                    max_x_axis_limit            = 20,
                    min_y_axis_limit            = -10,
                    max_y_axis_limit            = 10,
                    min_z_axis_limit            = -10,
                    max_z_axis_limit            = 10)    
    
    return 
 
    
def vehicle_setup(): 
    

    #------------------------------------------------------------------------------------------------------------------------------------
    #   Initialize the Vehicle
    #------------------------------------------------------------------------------------------------------------------------------------

    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Twin_Otter'

 
    # ################################################# Vehicle-level Properties ########################################################  

    # mass properties
    vehicle.mass_properties.max_takeoff   = 5670  # kg 
    vehicle.mass_properties.takeoff       = 5670  # kg 
    vehicle.mass_properties.max_zero_fuel = 5000  # kg # INCORRECT  
    vehicle.flight_envelope.ultimate_load        = 5.7
    vehicle.flight_envelope.limit_load           = 3.8 
    vehicle.reference_area                = 39 
    vehicle.passengers                    = 19
    vehicle.systems.control               = "fully powered"
    vehicle.systems.accessories           = "commuter"    
    
    cruise_speed                          = 130 * Units.kts
    altitude                              = 5000 * Units.feet
    atmo                                  = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    freestream                            = atmo.compute_values (0.)
    freestream0                           = atmo.compute_values (altitude)
    mach_number                           = (cruise_speed/freestream.speed_of_sound)[0][0] 
    vehicle.design_dynamic_pressure       = ( .5 *freestream0.density*(cruise_speed*cruise_speed))[0][0]
    vehicle.flight_envelope.design_mach_number            =  mach_number

         
    # ##########################################################  Wings ################################################################    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Main Wing
    #------------------------------------------------------------------------------------------------------------------------------------
    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.sweeps.quarter_chord             = 0.0 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 39
    wing.areas.wetted                     = 39 * 2
    wing.spans.projected                  = 19.81
    wing.chords.root                      = 2.03 
    wing.chords.tip                       = 2.03 
    wing.chords.mean_aerodynamic          = 2.03 
    wing.taper                            = wing.chords.root/wing.chords.tip 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 3. * Units.degree 
    wing.twists.tip                       = 0
    wing.origin                           = [[5.38, 0, 1.35]] 
    wing.aerodynamic_center               = [5.38 + 0.25 *wing.chords.root , 0, 1.35]  
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0  
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator # CHANGED
    airfoil                               = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                           = 'Clark_y' 
    airfoil.coordinate_file               = rel_path + 'Airfoils' + separator + 'Clark_y.txt'   # absolute path     
    cg_x                                  = wing.origin[0][0] + 0.25*wing.chords.mean_aerodynamic
    cg_z                                  = wing.origin[0][2] - 0.2*wing.chords.mean_aerodynamic
    vehicle.mass_properties.center_of_gravity = [[cg_x,   0.  ,  cg_z ]]  # SOURCE: Design and aerodynamic analysis of a twin-engine commuter aircraft

    # Wing Segments
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'inboard'
    segment.percent_span_location         = 0.0 
    segment.twist                         = 3. * Units.degree 
    segment.root_chord_percent            = 1. 
    segment.dihedral_outboard             = 0. * Units.degree 
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)
    
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0
    segment.root_chord_percent            = 1.0
    segment.dihedral_outboard             = 0.
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)
    
    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.963
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.16
    wing.append_control_surface(aileron)
    
    
    # Fill out more segment properties automatically
    wing = wing_segmented_planform(wing)           
    
    # add to vehicle
    vehicle.append_component(wing)


    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Horizontal Tail
    #------------------------------------------------------------------------------------------------------------------------------------    
    wing                                  = RCAIDE.Library.Components.Wings.Horizontal_Tail()
    wing.tag                              = 'horizontal_stabilizer' 
    wing.sweeps.quarter_chord             = 0.0 * Units.degree
    wing.thickness_to_chord               = 0.12 
    wing.areas.reference                  = 9.762
    wing.areas.wetted                     = 20.0
    wing.spans.projected                  = 6.29   
    wing.chords.root                      = 1.552 
    wing.chords.tip                       = 1.552 
    wing.chords.mean_aerodynamic          = 1.552  
    wing.taper                            = 0 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 0.0 * Units.degree
    wing.twists.tip                       = 0.0 * Units.degree 
    wing.origin                           = [[13.17 , 0 , 1.25]] 
    wing.aerodynamic_center               = [13.17 , 0 , 1.25]  
    wing.vertical                         = False
    wing.winglet_fraction                 = 0.0  
    wing.symmetric                        = True
    wing.high_lift                        = False 
    wing.dynamic_pressure_ratio           = 0.9

    # control surfaces -------------------------------------------
    elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.09
    elevator.span_fraction_end     = 0.92
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    wing.append_control_surface(elevator)
    
    # add to vehicle
    vehicle.append_component(wing)


    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Vertical Stabilizer
    #------------------------------------------------------------------------------------------------------------------------------------ 
    wing                                  = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag                              = 'vertical_stabilizer'     
    wing.sweeps.leading_edge              = 28.6 * Units.degree 
    wing.thickness_to_chord               = 0.12 
    wing.areas.reference                  = 8.753
    wing.areas.wetted                     = 18.0   
    wing.spans.projected                  = 3.5
    wing.chords.root                      = 2.975 
    wing.chords.tip                       = 1.514
    wing.chords.mean_aerodynamic          = 2.24 # incorrect 
    wing.taper                            = wing.chords.tip/wing.chords.root 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 0.0 * Units.degree
    wing.twists.tip                       = 0.0 * Units.degree 
    wing.origin                           = [[ 12.222 , 0 , 0.385 ]]  
    wing.aerodynamic_center               = [[ 12.222 + 0.25 * wing.chords.root, 0 , 0.385 ]]  
    wing.vertical                         = True 
    wing.symmetric                        = False
    wing.t_tail                           = False
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0
    

    rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
    rudder.tag                   = 'rudder'
    rudder.span_fraction_start   = 0.1
    rudder.span_fraction_end     = 0.963
    rudder.deflection            = 0.0 * Units.degrees
    rudder.chord_fraction        = 0.3
    
    wing.append_control_surface(rudder)    

    # add to vehicle
    vehicle.append_component(wing)

 
    # ##########################################################   Fuselage  ############################################################    
    fuselage = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.seats_abreast                      = 2.
    fuselage.number_coach_seats                 = 19
    fuselage.fineness.nose                      = 1.6
    fuselage.fineness.tail                      = 2.
    fuselage.lengths.nose                       = 2.95  
    fuselage.lengths.tail                       = 7.57
    fuselage.lengths.cabin                      = 4.62 
    fuselage.lengths.total                      = 15.77  
    fuselage.width                              = 1.75  
    fuselage.heights.maximum                    = 1.50  
    fuselage.heights.at_quarter_length          = 1.50  
    fuselage.heights.at_three_quarters_length   = 1.50  
    fuselage.heights.at_wing_root_quarter_chord = 1.50  
    fuselage.areas.side_projected               = fuselage.lengths.total *fuselage.heights.maximum  # estimate    
    fuselage.areas.wetted                       = 2 * np.pi * fuselage.width *  fuselage.lengths.total +  2 * np.pi * fuselage.width ** 2
    fuselage.areas.front_projected              =  np.pi * fuselage.width ** 2 
    fuselage.effective_diameter                 = 1.75 

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
    segment.percent_x_location                  = 0.005345402
    segment.percent_z_location                  = -0.027433333/ fuselage.lengths.total	 
    segment.height                              = 0.421666667
    segment.width                               = 0.106025757
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.019706071
    segment.percent_z_location                  = 6.66667E-05/ fuselage.lengths.total	 
    segment.height                              = 0.733333333
    segment.width                               = 0.61012023
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.054892307
    segment.percent_z_location                  = 0.009233333/ fuselage.lengths.total	 
    segment.height                              = 1.008333333
    segment.width                               = 1.009178159
    fuselage.Segments.append(segment)  

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.082575704 
    segment.percent_z_location                  = 0.0459 / fuselage.lengths.total	 
    segment.height                              = 1.228333333 
    segment.width                               = 1.275456588 
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.116879689
    segment.percent_z_location                  = 0.055066667/ fuselage.lengths.total	 
    segment.height                              = 1.393333333
    segment.width                               = 1.436068974
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.151124363 
    segment.percent_z_location                  = 0.091733333 / fuselage.lengths.total	 
    segment.height                              = 1.576666667 
    segment.width                               = 1.597041074  
    fuselage.Segments.append(segment) 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 0.172884077 
    segment.percent_z_location                  = 0.251733333 / fuselage.lengths.total	 
    segment.height                              = 1.953333333  
    segment.width                               = 1.75  
    fuselage.Segments.append(segment)
  

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_8'
    segment.percent_x_location                  = 0.194554824   
    segment.percent_z_location                  = 0.311733333/ fuselage.lengths.total	 	 
    segment.height                              = 2.09	 
    segment.width                               = 1.75 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_9'
    segment.percent_x_location                  =  0.479203019 
    segment.percent_z_location                  = 0.311733333 / fuselage.lengths.total	  
    segment.height                              = 2.09	 
    segment.width                               = 1.75 
    fuselage.Segments.append(segment)
    
    

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_10'
    segment.percent_x_location                  = 0.541657475 
    segment.percent_z_location                  = 0.311733333/ fuselage.lengths.total	 	 
    segment.height                              = 2.09	 
    segment.width                               = 2 
    fuselage.Segments.append(segment)
    

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_11'
    segment.percent_x_location                  = 0.716936232 
    segment.percent_z_location                  = 0.394233333/ fuselage.lengths.total	 	 
    segment.height                              = 1.558333333	 
    segment.width                               = 0.64
    fuselage.Segments.append(segment)
    

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_12'
    segment.percent_x_location                  = 1.0
    segment.percent_z_location                  = 0.440066667 / fuselage.lengths.total	 
    segment.height                              = 0.11	 
    segment.width                               = 0.05 
    fuselage.Segments.append(segment) 
          
 

    # add to vehicle
    vehicle.append_component(fuselage)

    # ------------------------------------------------------------------
    #   Landing gear
    # ------------------------------------------------------------------  
    landing_gear                                = RCAIDE.Library.Components.Landing_Gear.Landing_Gear()
    main_gear                                   = RCAIDE.Library.Components.Landing_Gear.Main_Landing_Gear()
    nose_gear                                   = RCAIDE.Library.Components.Landing_Gear.Nose_Landing_Gear()
    main_gear.strut_length                      = 12. * Units.inches  
    nose_gear.strut_length                      = 6. * Units.inches 
                                                
    landing_gear.main                           = main_gear
    landing_gear.nose                           = nose_gear
                                                
    #add to vehicle                             
    vehicle.landing_gear                        = landing_gear

    # ########################################################  Energy Network  #########################################################  
    net                                         = RCAIDE.Framework.Networks.Fuel()    

    #------------------------------------------------------------------------------------------------------------------------- 
    # Fuel Distrubition Line 
    #------------------------------------------------------------------------------------------------------------------------- 
    fuel_line                                       = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()  
 
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------    
    starboard_propulsor                              = RCAIDE.Library.Components.Propulsors.Turboprop()    
    starboard_propulsor.tag                          = 'starboard_propulsor' 
    starboard_propulsor.active_fuel_tanks            = ['fuel_tank']   
    starboard_propulsor.origin                       = [[3.5,2.8129,1.22 ]]
    starboard_propulsor.design_altitude              = 25000*Units.ft                                   # [-]         Design Altitude
    starboard_propulsor.design_mach_number           = 0.5                                              # [-]         Design Mach number
    starboard_propulsor.design_propeller_efficiency  = 0.86
    starboard_propulsor.design_thrust                = 23000 * Units.N                                  # [-]         Design Thrust 
    starboard_propulsor.working_fluid                = RCAIDE.Library.Attributes.Gases.Air()            
    starboard_propulsor.design_propeller_efficiency  = 0.83                                             # [-]         Design Propeller Efficiency
    starboard_propulsor.design_gearbox_efficiency    = 0.99                                             # [-]         Design Gearbox Efficiency
    
    # Ram inlet 
    ram                                              = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    ram.tag                                          = 'ram' 
    starboard_propulsor.ram                          = ram 
          
    # inlet nozzle          
    inlet_nozzle                                     = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                                 = 'inlet nozzle'                                       
    inlet_nozzle.pressure_ratio                      = 0.98
    inlet_nozzle.compressibility_effects             = False
    starboard_propulsor.inlet_nozzle                 = inlet_nozzle
                                                     
    # compressor                        
    compressor                                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    compressor.tag                                   = 'lpc'                   
    compressor.pressure_ratio                        = 10                   
    starboard_propulsor.compressor                   = compressor
    
    # combustor      
    combustor                                        = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                    = 'Comb'
    combustor.efficiency                             = 0.99                   
    combustor.turbine_inlet_temperature              = 1370                    
    combustor.pressure_ratio                         = 0.96                    
    combustor.fuel_data                              = RCAIDE.Library.Attributes.Propellants.Jet_A()  
    starboard_propulsor.combustor                    = combustor
        
    # high pressure turbine         
    high_pressure_turbine                            = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    high_pressure_turbine.tag                        ='hpt'
    high_pressure_turbine.mechanical_efficiency      = 0.99                       
    starboard_propulsor.high_pressure_turbine        = high_pressure_turbine 
        
    # low pressure turbine      
    low_pressure_turbine                             = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    low_pressure_turbine.tag                         ='lpt'
    low_pressure_turbine.mechanical_efficiency       = 0.99                      
    starboard_propulsor.low_pressure_turbine         = low_pressure_turbine
    
    # core nozzle    
    core_nozzle                                      = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    core_nozzle.tag                                  = 'core nozzle'          
    core_nozzle.pressure_ratio                       = 0.99
    starboard_propulsor.core_nozzle                  = core_nozzle
    
    # design starboard_propulsor
    design_turboprop(starboard_propulsor)
    

 
    #########################################################   Nacelles  ############################################################    
    nacelle                    = RCAIDE.Library.Components.Nacelles.Stack_Nacelle()
    nacelle.tag                = 'nacelle_1'
    nacelle.length             = 4
    nacelle.diameter           = 0.73480616 
    nacelle.areas.wetted       = 0.01*(2*np.pi*0.01/2)
    nacelle.origin             = [[3.5,2.8129,1]]
    nacelle.flow_through       = False  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.042687938 
    nac_segment.percent_z_location = 0.0284/ nacelle.length
    nac_segment.height             = 0.183333333 
    nac_segment.width              = 0.422484315 
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.143080714 
    nac_segment.percent_z_location = 0.046733333/ nacelle.length
    nac_segment.height             = 0.44	 
    nac_segment.width              = 0.685705173 
    nacelle.append_segment(nac_segment)  
     
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.170379029  
    nac_segment.percent_z_location = -0.154233333/ nacelle.length
    nac_segment.height             = 0.898333333	 
    nac_segment.width              = 0.73480616 
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.252189893  
    nac_segment.percent_z_location = -0.154233333/ nacelle.length
    nac_segment.height             = 1.008333333 
    nac_segment.width              = 0.736964445
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.383860821   
    nac_segment.percent_z_location = -0.072566667/ nacelle.length
    nac_segment.height             = 0.971666667 
    nac_segment.width              = 0.736964445 
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 0.551826736  
    nac_segment.percent_z_location = .055066667/ nacelle.length	
    nac_segment.height             = 0.77	 
    nac_segment.width              = 0.736964445  
    nacelle.append_segment(nac_segment)
    

    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_8'
    nac_segment.percent_x_location = 0.809871485   
    nac_segment.percent_z_location = 0.1284/ nacelle.length
    nac_segment.height             = 0.366666667 
    nac_segment.width              = 0.736964445 
    nacelle.append_segment(nac_segment) 
    

    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_9'
    nac_segment.percent_x_location = 1.0  
    nac_segment.percent_z_location = 0.201733333 / nacelle.length
    nac_segment.height             = 0.036666667	
    nac_segment.width              = 0.0  
    nacelle.append_segment(nac_segment) 
    
    starboard_propulsor.nacelle = nacelle  
 
    fuel_line.propulsors.append(starboard_propulsor)    

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Port Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------      
    # copy turboprop
    port_propulsor                                  = deepcopy(starboard_propulsor)
    port_propulsor.active_fuel_tanks                = ['fuel_tank'] 
    port_propulsor.tag                              = 'port_propulsor' 
    port_propulsor.origin                           = [[3.5, -2.8129,1.22 ]]  # change origin 
    port_propulsor.nacelle.tag                      = 'port_propulsor_nacelle' 
    port_propulsor.nacelle.origin                   = [[3.5, -2.8129,1.22 ]]
         
    # append propulsor to distribution line 
    fuel_line.propulsors.append(port_propulsor) 

    #------------------------------------------------------------------------------------------------------------------------- 
    #  Energy Source: Fuel Tank
    #------------------------------------------------------------------------------------------------------------------------- 
    # fuel tank
    fuel_tank                                        = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                                 = vehicle.wings.main_wing.origin  
    fuel_tank.fuel                                   = RCAIDE.Library.Attributes.Propellants.Jet_A1()   
    fuel_tank.fuel.mass_properties.mass              = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    fuel_tank.fuel.origin                            = vehicle.wings.main_wing.mass_properties.center_of_gravity      
    fuel_tank.fuel.mass_properties.center_of_gravity = vehicle.wings.main_wing.aerodynamic_center
    fuel_tank.volume                                 = fuel_tank.fuel.mass_properties.mass/fuel_tank.fuel.density   
    
    # apend fuel tank to dataclass of fuel tanks on fuel line 
    fuel_line.fuel_tanks.append(fuel_tank) 

    # Append fuel line to Network      
    net.fuel_lines.append(fuel_line)
    
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Avionics
    #------------------------------------------------------------------------------------------------------------------------------------ 
    Wuav                                        = 20. * Units.lbs
    avionics                                    = RCAIDE.Library.Components.Systems.Avionics()
    avionics.mass_properties.uninstalled        = Wuav
    vehicle.avionics                            = avionics     
    
    
    vehicle.append_energy_network(net)   

     
    return vehicle
  
# ---------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------
def configs_setup(vehicle):

    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs         = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config     = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'  
    configs.append(base_config)

 
    return configs



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
    aerodynamics.settings.number_of_spanwise_vortices   = 70
    aerodynamics.settings.number_of_chordwise_vortices  = 5   
    analyses.append(aerodynamics)
 
    # ------------------------------------------------------------------
    #  Energy
    energy = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle
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
    segment.analyses.extend( analyses.base ) 
    segment.altitude                                      = 1.0 * Units.km  
    segment.air_speed                                     = 77 * Units['m/s']
    segment.distance                                      = 1000 * Units.nmi   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
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
