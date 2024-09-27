# Twin_Otter.py
# 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
# RCAIDE imports 
import RCAIDE      
from RCAIDE.Framework.Core import Units  
from RCAIDE.Framework.Networks.Electric                 import Electric
#from RCAIDE.Framework.Networks.Thermal_Management.All_Electric_Thermal_Management_Network           import All_Electric_Thermal_Management_Network
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor      import design_propeller 
from RCAIDE.Library.Methods.Performance.estimate_stall_speed        import estimate_stall_speed 
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor   import design_motor 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion import nasa_motor
from RCAIDE.Library.Methods.Energy.Sources.Batteries.Common           import initialize_from_circuit_configuration
from RCAIDE.Library.Methods.Geometry.Planform       import wing_segmented_planform 
from RCAIDE.Library.Methods.Weights.Physics_Based_Buildups.Electric import compute_weight , converge_weight 
from RCAIDE.Library.Plots                                           import *  
from RCAIDE.Library.Methods.Thermal_Management.Heat_Exchangers.Cross_Flow_Heat_Exchanger        import design_cross_flow_heat_exchanger
from RCAIDE.Library.Methods.Thermal_Management.Batteries.Liquid_Cooled_Wavy_Channel  import design_wavy_channel    

# python imports 
import numpy as np 
from copy import deepcopy
import os 
import matplotlib.pyplot        as plt 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():

    BTMS_flag = True
    
    # vehicle data
    vehicle  = vehicle_setup(BTMS_flag)
    #plot_3d_vehicle(vehicle)
    
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
    
    return 
 
    
def vehicle_setup(BTMS_flag): 
    

    #------------------------------------------------------------------------------------------------------------------------------------
    #   Initialize the Vehicle
    #------------------------------------------------------------------------------------------------------------------------------------

    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Twin_Otter'

 
    # ################################################# Vehicle-level Properties ########################################################  

    # mass properties
    vehicle.mass_properties.max_takeoff   = 5670  # kg 
    vehicle.mass_properties.takeoff       = 5670  # kg 
    vehicle.mass_properties.max_zero_fuel = 5670  # kg 
    vehicle.envelope.ultimate_load        = 5.7
    vehicle.envelope.limit_load           = 3.8 
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
    vehicle.design_mach_number            =  mach_number

         
    # ##########################################################  Wings ################################################################    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Main Wing
    #------------------------------------------------------------------------------------------------------------------------------------
    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing' 
    wing.sweeps.quarter_chord             = 0.0 * Units.deg
    wing.thickness_to_chord               = 0.12
    wing.areas.reference                  = 39 
    wing.spans.projected                  = 19.81
    wing.chords.root                      = 2.03 
    wing.chords.tip                       = 2.03 
    wing.chords.mean_aerodynamic          = 2.03 
    wing.taper                            = wing.chords.root/wing.chords.tip 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 3. * Units.degree 
    wing.twists.tip                       = 0
    wing.origin                           = [[5.38, 0, 1, 35]] 
    wing.aerodynamic_center               = [[5.38 + 0.25 *wing.chords.root , 0, 1, 35]]  
    wing.vertical                         = False
    wing.symmetric                        = True
    wing.high_lift                        = True 
    wing.winglet_fraction                 = 0.0  
    wing.dynamic_pressure_ratio           = 1.0  
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator + '..' + separator + '..' + separator 
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
    segment.dihedral_outboard             = 3. * Units.degree 
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = 0.12
    segment.append_airfoil(airfoil)
    wing.append_segment(segment)
    
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0
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
    wing.sweeps.quarter_chord             = 0.0 * Units.degree
    wing.thickness_to_chord               = 0.12 
    wing.areas.reference                  = 9.762 
    wing.spans.projected                  = 6.29   
    wing.chords.root                      = 1.552 
    wing.chords.tip                       = 1.552 
    wing.chords.mean_aerodynamic          = 1.552  
    wing.taper                            = 0 
    wing.aspect_ratio                     = wing.spans.projected**2. / wing.areas.reference 
    wing.twists.root                      = 0.0 * Units.degree
    wing.twists.tip                       = 0.0 * Units.degree 
    wing.origin                           = [[13.17 , 0 , 1.25]] 
    wing.aerodynamic_center               = [[13.17 , 0 , 1.25]]  
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
    wing.sweeps.leading_edge              = 28.6 * Units.degree 
    wing.thickness_to_chord               = 0.12 
    wing.areas.reference                  = 8.753 
    wing.spans.projected                  = 3.9 
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

    # add to vehicle
    vehicle.append_component(wing)

 
    # ##########################################################   Fuselage  ############################################################    
    fuselage = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.seats_abreast                      = 2.
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
 
    #########################################################   Nacelles  ############################################################    
    nacelle                    = RCAIDE.Library.Components.Nacelles.Stack_Nacelle()
    nacelle.tag                = 'nacelle_1'
    nacelle.length             = 2
    nacelle.diameter           = 0.73480616 
    nacelle.areas.wetted       = 0.01*(2*np.pi*0.01/2)
    nacelle.origin             = [[2.81,3.34 , 1.22]]
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
    nac_segment.percent_z_location = 1.2284
    nac_segment.height             = 0.183333333 
    nac_segment.width              = 0.422484315 
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_3'
    nac_segment.percent_x_location = 0.143080714 
    nac_segment.percent_z_location = 1.246733333
    nac_segment.height             = 0.44	 
    nac_segment.width              = 0.685705173 
    nacelle.append_segment(nac_segment)  
     
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.170379029  
    nac_segment.percent_z_location = 1.054233333
    nac_segment.height             = 0.898333333	 
    nac_segment.width              = 0.73480616 
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.252189893  
    nac_segment.percent_z_location = 1.054233333
    nac_segment.height             = 1.008333333 
    nac_segment.width              = 0.736964445
    nacelle.append_segment(nac_segment)   
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.383860821   
    nac_segment.percent_z_location = 1.072566667
    nac_segment.height             = 0.971666667 
    nac_segment.width              = 0.736964445 
    nacelle.append_segment(nac_segment)  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 0.551826736  
    nac_segment.percent_z_location = 1.155066667	
    nac_segment.height             = 0.77	 
    nac_segment.width              = 0.736964445  
    nacelle.append_segment(nac_segment)
    

    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_8'
    nac_segment.percent_x_location = 0.809871485   
    nac_segment.percent_z_location = 1.2284
    nac_segment.height             = 0.366666667 
    nac_segment.width              = 0.736964445 
    nacelle.append_segment(nac_segment) 
    

    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_9'
    nac_segment.percent_x_location = 1.0  
    nac_segment.percent_z_location = 1.301733333 
    nac_segment.height             = 0.036666667	
    nac_segment.width              = 0.0  
    nacelle.append_segment(nac_segment)
    
    
    vehicle.append_component(nacelle)  

    nacelle_2          = deepcopy(nacelle)
    nacelle_2.tag      = 'nacelle_2'
    nacelle_2.origin   = [[ 2.81, -3.34 ,1.22]]
    vehicle.append_component(nacelle_2)    

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
    net                              = RCAIDE.Framework.Networks.Electric()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                              = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus()
    

    #------------------------------------------------------------------------------------------------------------------------------------           
    # Battery
    ##------------------------------------------------------------------------------------------------------------------------------------      
    bat_module                                                    = RCAIDE.Library.Components.Energy.Sources.Battery_Modules.Lithium_Ion_NMC()
    bat_module.electrical_configuration.series               = 10
    bat_module.electrical_configuration.parallel             = 210
    bat_module.cell.nominal_capacity                              = 3.8
    initialize_from_circuit_configuration(bat_module,module_weight_factor = 1.25)  
   
    #bat_module.module.power_split_ratio                               = 0.5
    bat_module.geometrtic_configuration.total              = bat_module.electrical_configuration.parallel*bat_module.electrical_configuration.series   #bat_module.pack.electrical_configuration.total
    bat_module.voltage                                     = bat_module.maximum_voltage #/bat_module.pack.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
    bat_module.geometrtic_configuration.normal_count       = 42
    bat_module.geometrtic_configuration.parallel_count     = 50
    #bat_module.capacity_Ah                                   = bat_module.cell.nominal_capacity * bat_module.electrical_configuration.parallel

    
    bat_module.nominal_capacity                   = bat_module.cell.nominal_capacity* bat_module.electrical_configuration.parallel

    for _ in range(12):
        bat_copy = deepcopy(bat_module)
        bus.battery_modules.append(bat_copy)
      

        
    bus.battery_module_electric_configuration =  'Series'
    bus.charging_c_rate                                          = 1

    bus.nominal_capacity = 0
    for battery_module in  bus.battery_modules:
        bus.voltage  +=   battery_module.voltage
        bus.nominal_capacity =  max(battery_module.nominal_capacity, bus.nominal_capacity)        
    
    
    #bus.charging_current                   = bus.nominal_capacity * bus.Charging_C_Rate 
    #bus.charging_time  =   bus.capacity_Ah /(bus.nominal_capacity * bus.Charging_C_Rate) *Units.hrs     


    
    #bat1                                                    = RCAIDE.Library.Components.Energy.Sources.Battery_Modules.Lithium_Ion_LFP()
    ##bat1.pack.electrical_configuration.series               = 120  
    ##bat1.pack.electrical_configuration.parallel             = 210   #250
    #bat1.cell.nominal_capacity                              = 3.8  
    ##initialize_from_circuit_configuration(bat1,module_weight_factor = 1.25)  
    #bat1.pack.number_of_modules                           = 12
    #bat1.module.geometrtic_configuration.total              = bat1.pack.electrical_configuration.total
    #bat1.module.voltage                                     = bat1.pack.maximum_voltage/bat1.pack.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
    #bat1.module.geometrtic_configuration.normal_count       = bat1.module.geometrtic_configuration.total/bat1.pack.number_of_modules / 50
    #bat1.module.geometrtic_configuration.parallel_count     = 50
    #bus.voltage                                            = bat1.pack.maximum_voltage
    
    #bus.batteries.append(bat1)
        
    

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    starboard_propulsor                              = RCAIDE.Library.Components.Propulsors.Electric_Rotor()  
    starboard_propulsor.tag                          = 'starboard_propulsor'
    starboard_propulsor.active_busses                = ['bus']   
  
    # Electronic Speed Controller       
    esc                                              = RCAIDE.Library.Components.Energy.Modulators.Electronic_Speed_Controller()
    esc.tag                                          = 'esc_1'
    esc.efficiency                                   = 0.95 
    esc.origin                                       = [[3.8,2.8129,1.22 ]]   
    starboard_propulsor.electronic_speed_controller  = esc   
     
    # Propeller              
    propeller                                        = RCAIDE.Library.Components.Propulsors.Converters.Propeller() 
    propeller.tag                                    = 'propeller_1'  
    propeller.tip_radius                             = 2.59
    propeller.number_of_blades                       = 3
    propeller.hub_radius                             = 10.    * Units.inches

    propeller.cruise.design_freestream_velocity      = 130 * Units.kts      
    speed_of_sound                                   = 343 
    propeller.cruise.design_tip_mach                 = 0.65
    propeller.cruise.design_angular_velocity         = propeller.cruise.design_tip_mach *speed_of_sound/propeller.tip_radius
    propeller.cruise.design_Cl                       = 0.7
    propeller.cruise.design_altitude                 = 8000. * Units.feet 
    propeller.cruise.design_thrust                   = 12500  
    propeller.clockwise_rotation                     = False
    propeller.variable_pitch                         = True  
    propeller.origin                                 = [[3.5,2.8129,1.22 ]]   
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
    motor                                            = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    motor.efficiency                                 = 0.98
    motor.origin                                     = [[4.0,2.8129,1.22 ]]
    
    
    motor.nominal_voltage                            = bus.voltage
    motor.no_load_current                            = 1
    motor.rotor_radius                               = propeller.tip_radius  
    motor.design_torque                              = propeller.cruise.design_torque 
    motor.angular_velocity                           = propeller.cruise.design_angular_velocity # Horse power of gas engine variant  750 * Units['hp']
    design_motor(motor)  
    motor.mass_properties.mass                       = nasa_motor(motor.design_torque) 
    starboard_propulsor.motor                        = motor 
 
    # append propulsor to distribution line 
    bus.propulsors.append(starboard_propulsor) 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Port Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    port_propulsor                             = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    port_propulsor.tag                         = "port_propulsor"
    port_propulsor.active_busses               = ['bus']   
            
    esc_2                                      = deepcopy(esc)
    esc_2.origin                               = [[3.8, -2.8129,1.22 ]]        
    port_propulsor.electronic_speed_controller = esc_2  

    propeller_2                                = deepcopy(propeller)
    propeller_2.tag                            = 'propeller_2' 
    propeller_2.origin                         =  [[3.5, -2.8129,1.22 ]]   
    propeller_2.clockwise_rotation             = False 
    port_propulsor.rotor                       = propeller_2  
              
    motor_2                                    = deepcopy(motor)
    motor_2.origin                             =  [[4.0, -2.8129,1.22 ]]        
    port_propulsor.motor                       = motor_2  
    
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
    

    ##------------------------------------------------------------------------------------------------------------------------------------  
    # Coolant Line
    #------------------------------------------------------------------------------------------------------------------------------------  
    coolant_line                      = RCAIDE.Library.Components.Energy.Distributors.Coolant_Line(bus)
    #coolant_line_1                      = RCAIDE.Library.Components.Energy.Distributors.Coolant_Line(bus)
    net.coolant_lines.append(coolant_line)
    #net.coolant_lines.append(coolant_line_1)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
     #Battery Thermal Management 
    #------------------------------------------------------------------------------------------------------------------------------------      
    #HAS                              = RCAIDE.Library.Components.Thermal_Management.Batteries.Air_Cooled()
    HAS                                                    = RCAIDE.Library.Components.Thermal_Management.Batteries.Liquid_Cooled_Wavy_Channel(coolant_line)
    #HAS2                                                    = RCAIDE.Library.Components.Thermal_Management.Batteries.Liquid_Cooled_Wavy_Channel(coolant_line) 
    HAS.design_altitude                                    = 2500. * Units.feet  
    atmosphere                                             = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976() 
    atmo_data                                              = atmosphere.compute_values(altitude = HAS.design_altitude)     
    HAS.coolant_inlet_temperature                          = atmo_data.temperature[0,0]  
    HAS.design_battery_operating_temperature               = 313
    HAS.design_heat_removed                                = 50000 /12 
    HAS                                                    = design_wavy_channel(HAS,bat_module) 
    
    
   
    for battery_module in bus.battery_modules:
        coolant_line.battery_modules[battery_module.tag].append(HAS)
    
                                            
    # Battery Heat Exchanger               
    HEX                                                    = RCAIDE.Library.Components.Thermal_Management.Heat_Exchangers.Cross_Flow_Heat_Exchanger() 
    HEX.design_altitude                                    = 2500. * Units.feet 
    HEX.inlet_temperature_of_cold_fluid                    = atmo_data.temperature[0,0]   
    HEX                                                    = design_cross_flow_heat_exchanger(HEX,coolant_line,bat_module)    
    coolant_line.heat_exchangers.append(HEX)
    
    
    # Reservoir for Battery TMS
    RES                                                    = RCAIDE.Library.Components.Thermal_Management.Reservoirs.Reservoir()

    coolant_line.reservoirs.append(RES)

    
    
    ##has1.tag =  'FOR nmc'
    ##has2                               = RCAIDE.Library.Components.Thermal_Management.Batteries.Air_Cooled()
    ##has2.tag =  'FOR lfp'    
    

    vehicle.append_energy_network(net)   
    
     
         
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################   Determine Vehicle Mass Properties Using Physic Based Methods  ################################ 
    #------------------------------------------------------------------------------------------------------------------------------------   
    #converge_weight(vehicle) 
    #breakdown = compute_weight(vehicle)
    #print(breakdown)  
     
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

    
    #config                                                = RCAIDE.Library.Components.Configs.Config(vehicle) 
    #config.tag                                            = 'no_hex_operation'
    #config_tms                                            = config.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.thermal_management_system
    #config_tms.heat_exchanger_system.percent_operation    = 0 
    #config_tms.heat_acquisition_system.percent_operation  = 0
    #config_tms.heat_exchanger_system.fan_operation        = False
    #configs.append(config)  
     
    
    #config                                                = RCAIDE.Library.Components.Configs.Config(vehicle) 
    #config.tag                                            = 'max_hex_operation'
    #config_tms                                            = config.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.thermal_management_system
    #config_tms.heat_exchanger_system.percent_operation    = 1 
    #config_tms.heat_acquisition_system.percent_operation  = 1
    #config_tms.heat_exchanger_system.fan_operation        = True
    #configs.append(config)  
     
 
    #config                                                = RCAIDE.Library.Components.Configs.Config(vehicle)   
    #config.tag                                            = 'hex_low_alt_climb_operation'
    #config_tms                                            = config.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.thermal_management_system
    #config_tms.heat_exchanger_system.percent_operation    = 0.5
    #config_tms.heat_acquisition_system.percent_operation  = 0.5
    #config_tms.heat_exchanger_system.fan_operation        = True
    #configs.append(config)     

    #config                                                = RCAIDE.Library.Components.Configs.Config(vehicle)
    #config.tag                                            = 'hex_high_alt_climb_operation'
    #config_tms                                            = config.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.thermal_management_system
    #config_tms.heat_exchanger_system.percent_operation    = 0.2
    #config_tms.heat_acquisition_system.percent_operation  = 0.5
    #config_tms.heat_exchanger_system.fan_operation        = True
    #configs.append(config)        

    #config                                                = RCAIDE.Library.Components.Configs.Config(vehicle) 
    #config.tag                                            = 'hex_cruise_operation'
    #config_tms                                            = config.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.thermal_management_system
    #config_tms.heat_exchanger_system.percent_operation    = 0.1
    #config_tms.heat_acquisition_system.percent_operation  = 0.8
    #config_tms.heat_exchanger_system.fan_operation        = False
    #configs.append(config)         
    

    #config                                               = RCAIDE.Library.Components.Configs.Config(vehicle) 
    #config.tag                                           = 'hex_descent_operation'
    #config_tms                                           = config.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.thermal_management_system
    #config_tms.heat_exchanger_system.percent_operation   = 0.1
    #config_tms.heat_acquisition_system.percent_operation = 0.8
    #config_tms.heat_exchanger_system.fan_operation       = False
    #configs.append(config)                   
 

    #config                                               = RCAIDE.Library.Components.Configs.Config(vehicle)
    #config.tag                                           = 'recharge'
    #config_tms                                           = config.networks.all_electric.busses.bus.batteries.lithium_ion_nmc.thermal_management_system
    #config_tms.heat_exchanger_system.percent_operation   = 1.0
    #config_tms.heat_acquisition_system.percent_operation = 1.0
    #config_tms.heat_exchanger_system.fan_operation       = True
    #configs.append(config)  
 
    return configs


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
    base_segment.temperature_deviation  = 2.5
    base_segment.state.numerics.number_of_control_points  = 16
  

    # VSTALL Calculation  
    vehicle        = analyses.base.aerodynamics.geometry
    vehicle_mass   = vehicle.mass_properties.max_takeoff
    reference_area = vehicle.reference_area 
    Vstall         = estimate_stall_speed(vehicle_mass,reference_area,altitude = 0.0,maximum_lift_coefficient = 1.2)
    
   
  
    

    ## ------------------------------------------------------------------
    ##   Taxi 
    ## ------------------------------------------------------------------      
    #segment = Segments.Ground.Idle(base_segment)
    #segment.tag = "Idle"  
    #segment.analyses.extend(analyses.no_hex_operation)  
    #segment.time                               = 0.5 * Units.hr
    #segment.initial_battery_state_of_charge    = 1.0  
    #mission.append_segment(segment)
     
    ## ------------------------------------------------------------------
    ##   Taxi 
    ## ------------------------------------------------------------------      
    #segment = Segments.Ground.Taxi(base_segment)
    #segment.tag = "Taxi"  
    #segment.analyses.extend( analyses.max_hex_operation)  
    #segment.velocity                                      = 25 * Units.knots  

    #segment.flight_dynamics.force_x                       = True   
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']]
    
    #mission.append_segment(segment)  
    
    ## ------------------------------------------------------------------
    ##   Takeoff
    ## ------------------------------------------------------------------      
    #segment = Segments.Ground.Takeoff(base_segment)
    #segment.tag = "Takeoff"  
    #segment.analyses.extend( analyses.max_hex_operation ) 
    #segment.velocity_end                                     = Vstall*1.2  
    #segment.friction_coefficient                             = 0.04   
    #segment.throttle                                         = 0.8   
    
    #segment.flight_dynamics.force_x                           = True 
    #segment.assigned_control_variables.elapsed_time.active               = True         
   
    #mission.append_segment(segment) 
  
    
    # ------------------------------------------------------------------
    #   Departure End of Runway Segment Flight 1 : 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Departure_End_of_Runway'       
    segment.analyses.extend( analyses.base )  
    #segment.analyses.extend( analyses.max_hex_operation )  
    segment.altitude_start                                = 0.0 * Units.feet
    segment.altitude_end                                  = 5
    0.0 * Units.feet
    segment.air_speed_start                               = Vstall *1.2  
    segment.air_speed_end                                 = Vstall *1.25
    segment.initial_battery_state_of_charge    = 1.0
            
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                  
       
    mission.append_segment(segment)
    
    # ------------------------------------------------------------------
    #   Initial Climb Area Segment Flight 1  
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Initial_CLimb_Area' 
    segment.analyses.extend( analyses.base )  
    #segment.analyses.extend( analyses.max_hex_operation )   
    segment.altitude_start                                = 50.0 * Units.feet
    segment.altitude_end                                  = 500.0 * Units.feet 
    segment.air_speed_end                                 = Vstall *1.3 
    segment.climb_rate                                    = 600 * Units['ft/min']   
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                  
          
    mission.append_segment(segment)  
   
             
    # ------------------------------------------------------------------
    #   Climb Segment Flight 1 
    # ------------------------------------------------------------------ 
    segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    segment.tag = 'Climb_1'        
    segment.analyses.extend( analyses.base )  
    #segment.analyses.extend( analyses.hex_low_alt_climb_operation )      
    segment.altitude_start                                = 500.0 * Units.feet
    segment.altitude_end                                  = 2500 * Units.feet  
    segment.air_speed_end                                 = 120 * Units.kts 
    segment.climb_rate                                    = 500* Units['ft/min']  
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    segment.assigned_control_variables.body_angle.active             = True                 
           
    mission.append_segment(segment)
    
        
    ## ------------------------------------------------------------------
    ##   Climb 1 : constant Speed, constant rate segment 
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    #segment.tag = "Climb_2"
    #segment.analyses.extend( analyses.base )  
    ##segment.analyses.extend( analyses.hex_high_alt_climb_operation)
    #segment.altitude_start                                = 2500.0  * Units.feet
    #segment.altitude_end                                  = 5000   * Units.feet  
    #segment.air_speed_end                                 = 130 * Units.kts 
    #segment.climb_rate                                    = 700.034 * Units['ft/min']   
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                 
            
    #mission.append_segment(segment)

    ## ------------------------------------------------------------------
    ##   Cruise Segment: constant Speed, constant altitude
    ## ------------------------------------------------------------------ 
    #segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    #segment.tag = "Cruise" 
    #segment.analyses.extend( analyses.base )  
    ##segment.analyses.extend(analyses.hex_cruise_operation) 
    #segment.altitude                                      = 5000   * Units.feet 
    #segment.air_speed                                     = 130 * Units.kts
    #segment.distance                                      = 20.   * Units.nautical_mile  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                  
          
    #mission.append_segment(segment)    


    ## ------------------------------------------------------------------
    ##   Descent Segment Flight 1   
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    #segment.tag = "Decent"  
    #segment.analyses.extend( analyses.base )  
    ##segment.analyses.extend( analyses.hex_descent_operation )       
    #segment.altitude_start                                = 5000   * Units.feet 
    #segment.altitude_end                                  = 1000 * Units.feet  
    #segment.air_speed_end                                 = 100 * Units['mph']   
    #segment.climb_rate                                    = -200 * Units['ft/min']  
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                 
          
    #mission.append_segment(segment)   
               
    ## ------------------------------------------------------------------
    ##  Downleg_Altitude Segment Flight 1 
    ## ------------------------------------------------------------------

    #segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    #segment.tag = 'Downleg'
    #segment.analyses.extend( analyses.base )  
    ##segment.analyses.extend(analyses.hex_descent_operation)  
    #segment.air_speed                                     = 100 * Units['mph']   
    #segment.distance                                      = 6000 * Units.feet 
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                   
            
    #mission.append_segment(segment)     
    
    ### ------------------------------------------------------------------
    ###  Reserve Climb 
    ### ------------------------------------------------------------------ 
    ##segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment) 
    ##segment.tag = 'Reserve_Climb'        
    ##segment.analyses.extend( analyses.hex_low_alt_climb_operation)      
    ##segment.altitude_end                                  = 5000 * Units.feet
    ##segment.air_speed                                     = 120 * Units['mph']
    ##segment.climb_rate                                    = 500* Units['ft/min']  
    
    ### define flight dynamics to model 
    ##segment.flight_dynamics.force_x                       = True  
    ##segment.flight_dynamics.force_z                       = True     
    
    ### define flight controls 
    ##segment.assigned_control_variables.throttle.active               = True           
    ##segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    ##segment.assigned_control_variables.body_angle.active             = True                
        
    ##mission.append_segment(segment)
    
    ### ------------------------------------------------------------------
    ###  Researve Cruise Segment 
    ### ------------------------------------------------------------------ 
    ##segment = Segments.Cruise.Constant_Speed_Constant_Altitude_Loiter(base_segment) 
    ##segment.tag = 'Reserve_Cruise'  
    ##segment.analyses.extend(analyses.hex_cruise_operation)  
    ##segment.altitude                                      = 5000 * Units.feet
    ##segment.air_speed                                     = 130 * Units.kts
    ##segment.time                                          = 60*30 * Units.sec  
    
    ### define flight dynamics to model 
    ##segment.flight_dynamics.force_x                       = True  
    ##segment.flight_dynamics.force_z                       = True     
    
    ### define flight controls 
    ##segment.assigned_control_variables.throttle.active               = True           
    ##segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    ##segment.assigned_control_variables.body_angle.active             = True                  
       
    ##mission.append_segment(segment)     
    
    ### ------------------------------------------------------------------
    ###  Researve Descent
    ### ------------------------------------------------------------------ 
    ##segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment) 
    ##segment.tag = 'Reserve_Descent'
    ##segment.analyses.extend( analyses.hex_descent_operation)    
    ##segment.altitude_end                                  = 1000 * Units.feet 
    ##segment.air_speed                                     = 110 * Units['mph']
    ##segment.descent_rate                                  = 300 * Units['ft/min']   
    
    ### define flight dynamics to model 
    ##segment.flight_dynamics.force_x                       = True  
    ##segment.flight_dynamics.force_z                       = True     
    
    ### define flight controls 
    ##segment.assigned_control_variables.throttle.active               = True           
    ##segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    ##segment.assigned_control_variables.body_angle.active             = True                
    ##mission.append_segment(segment)  

    
    ## ------------------------------------------------------------------
    ##  Baseleg Segment Flight 1  
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment)
    #segment.tag = 'Baseleg'
    #segment.analyses.extend( analyses.base )  
    ##segment.analyses.extend( analyses.hex_descent_operation)   
    #segment.altitude_start                                = 1000 * Units.feet
    #segment.altitude_end                                  = 500.0 * Units.feet
    #segment.air_speed_end                                 = 90 * Units['mph']  
    #segment.climb_rate                                    = -350 * Units['ft/min'] 
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    #mission.append_segment(segment) 

    ## ------------------------------------------------------------------
    ##  Final Approach Segment Flight 1  
    ## ------------------------------------------------------------------ 
    #segment = Segments.Climb.Linear_Speed_Constant_Rate(base_segment) 
    #segment.tag = 'Final_Approach'
    #segment.analyses.extend( analyses.base )      
    ##segment.analyses.extend( analyses.hex_descent_operation)      
    #segment.altitude_start                                = 500.0 * Units.feet
    #segment.altitude_end                                  = 00.0 * Units.feet
    #segment.air_speed_end                                 = 80 * Units['mph']  
    #segment.climb_rate                                    = -300 * Units['ft/min']   
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_propulsor','port_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                      
    #mission.append_segment(segment)  


    ## ------------------------------------------------------------------
    ##   Landing  
    ## ------------------------------------------------------------------  
    #segment = Segments.Ground.Landing(base_segment)
    #segment.tag = "Landing"   
    #segment.analyses.extend( analyses.base )
    ##segment.analyses.extend( analyses.hex_descent_operation)  
    #segment.velocity_end                                     = Vstall*0.1   
    #segment.assigned_control_variables.elapsed_time.active           = True  
    #segment.assigned_control_variables.elapsed_time.initial_guess_values  = [[30.]]  
    
    #mission.append_segment(segment)
    
    ## ------------------------------------------------------------------
    ##  Charge Segment: 
    ## ------------------------------------------------------------------     
    ## Charge Model 
    #segment                               = Segments.Ground.Battery_Recharge(base_segment)     
    #segment.tag  = 'Charge_Day'
    #segment.analyses.extend( analyses.base )
    ##segment.state.numerics.number_of_control_points  = 64 
    ##segment.analyses.extend(analyses.recharge)
    ##segment.time =  vehicle.networks.electric.busses.bus.charging_time
    ##segment.current = vehicle.networks.electric.busses.bus.charging_current
    
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


# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results):
    
    plot_propulsor_throttles(results)
    
    plot_flight_conditions(results) 
    
    plot_aerodynamic_forces(results)

    plot_aerodynamic_coefficients(results)  
    
    plot_aircraft_velocities(results)
    
    plot_bus_conditions(results)
    
    plot_battery_cell_conditions(results)
    
    plot_thermal_management_component(results)

    plot_battery_degradation(results)

    plot_rotor_conditions(results) 

    plot_electric_propulsor_efficiencies(results) 
    
    plot_battery_pack_C_rates(results)

    return


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
    
    # Calculate extra drag from landing gear: 
    main_wheel_width  = 4. * Units.inches
    main_wheel_height = 12. * Units.inches
    nose_gear_height  = 10. * Units.inches
    nose_gear_width   = 4. * Units.inches 
    total_wheel       = 2*main_wheel_width*main_wheel_height + nose_gear_width*nose_gear_height 
    main_gear_strut_height = 2. * Units.inches
    main_gear_strut_length = 24. * Units.inches
    nose_gear_strut_height = 12. * Units.inches
    nose_gear_strut_width  = 2. * Units.inches 
    total_strut = 2*main_gear_strut_height*main_gear_strut_length + nose_gear_strut_height*nose_gear_strut_width 
    drag_area = 1.4*( total_wheel + total_strut)
    
    
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.geometry                            = vehicle
    aerodynamics.settings.drag_coefficient_increment = drag_area/vehicle.reference_area
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle  = vehicle 
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


def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses



# ----------------------------------------------------------------------        
#   Call Main
# ----------------------------------------------------------------------    

if __name__ == '__main__':
    main()
    plt.show()
