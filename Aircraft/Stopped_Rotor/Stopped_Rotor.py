''' 
# Stopped_Rotor_EVTOL.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units, Data   
from RCAIDE.Framework.Networks.All_Electric_Network                            import All_Electric_Network 
from RCAIDE.Library.Methods.Performance.estimate_cruise_drag                   import estimate_cruise_drag
from RCAIDE.Library.Methods.Geometry.Planform                                  import segment_properties,wing_segmented_planform    
from RCAIDE.Library.Methods.Energy.Sources.Battery.Common                      import initialize_from_circuit_configuration 
from RCAIDE.Library.Methods.Weights.Correlation_Buildups.Propulsion            import nasa_motor
from RCAIDE.Library.Methods.Propulsors.Converters.DC_Motor                     import design_motor
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor                        import design_propeller ,design_lift_rotor 
from RCAIDE.Library.Methods.Weights.Physics_Based_Buildups.Electric            import compute_weight , converge_weight 
from RCAIDE.Library.Plots                                                      import *       
 
import os
import numpy as np 
from copy import deepcopy
import matplotlib.pyplot as plt 
import  pickle
# ----------------------------------------------------------------------------------------------------------------------
#  REGRESSION
# ----------------------------------------------------------------------------------------------------------------------  
def main():           
         
    # vehicle data
    new_geometry =  True 
    if new_geometry :
        vehicle  = vehicle_setup()
        save_aircraft_geometry(vehicle , 'Stopped_Rotor')
    else: 
        vehicle = load_aircraft_geometry('Stopped_Rotor')

    # plot vehicle 
    plot_3d_vehicle(vehicle, 
                    min_x_axis_limit            = -5,
                    max_x_axis_limit            = 15,
                    min_y_axis_limit            = -10,
                    max_y_axis_limit            = 10,
                    min_z_axis_limit            = -10,
                    max_z_axis_limit            = 10,
                    show_figure                 = False 
                    )           

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
     
    return
 
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
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics          = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)   

    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Framework.Analyses.Energy.Energy()
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

    # done!
    return analyses    



# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup() :
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------  
     
         
    vehicle               = RCAIDE.Vehicle()
    vehicle.tag           = 'Stopped_Rotor'
    vehicle.configuration = 'eVTOL'
     
    #------------------------------------------------------------------------------------------------------------------------------------
    # ################################################# Vehicle-level Properties #####################################################  
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # mass properties 
    vehicle.mass_properties.max_takeoff       = 2700 
    vehicle.mass_properties.takeoff           = vehicle.mass_properties.max_takeoff
    vehicle.mass_properties.operating_empty   = vehicle.mass_properties.max_takeoff
    vehicle.envelope.ultimate_load            = 5.7   
    vehicle.envelope.limit_load               = 3.  
    vehicle.passengers                        = 5 
        
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################  Wings ################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------ 
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Main Wing
    #------------------------------------------------------------------------------------------------------------------------------------
    wing                          = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                      = 'main_wing'  
    wing.aspect_ratio             = 8.95198  # will  be overwritten
    wing.sweeps.quarter_chord     = 0.0  
    wing.thickness_to_chord       = 0.14 
    wing.taper                    = 0.292
    wing.spans.projected          = 11.82855
    wing.chords.root              = 1.75
    wing.total_length             = 1.75
    wing.chords.tip               = 1.0
    wing.chords.mean_aerodynamic  = 0.8
    wing.dihedral                 = 0.0  
    wing.areas.reference          = 15.629
    wing.twists.root              = 4. * Units.degrees
    wing.twists.tip               = 0. 
    wing.origin                   = [[1.5, 0., 0.991]]
    wing.aerodynamic_center       = [ 1.567, 0., 0.991]    
    wing.winglet_fraction         = 0.0  
    wing.symmetric                = True
    wing.vertical                 = False
    
    ospath                        = os.path.abspath(__file__) 
    separator                     = os.path.sep
    rel_path                      = os.path.dirname(ospath) + separator + '..' +  separator
    airfoil                       = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.coordinate_file       = rel_path + 'Airfoils' + separator + 'NACA_63_412.txt'
    
    # Segment                                  
    segment                       = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'Section_1'   
    segment.percent_span_location = 0.0
    segment.twist                 = 4. * Units.degrees 
    segment.root_chord_percent    = 1. 
    segment.dihedral_outboard     = 8 * Units.degrees
    segment.sweeps.quarter_chord  = 0.9  * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)               
    
    # Segment                                   
    segment                       = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'Section_2'    
    segment.percent_span_location = 3.5/wing.spans.projected
    segment.twist                 = 3. * Units.degrees 
    segment.root_chord_percent    = 1.4000/1.7500
    segment.dihedral_outboard     = 0.0 * Units.degrees
    segment.sweeps.quarter_chord  = 1.27273 * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)               
     
    # Segment                                  
    segment                       = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'Section_3'   
    segment.percent_span_location = 11.3/wing.spans.projected 
    segment.twist                 = 2.0 * Units.degrees 
    segment.root_chord_percent    = 1.000/1.7500
    segment.dihedral_outboard     = 35.000* Units.degrees 
    segment.sweeps.quarter_chord  = 45.000* Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)     
    
    # Segment                                  
    segment                       = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'Section_4'   
    segment.percent_span_location = 11.6/wing.spans.projected 
    segment.twist                 = 0.0 * Units.degrees 
    segment.root_chord_percent    = 0.9/1.7500
    segment.dihedral_outboard     = 60. * Units.degrees 
    segment.sweeps.quarter_chord  = 70.0 * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)  
    
    # Segment                                  
    segment                       = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                   = 'Section_5'   
    segment.percent_span_location = 1.0
    segment.twist                 = 0.0 * Units.degrees 
    segment.root_chord_percent    = 0.35/1.7500
    segment.dihedral_outboard     = 0  * Units.degrees 
    segment.sweeps.quarter_chord  = 0  * Units.degrees 
    segment.thickness_to_chord    = 0.16  
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)                 
    
    
    # compute reference properties 
    wing_segmented_planform(wing, overwrite_reference = True ) 
    wing = segment_properties(wing)
    vehicle.reference_area        = wing.areas.reference  
    wing.areas.wetted             = wing.areas.reference  * 2 
    wing.areas.exposed            = wing.areas.reference  * 2  
        
    # add to vehicle 
    vehicle.append_component(wing)   
    
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    #   Horizontal Tail
    #------------------------------------------------------------------------------------------------------------------------------------
    wing                          = RCAIDE.Library.Components.Wings.Wing()
    wing.tag                      = 'horizontal_tail'  
    wing.aspect_ratio             = 3.04444
    wing.sweeps.quarter_chord     = 17. * Units.degrees
    wing.thickness_to_chord       = 0.12 
    wing.spans.projected          = 2.71805
    wing.chords.root              = 0.94940
    wing.total_length             = 0.94940
    wing.chords.tip               = 0.62731 
    wing.chords.mean_aerodynamic  = 0.809 
    wing.dihedral                 = 20 *Units.degrees
    wing.taper                    = wing.chords.tip / wing.chords.root 
    wing.areas.reference          = 2.14279
    wing.areas.wetted             = 2.14279   * 2
    wing.areas.exposed            = 2.14279   * 2
    wing.twists.root              = 0.0
    wing.twists.tip               = 0.0
    wing.origin                   = [[  5.374 ,0.0 ,  0.596]]
    wing.aerodynamic_center       = [   5.374, 0.0,   0.596] 
    wing.winglet_fraction         = 0.0 
    wing.symmetric                = True    
    
    # add to vehicle 
    vehicle.append_component(wing)     
      
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################   Fuselage  ############################################################   
    #------------------------------------------------------------------------------------------------------------------------------------ 
    fuselage                                    = RCAIDE.Library.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage' 
    fuselage.seats_abreast                      = 2.  
    fuselage.seat_pitch                         = 3.  
    fuselage.fineness.nose                      = 0.88   
    fuselage.fineness.tail                      = 1.13   
    fuselage.lengths.nose                       = 0.5  
    fuselage.lengths.tail                       = 1.5
    fuselage.lengths.cabin                      = 4.46 
    fuselage.lengths.total                      = 6.46
    fuselage.width                              = 5.85 * Units.feet      # change 
    fuselage.heights.maximum                    = 4.65 * Units.feet      # change 
    fuselage.heights.at_quarter_length          = 3.75 * Units.feet      # change 
    fuselage.heights.at_wing_root_quarter_chord = 4.65 * Units.feet      # change 
    fuselage.heights.at_three_quarters_length   = 4.26 * Units.feet      # change 
    fuselage.areas.wetted                       = 236. * Units.feet**2   # change 
    fuselage.areas.front_projected              = 0.14 * Units.feet**2   # change 
    fuselage.effective_diameter                 = 1.276     # change 
    fuselage.differential_pressure              = 0. 
    
    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_0'    
    segment.percent_x_location                  = 0.0 
    segment.percent_z_location                  = 0.     # change  
    segment.height                              = 0.049 
    segment.width                               = 0.032 
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                             
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'   
    segment.percent_x_location                  = 0.10912/fuselage.lengths.total 
    segment.percent_z_location                  = 0.00849
    segment.height                              = 0.481 
    segment.width                               = 0.553 
    fuselage.Segments.append(segment)           
                                                
    # Segment                                             
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'   
    segment.percent_x_location                  = 0.47804/fuselage.lengths.total
    segment.percent_z_location                  = 0.02874
    segment.height                              = 1.00
    segment.width                               = 0.912 
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                            
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'   
    segment.percent_x_location                  = 0.161  
    segment.percent_z_location                  = 0.04348  
    segment.height                              = 1.41
    segment.width                               = 1.174  
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                             
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'   
    segment.percent_x_location                  = 0.284 
    segment.percent_z_location                  = 0.05435 
    segment.height                              = 1.62
    segment.width                               = 1.276  
    fuselage.Segments.append(segment)              
                                                
    # Segment                                             
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'   
    segment.percent_x_location                  = 3.43026/fuselage.lengths.total
    segment.percent_z_location                  = 0.31483/fuselage.lengths.total 
    segment.height                              = 1.409
    segment.width                               = 1.121 
    fuselage.Segments.append(segment)                     
                                                
    # Segment                                             
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  = 4.20546/fuselage.lengths.total
    segment.percent_z_location                  = 0.32216/fuselage.lengths.total
    segment.height                              = 1.11
    segment.width                               = 0.833
    fuselage.Segments.append(segment)                  
                                                
    # Segment                                             
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'   
    segment.percent_x_location                  = 4.99358/fuselage.lengths.total
    segment.percent_z_location                  = 0.37815/fuselage.lengths.total
    segment.height                              = 0.78
    segment.width                               = 0.512 
    fuselage.Segments.append(segment)                  
                                                
    # Segment                                             
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_8'   
    segment.percent_x_location                  = 1.
    segment.percent_z_location                  = 0.55/fuselage.lengths.total
    segment.height                              = 0.195  
    segment.width                               = 0.130 
    fuselage.Segments.append(segment)                   
                                                
    vehicle.append_component(fuselage) 
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##########################################################  Booms  ################################################################  
    #------------------------------------------------------------------------------------------------------------------------------------          
    boom                                    = RCAIDE.Library.Components.Booms.Boom()
    boom.tag                                = 'boom_1r'
    boom.configuration                      = 'boom'  
    boom.origin                             = [[   0.036, 1.950,  1]]  
    boom.seats_abreast                      = 0.  
    boom.seat_pitch                         = 0.0 
    boom.fineness.nose                      = 0.950   
    boom.fineness.tail                      = 1.029   
    boom.lengths.nose                       = 0.2 
    boom.lengths.tail                       = 0.2
    boom.lengths.cabin                      = 4.15
    boom.lengths.total                      = 4.2
    boom.width                              = 0.15 
    boom.heights.maximum                    = 0.15  
    boom.heights.at_quarter_length          = 0.15  
    boom.heights.at_three_quarters_length   = 0.15 
    boom.heights.at_wing_root_quarter_chord = 0.15 
    boom.areas.wetted                       = 0.018
    boom.areas.front_projected              = 0.018 
    boom.effective_diameter                 = 0.15  
    boom.differential_pressure              = 0.  
    boom.symmetric                          = True 
    boom.index                              = 1
    
    # Segment  
    segment                           = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                       = 'segment_1'   
    segment.percent_x_location        = 0.
    segment.percent_z_location        = 0.0 
    segment.height                    = 0.05  
    segment.width                     = 0.05   
    boom.Segments.append(segment)           
    
    # Segment                                   
    segment                           = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                       = 'segment_2'   
    segment.percent_x_location        = 0.03
    segment.percent_z_location        = 0. 
    segment.height                    = 0.15 
    segment.width                     = 0.15 
    boom.Segments.append(segment) 
    
    # Segment                                   
    segment                           = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                       = 'segment_3'    
    segment.percent_x_location        = 0.97
    segment.percent_z_location        = 0. 
    segment.height                    = 0.15
    segment.width                     = 0.15
    boom.Segments.append(segment)           
    
    # Segment                                  
    segment                           = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                       = 'segment_4'   
    segment.percent_x_location        = 1.   
    segment.percent_z_location        = 0.   
    segment.height                    = 0.05   
    segment.width                     = 0.05   
    boom.Segments.append(segment)           
    
    # add to vehicle
    vehicle.append_component(boom)   
    
    # add left long boom 
    boom              = deepcopy(vehicle.booms.boom_1r)
    boom.origin[0][1] = -boom.origin[0][1]
    boom.tag          = 'boom_1l' 
    vehicle.append_component(boom)         
     
    # add left long boom 
    boom              = deepcopy(vehicle.booms.boom_1r)
    boom.origin       = [[     0.110,    4.891,   1.050]] 
    boom.tag          = 'boom_2r' 
    boom.lengths.total                      = 4.16
    vehicle.append_component(boom)  
     
    # add inner left boom 
    boom              = deepcopy(vehicle.booms.boom_1r)
    boom.origin       = [[     0.110, -  4.891,    1.050 ]]   
    boom.lengths.total                      = 4.16
    boom.tag          = 'boom_2l' 
    vehicle.append_component(boom)      
      
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################   Determine Vehicle Mass Properties Using Physic Based Methods  ################################ 
    #------------------------------------------------------------------------------------------------------------------------------------     
    sys                            = RCAIDE.Library.Components.Systems.System()
    sys.mass_properties.mass       = 5 # kg   
    vehicle.append_component(sys)    
   
    #------------------------------------------------------------------------------------------------------------------------------------
    # ########################################################  Energy Network  ######################################################### 
    #------------------------------------------------------------------------------------------------------------------------------------
    # define network
    network                                                = All_Electric_Network()   
    
    #==================================================================================================================================== 
    # Forward Bus
    #====================================================================================================================================  
    cruise_bus                                             = RCAIDE.Library.Components.Energy.Distribution.Electrical_Bus() 
    cruise_bus.tag                                         = 'cruise_bus'
     
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus Battery
    #------------------------------------------------------------------------------------------------------------------------------------ 
    bat                                                    = RCAIDE.Library.Components.Energy.Batteries.Lithium_Ion_NMC() 
    bat.tag                                                = 'cruise_bus_battery'
    bat.pack.electrical_configuration.series               = 140   
    bat.pack.electrical_configuration.parallel             = 60
    initialize_from_circuit_configuration(bat)  
    bat.module.number_of_modules                           = 14 
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module.voltage                                     = bat.pack.maximum_voltage/bat.module.number_of_modules 
    bat.module.geometrtic_configuration.normal_count       = 25
    bat.module.geometrtic_configuration.parallel_count     = 40 
    cruise_bus.voltage                                     =  bat.pack.maximum_voltage  
    cruise_bus.batteries.append(bat)      
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Forward Bus Propulsors  
    #------------------------------------------------------------------------------------------------------------------------------------       
    # Define Forward Propulsor Container 
    cruise_propulsor_1                                     = RCAIDE.Library.Components.Propulsors.Electric_Rotor()
    cruise_propulsor_1.tag                                 = 'cruise_propulsor_1' 
    cruise_propulsor_1.active_batteries                    = ['cruise_bus_battery']   
                 
    # Electronic Speed Controller                     
    propeller_esc                                          = RCAIDE.Library.Components.Propulsors.Modulators.Electronic_Speed_Controller() 
    propeller_esc.efficiency                               = 0.95  
    propeller_esc.origin                                   = [[6.583, 1.300,  1.092 ]] 
    propeller_esc.tag                                      = 'propeller_esc_1' 
    cruise_propulsor_1.electronic_speed_controller= propeller_esc      
    
    # Propeller 
    g                                                      = 9.81                                   # gravitational acceleration 
    speed_of_sound                                         = 340                                    # speed of sound  
    Hover_Load                                             = vehicle.mass_properties.takeoff*g      # hover load   
            
    propeller                                              = RCAIDE.Library.Components.Propulsors.Converters.Propeller()
    propeller.number_of_blades                             = 3
    propeller.tag                                          = 'propeller_1'  
    propeller.origin                                       = [[6.583, 1.300,  1.092 ]] 
    propeller.tip_radius                                   = 1.15  
    propeller.hub_radius                                   = 0.1 * propeller.tip_radius  
    propeller.cruise.design_freestream_velocity            = 130.* Units['mph'] 
    propeller.cruise.design_tip_mach                       = 0.65
    propeller.cruise.design_angular_velocity               = propeller.cruise.design_tip_mach *speed_of_sound/propeller.tip_radius
    propeller.cruise.design_Cl                             = 0.7
    propeller.cruise.design_altitude                       = 1500 * Units.feet
    propeller.cruise.design_thrust                         = 2500
    propeller.rotation                                     = 1
    propeller.variable_pitch                               = True  
    airfoil                                                = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.coordinate_file                                = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                                    = [rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
    propeller.append_airfoil(airfoil)                     
    propeller.airfoil_polar_stations                       = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  
    propeller                                              = design_propeller(propeller)   
    cruise_propulsor_1.rotor                               = propeller    
                
    # Propeller Motor              
    propeller_motor                                        = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    propeller_motor.efficiency                             = 0.95
    propeller_motor.tag                                    = 'propeller_motor_1'  
    propeller_motor.origin                                 = [[6.583, 1.300,  1.092 ]] 
    propeller_motor.nominal_voltage                        = bat.pack.maximum_voltage 
    propeller_motor.origin                                 = propeller.origin
    propeller_motor.propeller_radius                       = propeller.tip_radius 
    propeller_motor.no_load_current                        = 0.001
    propeller_motor.wing_mounted                           = True 
    propeller_motor.wing_tag                               = 'horizontal_tail'
    propeller_motor.rotor_radius                           = propeller.tip_radius
    propeller_motor.design_torque                          = propeller.cruise.design_torque
    propeller_motor.angular_velocity                       = propeller.cruise.design_angular_velocity/propeller_motor.gear_ratio  
    propeller_motor                                        = design_motor(propeller_motor)  
    propeller_motor.mass_properties.mass                   = nasa_motor(propeller_motor.design_torque)  
    cruise_propulsor_1.motor                               = propeller_motor 
      
    # rear propeller nacelle 
    propeller_nacelle                = RCAIDE.Library.Components.Nacelles.Stack_Nacelle()
    propeller_nacelle.tag            = 'propeller_nacelle'
    propeller_nacelle.length         = 1.24
    propeller_nacelle.diameter       = 0.4    
    propeller_nacelle.origin        = [[5.583,  1.300 ,    1.092]]
    propeller_nacelle.flow_through   = False  
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_1'
    nac_segment.percent_x_location = 0.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.10/propeller_nacelle.length
    nac_segment.height             = 0.2
    nac_segment.width              = 0.2
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_2'
    nac_segment.percent_x_location = 0.15 /propeller_nacelle.length 
    nac_segment.height             = 0.25
    nac_segment.width              = 0.25
    propeller_nacelle.append_segment(nac_segment)    
    
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_4'
    nac_segment.percent_x_location = 0.2/propeller_nacelle.length  
    nac_segment.height             = 0.3
    nac_segment.width              = 0.3
    propeller_nacelle.append_segment(nac_segment)    
    
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_5'
    nac_segment.percent_x_location = 0.25/propeller_nacelle.length  
    nac_segment.height             = 0.35
    nac_segment.width              = 0.35
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_6'
    nac_segment.percent_x_location = 0.5/propeller_nacelle.length 
    nac_segment.height             = 0.4
    nac_segment.width              = 0.4
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_7'
    nac_segment.percent_x_location = 0.75/propeller_nacelle.length
    nac_segment.height             = 0.35
    nac_segment.width              = 0.35
    propeller_nacelle.append_segment(nac_segment)        
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_8'
    nac_segment.percent_x_location = 0.98/propeller_nacelle.length
    nac_segment.height             = 0.3
    nac_segment.width              = 0.3
    propeller_nacelle.append_segment(nac_segment)    
    
    nac_segment                    = RCAIDE.Library.Components.Nacelles.Segment()
    nac_segment.tag                = 'segment_9'
    nac_segment.percent_x_location = 1.0  
    nac_segment.height             = 0.0
    nac_segment.width              = 0.0
    propeller_nacelle.append_segment(nac_segment)      
    cruise_propulsor_1.nacelle = propeller_nacelle   
    cruise_bus.propulsors.append(cruise_propulsor_1)
      
    # make and append copy of forward propulsor (efficient coding)    
    cruise_propulsor_2                             = deepcopy(cruise_propulsor_1)
    cruise_propulsor_2.tag                         = 'cruise_propulsor_2' 
    cruise_propulsor_2.rotor.origin                = [[6.583, -1.300,  1.092 ]]   
    propeller_nacelle_2                            = deepcopy(propeller_nacelle)
    propeller_nacelle_2.tag                        = 'propeller_nacelle_2' 
    propeller_nacelle_2.origin                     = [[5.583, - 1.300,     1.092]]
    cruise_propulsor_2.nacelle                     = propeller_nacelle_2
    cruise_bus.propulsors.append(cruise_propulsor_2)    
        
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Additional Bus Loads
    #------------------------------------------------------------------------------------------------------------------------------------     
    # Payload   
    payload                        = RCAIDE.Library.Components.Systems.Avionics()
    payload.power_draw             = 10. # Watts 
    payload.mass_properties.mass   = 1.0 * Units.kg
    cruise_bus.payload             = payload 
    
    # Avionics   
    avionics                       = RCAIDE.Library.Components.Systems.Avionics()
    avionics.power_draw            = 10. # Watts  
    avionics.mass_properties.mass  = 1.0 * Units.kg
    cruise_bus.avionics            = avionics    

    # append forward bus
    network.busses.append(cruise_bus)    
    
        
    #==================================================================================================================================== 
    # Lift Bus 
    #====================================================================================================================================          
    lift_bus                                               = RCAIDE.Library.Components.Energy.Distribution.Electrical_Bus()
    lift_bus.tag                                           = 'lift_bus' 

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus Battery
    #------------------------------------------------------------------------------------------------------------------------------------ 
    bat                                                    = RCAIDE.Library.Components.Energy.Batteries.Lithium_Ion_NMC() 
    bat.tag                                                = 'lift_bus_battery'
    bat.pack.electrical_configuration.series               = 140   
    bat.pack.electrical_configuration.parallel             = 20
    initialize_from_circuit_configuration(bat)  
    bat.module.number_of_modules                           = 14 
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module.voltage                                     = bat.pack.maximum_voltage/bat.module.number_of_modules 
    bat.module.geometrtic_configuration.normal_count       = 25
    bat.module.geometrtic_configuration.parallel_count     = 40 
    lift_bus.voltage                                       =  bat.pack.maximum_voltage  
    lift_bus.batteries.append(bat)      
    

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Lift Propulsors 
    #------------------------------------------------------------------------------------------------------------------------------------    
     
    # Define Lift Propulsor Container 
    lift_propulsor_1                                       = RCAIDE.Library.Components.Propulsors.Electric_Rotor()
    lift_propulsor_1.tag                                   = 'lift_propulsor_1'     
    lift_propulsor_1.active_batteries                      = ['lift_bus_battery']          
              
    # Electronic Speed Controller           
    lift_rotor_esc                                         = RCAIDE.Library.Components.Propulsors.Modulators.Electronic_Speed_Controller()
    lift_rotor_esc.efficiency                              = 0.95    
    lift_rotor_esc.tag                                     = 'lift_rotor_esc_1' 
    lift_rotor_esc.origin                                  = [[-0.073 ,  1.950 , 1.2]] 
    lift_propulsor_1.electronic_speed_controller           = lift_rotor_esc 
           
    # Lift Rotor Design              
    lift_rotor                                             = RCAIDE.Library.Components.Propulsors.Converters.Lift_Rotor()   
    lift_rotor.tag                                         = 'lift_rotor_1'  
    lift_rotor.origin                                      = [[-0.073 ,  1.950 , 1.2]] 
    lift_rotor.active                                      = True          
    lift_rotor.orientation_euler_angles                    = [10.0*Units.degrees,np.pi/2.,0.]   # vector of angles defining default orientation of rotor
    lift_rotor.tip_radius                                  = 2.8/2
    lift_rotor.hub_radius                                  = 0.1 
    lift_rotor.number_of_blades                            = 3     
    lift_rotor.hover.design_altitude                       = 40 * Units.feet  
    lift_rotor.hover.design_thrust                         = Hover_Load/8
    lift_rotor.hover.design_freestream_velocity            = np.sqrt(lift_rotor.hover.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
    lift_rotor.oei.design_altitude                         = 40 * Units.feet  
    lift_rotor.oei.design_thrust                           = Hover_Load/7  
    lift_rotor.oei.design_freestream_velocity              = np.sqrt(lift_rotor.oei.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
    airfoil                                                = RCAIDE.Library.Components.Airfoils.Airfoil()   
    airfoil.coordinate_file                                = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
    airfoil.polar_files                                    = [rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                             rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                              rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
    lift_rotor.append_airfoil(airfoil)                         
    lift_rotor.airfoil_polar_stations                      = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  
    lift_rotor                                             = design_lift_rotor(lift_rotor) 
    lift_propulsor_1.rotor                                 = lift_rotor      
    
    #------------------------------------------------------------------------------------------------------------------------------------               
    # Lift Rotor Motor  
    #------------------------------------------------------------------------------------------------------------------------------------    
    lift_rotor_motor                                       = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    lift_rotor_motor.efficiency                            = 0.9
    lift_rotor_motor.nominal_voltage                       = bat.pack.maximum_voltage*3/4  
    lift_rotor_motor.origin                                = [[-0.073 ,  1.950 , 1.2]]
    lift_rotor_motor.propeller_radius                      = lift_rotor.tip_radius
    lift_rotor_motor.tag                                   = 'lift_rotor_motor_1' 
    lift_rotor_motor.no_load_current                       = 0.01  
    lift_rotor_motor.wing_mounted                          = True 
    lift_rotor_motor.wing_tag                              = 'main_wing'
    lift_rotor_motor.rotor_radius                          = lift_rotor.tip_radius
    lift_rotor_motor.design_torque                         = lift_rotor.hover.design_torque
    lift_rotor_motor.angular_velocity                      = lift_rotor.hover.design_angular_velocity/lift_rotor_motor.gear_ratio  
    lift_rotor_motor                                       = design_motor(lift_rotor_motor)
    lift_rotor_motor.mass_properties.mass                  = nasa_motor(lift_rotor_motor.design_torque)     
    lift_propulsor_1.motor                                 = lift_rotor_motor
    


    #------------------------------------------------------------------------------------------------------------------------------------               
    # Lift Rotor Nacelle
    #------------------------------------------------------------------------------------------------------------------------------------     
    nacelle                           = RCAIDE.Library.Components.Nacelles.Nacelle()
    nacelle.tag                       = 'rotor_nacelle'
    nacelle.length                    = 0.45
    nacelle.diameter                  = 0.3
    nacelle.orientation_euler_angles  = [0,-90*Units.degrees,0.]    
    nacelle.flow_through              = False   
    nacelle.tag                       = 'rotor_nacelle_1'
    nacelle.origin                    = [[  -0.073,  1.950, 1.2]]
    lift_propulsor_1.nacelle          =  nacelle 
    lift_bus.propulsors.append(lift_propulsor_1)   


    # make and append copy of lift propulsor (efficient coding)    
    lift_propulsor_2                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_2.tag                                   = 'lift_propulsor_2' 
    lift_propulsor_2.rotor.origin                          = [[-0.073  , -1.950  , 1.2]]  
    lift_propulsor_2.rotor.orientation_euler_angle         = [-10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_2.motor.origin                          = [[-0.073  , -1.950  , 1.2]] 
    lift_propulsor_2.rotor.origin                          = [[-0.073  , -1.950  , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_2' 
    rotor_nacelle.origin                                   = [[ -0.073, -1.950, 1.2]]
    lift_propulsor_2.nacelle                               = rotor_nacelle  
    lift_bus.propulsors.append(lift_propulsor_2)    
        

    lift_propulsor_3                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_3.tag                                   = 'lift_propulsor_3' 
    lift_propulsor_3.rotor.origin                          = [[ 4.440 ,  1.950 , 1.2]]  
    lift_propulsor_3.rotor.orientation_euler_angle         = [10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_3.motor.origin                          = [[ 4.440 ,  1.950 , 1.2]] 
    lift_propulsor_3.rotor.origin                          = [[ 4.440 ,  1.950 , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_3'  
    rotor_nacelle.origin                                   = [[   4.413,   1.950 ,1.2]]
    lift_propulsor_3.nacelle                               = rotor_nacelle  
    lift_bus.propulsors.append(lift_propulsor_3)    
    

    lift_propulsor_4                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_4.tag                                   = 'lift_propulsor_4' 
    lift_propulsor_4.rotor.origin                          = [[ 4.440  , -1.950  , 1.2]]  
    lift_propulsor_4.rotor.orientation_euler_angle         = [-10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_4.motor.origin                          = [[ 4.440  , -1.950  , 1.2]] 
    lift_propulsor_4.rotor.origin                          = [[ 4.440  , -1.950  , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_4' 
    rotor_nacelle.origin                                   = [[   4.413, -1.950, 1.2]]
    lift_propulsor_4.nacelle                               = rotor_nacelle  
    vehicle.append_component(rotor_nacelle) 
    

    lift_propulsor_5                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_5.tag                                   = 'lift_propulsor_5' 
    lift_propulsor_5.rotor.origin                          = [[ 0.219 ,  4.891 , 1.2]]  
    lift_propulsor_5.rotor.orientation_euler_angle         = [10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_5.motor.origin                          = [[ 0.219 ,  4.891 , 1.2]] 
    lift_propulsor_5.rotor.origin                          = [[ 0.219 ,  4.891 , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_5'  
    rotor_nacelle.origin                                   = [[   0.219 ,   4.891 , 1.2]] 
    lift_propulsor_5.nacelle                               = rotor_nacelle   
    lift_bus.propulsors.append(lift_propulsor_5)    
    

    lift_propulsor_6                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_6.tag                                   = 'lift_propulsor_6' 
    lift_propulsor_6.rotor.origin                          = [[ 0.219  , - 4.891 , 1.2]]  
    lift_propulsor_6.rotor.orientation_euler_angle         = [-10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_6.motor.origin                          = [[ 0.219  , - 4.891 , 1.2]] 
    lift_propulsor_6.rotor.origin                          = [[ 0.219  , - 4.891 , 1.2]] 
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_6'  
    rotor_nacelle.origin                                   = [[   0.219 , -  4.891 ,1.2]]
    lift_propulsor_6.nacelle                               = rotor_nacelle   
    lift_bus.propulsors.append(lift_propulsor_6)    
    
    

    lift_propulsor_7                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_7.tag                                   = 'lift_propulsor_7' 
    lift_propulsor_7.rotor.origin                          = [[ 4.196 ,  4.891 , 1.2]]  
    lift_propulsor_7.rotor.orientation_euler_angle         = [10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_7.motor.origin                          = [[ 4.196 ,  4.891 , 1.2]] 
    lift_propulsor_7.rotor.origin                          = [[ 4.196 ,  4.891 , 1.2]]   
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_7'  
    rotor_nacelle.origin                                   = [[  4.196 ,   4.891 ,1.2]]
    lift_propulsor_7.nacelle                               = rotor_nacelle  
    lift_bus.propulsors.append(lift_propulsor_7)    
    

    lift_propulsor_8                                       = deepcopy(lift_propulsor_1)
    lift_propulsor_8.tag                                   = 'lift_propulsor_8' 
    lift_propulsor_8.rotor.origin                          = [[ 4.196  , - 4.891 , 1.2]]  
    lift_propulsor_8.rotor.orientation_euler_angle         = [-10.0* Units.degrees,np.pi/2.,0.]
    lift_propulsor_8.motor.origin                          = [[ 4.196  , - 4.891 , 1.2]] 
    lift_propulsor_8.rotor.origin                          = [[ 4.196  , - 4.891 , 1.2]]  
    rotor_nacelle                                          = deepcopy(nacelle)
    rotor_nacelle.tag                                      = 'rotor_nacelle_8' 
    rotor_nacelle.origin                                   = [[   4.196, -  4.891 ,1.2]]
    lift_propulsor_8.nacelle                               = rotor_nacelle  
    lift_bus.propulsors.append(lift_propulsor_8)        
     

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Additional Bus Loads
    #------------------------------------------------------------------------------------------------------------------------------------            
    # Payload   
    payload                                                 = RCAIDE.Library.Components.Systems.Avionics()
    payload.power_draw                                      = 10. # Watts 
    payload.mass_properties.mass                            = 1.0 * Units.kg
    lift_bus.payload                                        = payload 
                             
    # Avionics                            
    avionics                                                = RCAIDE.Library.Components.Systems.Avionics()
    avionics.power_draw                                     = 10. # Watts  
    avionics.mass_properties.mass                           = 1.0 * Units.kg
    lift_bus.avionics                                       = avionics    

   
    network.busses.append(lift_bus)       
        
    # append energy network 
    vehicle.append_energy_network(network) 
    
    #------------------------------------------------------------------------------------------------------------------------------------
    # ##################################   Determine Vehicle Mass Properties Using Physic Based Methods  ################################ 
    #------------------------------------------------------------------------------------------------------------------------------------   
    converge_weight(vehicle) 
    breakdown = compute_weight(vehicle)
    print(breakdown) 
     
    return vehicle


# ---------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle): 

    configs = RCAIDE.Library.Components.Configs.Config.Container()

    base_config                                                       = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag                                                   = 'base'     
    configs.append(base_config) 

    forward_config                                                    = RCAIDE.Library.Components.Configs.Config(vehicle)
    forward_config.tag                                                = 'forward_flight'  
    forward_config.networks.all_electric.busses['lift_bus'].active    = False  
    configs.append(forward_config)  

    transition_config                                                 = RCAIDE.Library.Components.Configs.Config(vehicle)
    transition_config.tag                                             = 'transition_flight'    
    configs.append(transition_config)
    

    vertical_config                                                   = RCAIDE.Library.Components.Configs.Config(vehicle)
    vertical_config.tag                                               = 'vertical_flight'  
    vertical_config.networks.all_electric.busses['cruise_bus'].active = False  
    configs.append(vertical_config)   
     
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
    
    # base segment
    base_segment = Segments.Segment()       
      
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'baseline_mission' 
    
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    # base segment           
    base_segment  = Segments.Segment()   
    
    # VSTALL Calculation   
    Vstall         = 48.3144   
     
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Vertical Climb 
    #------------------------------------------------------------------------------------------------------------------------------------  
    segment     = Segments.Vertical_Flight.Climb(base_segment)
    segment.tag = "Vertical_Climb"   
    segment.analyses.extend( analyses.vertical_flight )  
    segment.altitude_start                                = 0.0  * Units.ft  
    segment.altitude_end                                  = 200.  * Units.ft   
    segment.initial_battery_state_of_charge               = 1.0 
    segment.climb_rate                                    = 500. * Units['ft/min']   
            
    # define flight dynamics to model  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True           
    segment.flight_controls.throttle.assigned_propulsors  = [['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                              'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']] 
     
    segment.initial_battery_state_of_charge              = 0.89 
    segment.initial_battery_resistance_growth_factor     = 1
    segment.initial_battery_capacity_fade_factor         = 1             
       
    mission.append_segment(segment)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Low-Speed Transition
    #------------------------------------------------------------------------------------------------------------------------------------  
 
    segment                                               = Segments.Transition.Constant_Acceleration_Constant_Pitchrate_Constant_Altitude(base_segment)
    segment.tag                                           = "Low_Speed_Transition"  
    segment.analyses.extend( analyses.transition_flight ) 
              
    segment.air_speed_end                                 = 0.75 * Vstall
    segment.acceleration                                  = 1.5
    segment.pitch_initial                                 = 0.0 * Units.degrees
    segment.pitch_final                                   = 2.  * Units.degrees   

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.flight_controls.throttle.active               = True
    segment.flight_controls.throttle.initial_guess_values = [[1.0],[1.0]]
    segment.flight_controls.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2'],
                                                             ['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                            'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']]
    mission.append_segment(segment) 
    
    ##------------------------------------------------------------------------------------------------------------------------------------  
    ## High-Speed Climbing Transition 
    ##------------------------------------------------------------------------------------------------------------------------------------  
    #segment                                               = Segments.Transition.Constant_Acceleration_Constant_Angle_Linear_Climb(base_segment)
    #segment.tag                                           = "High_Speed_Climbing_Transition" 
    #segment.analyses.extend( analyses.transition_flight)    
    #segment.altitude_start                                = 200.0 * Units.ft   
    #segment.altitude_end                                  = 500.0 * Units.ft   
    #segment.air_speed_start                               = 0.75   * Vstall
    #segment.climb_angle                                   = 3     * Units.degrees   
    #segment.acceleration                                  = 0.25  * Units['m/s/s'] 
    #segment.pitch_initial                                 = 2.    * Units.degrees 
    #segment.pitch_final                                   = 7.    * Units.degrees   

    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.flight_controls.throttle.active               = True           
    #segment.flight_controls.throttle.assigned_propulsors  = [['cruise_propulsor_1','cruise_propulsor_2'],
                                                             #['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                            #'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']]
    #mission.append_segment(segment) 
   
    ##------------------------------------------------------------------------------------------------------------------------------------ 
    ## Vertical Descent 
    ##------------------------------------------------------------------------------------------------------------------------------------ 
    #segment                                               = Segments.Vertical_Flight.Descent(base_segment)
    #segment.tag                                           = "Vertical_Descent" 
    #segment.analyses.extend( analyses.vertical_flight)     
    #segment.altitude_start                                = 300.0 * Units.ft   
    #segment.altitude_end                                  = 0.   * Units.ft  
    #segment.descent_rate                                  = 300. * Units['ft/min']  
    
    ## define flight dynamics to model  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.flight_controls.throttle.active               = True           
    #segment.flight_controls.throttle.assigned_propulsors  = [['lift_propulsor_1','lift_propulsor_2','lift_propulsor_3','lift_propulsor_4',
                                                              #'lift_propulsor_5','lift_propulsor_6','lift_propulsor_7','lift_propulsor_8']] 
            
    #mission.append_segment(segment)  
  

    ## ------------------------------------------------------------------
    ##  Charge Segment: 
    ## ------------------------------------------------------------------  
    ## Charge Model 
    #segment                                         = Segments.Ground.Battery_Recharge(base_segment)     
    #segment.analyses.extend(analyses.base)              
    #segment.tag                                     = 'Recharge' 
    #segment.time                                    = 1 * Units.hr
    #segment.current                                 = 100  
    #mission.append_segment(segment)                

    return mission

def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions

def plot_results(results):

    # Plots fligh conditions 
    plot_flight_conditions(results) 
    
    # Plot arcraft trajectory
    plot_flight_trajectory(results)   

    plot_propulsor_throttles(results)
    
    # Plot Aircraft Electronics
    plot_battery_pack_conditions(results) 
    plot_battery_temperature(results)
    plot_battery_cell_conditions(results) 
    plot_battery_pack_C_rates(results)
    plot_battery_degradation(results) 
    
    # Plot Propeller Conditions 
    plot_rotor_conditions(results) 
    plot_disc_and_power_loading(results)
    
    # Plot Electric Motor and Propeller Efficiencies 
    plot_electric_propulsor_efficiencies(results)  
      
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