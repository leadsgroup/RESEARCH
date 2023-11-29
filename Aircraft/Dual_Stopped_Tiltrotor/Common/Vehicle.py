''' 
# Vehicle.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M. Clarke 

'''
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Core import Units, Data   
import pickle
from RCAIDE.Visualization.Performance.Aerodynamics.Vehicle import *  
from RCAIDE.Visualization.Performance.Mission import *    
from RCAIDE.Visualization.Performance.Energy.Battery import *   
from RCAIDE.Visualization.Performance.Noise import *  
from RCAIDE.Visualization.Geometry import *
from RCAIDE.Components.Energy.Networks.Battery_Electric_Rotor                 import Battery_Electric_Rotor
from RCAIDE.Methods.Performance.estimate_cruise_drag                          import estimate_cruise_drag
from RCAIDE.Methods.Power.Battery.Sizing                                      import initialize_from_mass 
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform                         import segment_properties
from RCAIDE.Methods.Power.Battery.Sizing                                      import initialize_from_circuit_configuration 
from RCAIDE.Methods.Weights.Correlations.Propulsion                           import nasa_motor
from RCAIDE.Methods.Propulsion.electric_motor_sizing                          import size_optimal_motor
from RCAIDE.Methods.Propulsion                                                import propeller_design ,lift_rotor_design , prop_rotor_design
from RCAIDE.Methods.Weights.Buildups.eVTOL.empty                              import empty
from RCAIDE.Methods.Center_of_Gravity.compute_component_centers_of_gravity    import compute_component_centers_of_gravity
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform.wing_segmented_planform import wing_segmented_planform 
from RCAIDE.Methods.Weights.Buildups.eVTOL.converge_evtol_weight              import converge_evtol_weight  
 
import os
import numpy as np 
from copy import deepcopy 

# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup(resize_aircraft,vehicle_name = 'Dual_Stopped_Tiltrotor_CRM') :
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------  
    
    if resize_aircraft:
         
        vehicle               = RCAIDE.Vehicle()
        vehicle.tag           = vehicle_name
        vehicle.configuration = 'eVTOL'
        
        # ------------------------------------------------------------------
        #   Vehicle-level Properties
        # ------------------------------------------------------------------    
        # mass properties 
        vehicle.mass_properties.max_takeoff       = 2200 
        vehicle.mass_properties.takeoff           = vehicle.mass_properties.max_takeoff
        vehicle.mass_properties.operating_empty   = vehicle.mass_properties.max_takeoff
        vehicle.envelope.ultimate_load            = 5.7   
        vehicle.envelope.limit_load               = 3.  
        vehicle.passengers                        = 6 
           
        # ------------------------------------------------------------------    
        # WINGS                                    
        # ------------------------------------------------------------------  
        wing                          = RCAIDE.Components.Wings.Main_Wing()
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
        rel_path                      = ospath.split( 'Dual_Stopped_Tiltrotor' + separator + 'Common')[0]
        
        airfoil                       = RCAIDE.Components.Airfoils.Airfoil()
        airfoil.coordinate_file       = rel_path + 'Airfoils' + separator + 'NACA_63_412.txt'
        
        # Segment                                  
        segment                       = RCAIDE.Components.Wings.Segment()
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
        segment                       = RCAIDE.Components.Wings.Segment()
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
        segment                       = RCAIDE.Components.Wings.Segment()
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
        segment                       = RCAIDE.Components.Wings.Segment()
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
        segment                       = RCAIDE.Components.Wings.Segment()
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
        
        
        # WING PROPERTIES 
        wing                          = RCAIDE.Components.Wings.Wing()
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
         
         
        # ------------------------------------------------------------------
        #   Vertical Stabilizer
        # ------------------------------------------------------------------
    
        wing = RCAIDE.Components.Wings.Vertical_Tail()
        wing.tag = 'vertical_stabilizer'
    
    
        wing.aspect_ratio             = 1.87749
        wing.sweeps.quarter_chord     = 40 * Units.degrees
        wing.thickness_to_chord       = 0.12
        wing.spans.projected          = 1.78701
        wing.chords.root              = 1.79192
        wing.total_length             = 1.79192
        wing.chords.tip               = 0.44578
        wing.taper                    = wing.chords.tip / wing.chords.root 
        wing.chords.mean_aerodynamic  = 0.809  
        wing.areas.reference          = 1.70089
        wing.areas.wetted             = 1.70089 * 2
        wing.areas.exposed            = 1.70089 * 2
        wing.dihedral                 = 45.0 * Units.degrees 
        wing.twists.root              = 0.0 * Units.degrees
        wing.twists.tip               = 0.0 * Units.degrees
        wing.origin                   = [[  4.613,0.0 , 0.596]]
        wing.aerodynamic_center       = [0, 0.0, 0]  
    
        wing.vertical                = True
        wing.symmetric               = False
        wing.t_tail                  = False
    
        wing.dynamic_pressure_ratio  = 1.0
     
        # add to vehicle
        vehicle.append_component(wing) 
         
        # ---------------------------------------------------------------   
        # FUSELAGE                
        # ---------------------------------------------------------------   
        # FUSELAGE PROPERTIES
        fuselage                                    = RCAIDE.Components.Fuselages.Fuselage()
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
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment() 
        segment.tag                                 = 'segment_0'    
        segment.percent_x_location                  = 0.0 
        segment.percent_z_location                  = 0.     # change  
        segment.height                              = 0.049 
        segment.width                               = 0.032 
        fuselage.Segments.append(segment)                     
                                                    
        # Segment                                             
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_1'   
        segment.percent_x_location                  = 0.10912/fuselage.lengths.total 
        segment.percent_z_location                  = 0.00849
        segment.height                              = 0.481 
        segment.width                               = 0.553 
        fuselage.Segments.append(segment)           
                                                    
        # Segment                                             
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_2'   
        segment.percent_x_location                  = 0.47804/fuselage.lengths.total
        segment.percent_z_location                  = 0.02874
        segment.height                              = 1.00
        segment.width                               = 0.912 
        fuselage.Segments.append(segment)                     
                                                    
        # Segment                                            
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_3'   
        segment.percent_x_location                  = 0.161  
        segment.percent_z_location                  = 0.04348  
        segment.height                              = 1.41
        segment.width                               = 1.174  
        fuselage.Segments.append(segment)                     
                                                    
        # Segment                                             
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_4'   
        segment.percent_x_location                  = 0.284 
        segment.percent_z_location                  = 0.05435 
        segment.height                              = 1.62
        segment.width                               = 1.276  
        fuselage.Segments.append(segment)              
                                                    
        # Segment                                             
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_5'   
        segment.percent_x_location                  = 3.43026/fuselage.lengths.total
        segment.percent_z_location                  = 0.31483/fuselage.lengths.total 
        segment.height                              = 1.409
        segment.width                               = 1.121 
        fuselage.Segments.append(segment)                     
                                                    
        # Segment                                             
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_6'   
        segment.percent_x_location                  = 4.20546/fuselage.lengths.total
        segment.percent_z_location                  = 0.32216/fuselage.lengths.total
        segment.height                              = 1.11
        segment.width                               = 0.833
        fuselage.Segments.append(segment)                  
                                                    
        # Segment                                             
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_7'   
        segment.percent_x_location                  = 4.99358/fuselage.lengths.total
        segment.percent_z_location                  = 0.37815/fuselage.lengths.total
        segment.height                              = 0.78
        segment.width                               = 0.512 
        fuselage.Segments.append(segment)                  
                                                    
        # Segment                                             
        segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                                 = 'segment_8'   
        segment.percent_x_location                  = 1.
        segment.percent_z_location                  = 0.55/fuselage.lengths.total
        segment.height                              = 0.195  
        segment.width                               = 0.130 
        fuselage.Segments.append(segment)                   
                                                    
        vehicle.append_component(fuselage) 
        
        #-------------------------------------------------------------------
        # INNER BOOMS   
        #-------------------------------------------------------------------   
        boom                                    = RCAIDE.Components.Fuselages.Fuselage()
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
        segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment() 
        segment.tag                       = 'segment_1'   
        segment.percent_x_location        = 0.
        segment.percent_z_location        = 0.0 
        segment.height                    = 0.05  
        segment.width                     = 0.05   
        boom.Segments.append(segment)           
        
        # Segment                                   
        segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                       = 'segment_2'   
        segment.percent_x_location        = 0.03
        segment.percent_z_location        = 0. 
        segment.height                    = 0.15 
        segment.width                     = 0.15 
        boom.Segments.append(segment) 
        
        # Segment                                   
        segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                       = 'segment_3'    
        segment.percent_x_location        = 0.97
        segment.percent_z_location        = 0. 
        segment.height                    = 0.15
        segment.width                     = 0.15
        boom.Segments.append(segment)           
        
        # Segment                                  
        segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
        segment.tag                       = 'segment_4'   
        segment.percent_x_location        = 1.   
        segment.percent_z_location        = 0.   
        segment.height                    = 0.05   
        segment.width                     = 0.05   
        boom.Segments.append(segment)           
        
        # add to vehicle
        vehicle.append_component(boom)   
        
        # add left long boom 
        boom              = deepcopy(vehicle.fuselages.boom_1r)
        boom.origin[0][1] = -boom.origin[0][1]
        boom.tag          = 'boom_1l' 
        vehicle.append_component(boom)         
         
        # add left long boom 
        boom              = deepcopy(vehicle.fuselages.boom_1r)
        boom.origin       = [[     0.110,    4.891,   1.050]] 
        boom.tag          = 'boom_2r' 
        boom.lengths.total                      = 4.16
        vehicle.append_component(boom)  
         
        # add inner left boom 
        boom              = deepcopy(vehicle.fuselages.boom_1r)
        boom.origin       = [[     0.110, -  4.891,    1.050 ]]   
        boom.lengths.total                      = 4.16
        boom.tag          = 'boom_2l' 
        vehicle.append_component(boom)      
        
     
        # ------------------------------------------------------------------
        #   Nacelles
        # ------------------------------------------------------------------ 
        nacelle                = RCAIDE.Components.Nacelles.Nacelle()
        nacelle.tag            = 'rotor_nacelle'
        nacelle.length         = 0.45
        nacelle.diameter       = 0.3
        nacelle.orientation_euler_angles  = [0,-90*Units.degrees,0.]    
        nacelle.flow_through   = False  
        
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_1'
        nac_segment.percent_x_location = 0.0  
        nac_segment.height             = 0.0
        nac_segment.width              = 0.0
        nacelle.append_segment(nac_segment)    
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_2'
        nac_segment.percent_x_location = 0.05 
        nac_segment.height             = 0.1
        nac_segment.width              = 0.1
        nacelle.append_segment(nac_segment)    
        
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_3'
        nac_segment.percent_x_location = 0.15 
        nac_segment.height             = 0.2
        nac_segment.width              = 0.2
        nacelle.append_segment(nac_segment)        
    
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_4'
        nac_segment.percent_x_location = 0.25  
        nac_segment.height             = 0.25
        nac_segment.width              = 0.25
        nacelle.append_segment(nac_segment)     
        
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_5'
        nac_segment.percent_x_location = 0.25  
        nac_segment.height             = 0.25
        nac_segment.width              = 0.25
        nacelle.append_segment(nac_segment)    
        
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_6'
        nac_segment.percent_x_location = 0.5 
        nac_segment.height             = 0.3
        nac_segment.width              = 0.3
        nacelle.append_segment(nac_segment)    
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_7'
        nac_segment.percent_x_location = 0.75
        nac_segment.height             = 0.25
        nac_segment.width              = 0.25
        nacelle.append_segment(nac_segment)        
    
        nac_segment                    = RCAIDE.Components.Lofted_Body_Segment.Segment()
        nac_segment.tag                = 'segment_8'
        nac_segment.percent_x_location = 1.0
        nac_segment.height             = 0.2
        nac_segment.width              = 0.2
        nacelle.append_segment(nac_segment)      
     
        prop_rotor_nacelle_origins   = [[  -0.073,  1.950, 0.850] ,[ -0.073, -1.950, 0.850], 
                               [   0.219 ,   4.891 , 0.950] ,[   0.219 , -  4.891 ,0.950]]
     
        lift_rotor_nacelle_origins   = [ [   4.413,   1.950 ,0.850] ,[   4.413, -1.950, 0.850], 
                               [  4.196 ,   4.891 ,0.950] ,[   4.196, -  4.891 ,0.950]] 
     
        for ii in range(4):
            rotor_nacelle          = deepcopy(nacelle)
            rotor_nacelle.tag      = 'prop_rotor_nacelle_' + str(ii+1) 
            rotor_nacelle.origin   = [prop_rotor_nacelle_origins[ii]]
            vehicle.append_component(rotor_nacelle)       
    
        for ii in range(4):
            rotor_nacelle          = deepcopy(nacelle)
            rotor_nacelle.tag      = 'lift_rotor_nacelle_' + str(ii+1) 
            rotor_nacelle.origin   = [lift_rotor_nacelle_origins[ii]]
            vehicle.append_component(rotor_nacelle)   
                
        #------------------------------------------------------------------
        # network
        #------------------------------------------------------------------
        net                              = Battery_Electric_Rotor()
        net.rotor_group_indexes          = [0,0,0,0,1,1,1,1]
        net.motor_group_indexes          = [0,0,0,0,1,1,1,1]  
        net.esc_group_indexes            = [0,0,0,0,1,1,1,1]     
        net.active_propulsor_groups      = [True,True]
    
        #------------------------------------------------------------------
        # Design Battery
        #------------------------------------------------------------------    
        bat                                                    = RCAIDE.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion_LiNiMnCoO2_18650() 
        bat.pack.electrical_configuration.series               = 140   
        bat.pack.electrical_configuration.parallel             = 100
        initialize_from_circuit_configuration(bat)  
        bat.module_config.number_of_modules                    = 14 
        bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
        bat.module_config.voltage                              = bat.pack.max_voltage/bat.module_config.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
        bat.module.geometrtic_configuration.normal_count       = 25
        bat.module.geometrtic_configuration.parallel_count     = 40
        net.battery                                            = bat    
        net.voltage                                            = bat.pack.max_voltage
    
    
        # --------------------------------------------------------------
        # Lift/Forward Cruise Propulsor System 
        # -------------------------------------------------------------- 
        g               = 9.81                                   # gravitational acceleration  
        Drag            = estimate_cruise_drag(vehicle,altitude = 1500. * Units.ft,speed= 130.* Units['mph'] ,lift_coefficient = 0.5 ,profile_drag = 0.06)
        Hover_Load      = vehicle.mass_properties.takeoff*g      # hover load  
        
        # 1. Electronic Speed Controller  
        prop_rotor_esc              = RCAIDE.Components.Energy.Distributors.Electronic_Speed_Controller()
        prop_rotor_esc.efficiency   = 0.95
        for i in range(4):
            prop_rotor_ESC          = deepcopy(prop_rotor_esc)
            prop_rotor_ESC.tag      = 'prop_rotor_esc' + str(i + 1)  
            net.electronic_speed_controllers.append(prop_rotor_ESC) 
        
        # 2. Lift Rotors   
        prop_rotor                                   = RCAIDE.Components.Energy.Converters.Prop_Rotor()  
        prop_rotor.tip_radius                        = 2.8/2
        prop_rotor.hub_radius                        = 0.15 * prop_rotor.tip_radius   
        prop_rotor.number_of_blades                  = 4      
        
        prop_rotor.hover.design_altitude             = 40 * Units.feet 
        prop_rotor.hover.design_thrust               = Hover_Load/8
        prop_rotor.hover.design_freestream_velocity  = np.sqrt(prop_rotor.hover.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2)))    
        
        prop_rotor.cruise.design_altitude             = 1500 * Units.feet  
        prop_rotor.cruise.design_thrust               = Drag*1.1/8 
        prop_rotor.cruise.design_freestream_velocity  = 130*Units['mph']   
        
        prop_rotor.oei.design_altitude               = 40 * Units.feet  
        prop_rotor.oei.design_thrust                 = Hover_Load/6  
        prop_rotor.oei.design_freestream_velocity    = np.sqrt(prop_rotor.oei.design_thrust/(2*1.2*np.pi*(prop_rotor.tip_radius**2))) 
        
        airfoil                                = RCAIDE.Components.Airfoils.Airfoil()
        ospath    = os.path.abspath(__file__)
        separator = os.path.sep
        rel_path  = os.path.dirname(ospath) + separator     
        airfoil.coordinate_file                     = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
        airfoil.polar_files                         = [rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                      rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
        prop_rotor.append_airfoil(airfoil)          
        prop_rotor.airfoil_polar_stations            = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]    
        prop_rotor.variable_pitch                    = True    
        prop_rotor                                   = prop_rotor_design(prop_rotor)   
          
        prop_rotor_origins                          = [[ -0.073, 1.950,  1.2] ,  [  -0.073, -1.950 ,  1.2], 
                                                       [ 0.219 ,  4.891 , 1.2]  ,[ 0.219 , - 4.891 , 1.2], ]
        
        # Appending rotors with different origins
        rotations     = [-1,1,-1,1] 
        angle_offsets = np.random.rand(8)*(np.pi)    
        for i in range(4):
            pr                        = deepcopy(prop_rotor)
            pr.tag                    = 'prop_rotor_' + str(i)
            pr.rotation               = rotations[i]
            pr.origin                 = [prop_rotor_origins[i]]
            pr.phase_offset_angle     = angle_offsets[i]
            net.rotors.append(pr)    
            
             
        # 3. Lift Rotor Motors 
        prop_rotor_motor                         = RCAIDE.Components.Energy.Converters.Motor()
        prop_rotor_motor.efficiency              = 0.9
        prop_rotor_motor.nominal_voltage         = bat.pack.max_voltage*3/4  
        prop_rotor_motor.origin                  = prop_rotor.origin 
        prop_rotor_motor.propeller_radius        = prop_rotor.tip_radius   
        prop_rotor_motor.no_load_current         = 0.01  
        prop_rotor_motor.rotor_radius            = prop_rotor.tip_radius
        prop_rotor_motor.design_torque           = prop_rotor.hover.design_torque
        prop_rotor_motor.angular_velocity        = prop_rotor.hover.design_angular_velocity/prop_rotor_motor.gear_ratio  
        prop_rotor_motor                         = size_optimal_motor(prop_rotor_motor)
        prop_rotor_motor.mass_properties.mass    = nasa_motor(prop_rotor_motor.design_torque)    
    
        # Appending motors with different origins
        for i in range(4):
            pr_motor        = deepcopy(prop_rotor_motor)
            pr_motor.tag    = 'prop_rotor_motor_' + str(i+1)
            pr_motor.origin = [prop_rotor_origins[i]] 
            net.motors.append(pr_motor) 
    
    
        # --------------------------------------------------------------
        # Lift Propulsor System 
        # -------------------------------------------------------------- 
        # 1. Electronic Speed Controller  
        lift_rotor_esc              = RCAIDE.Components.Energy.Distributors.Electronic_Speed_Controller()
        lift_rotor_esc.efficiency   = 0.95
        for i in range(4):
            lift_rotor_ESC          = deepcopy(lift_rotor_esc)
            lift_rotor_ESC.tag      = 'lift_rotor_esc' + str(i)  
            net.electronic_speed_controllers.append(lift_rotor_ESC) 
        
        # 2. Lift Rotors   
        lift_rotor                                   = RCAIDE.Components.Energy.Converters.Lift_Rotor()  
        lift_rotor.tip_radius                        = 2.8/2
        lift_rotor.hub_radius                        = 0.15 * lift_rotor.tip_radius   
        lift_rotor.number_of_blades                  = 3     
        lift_rotor.hover.design_altitude             = 40 * Units.feet  
        lift_rotor.hover.design_thrust               = Hover_Load/8
        lift_rotor.hover.design_freestream_velocity  = np.sqrt(lift_rotor.hover.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
        lift_rotor.oei.design_altitude               = 40 * Units.feet  
        lift_rotor.oei.design_thrust                 = Hover_Load/6  
        lift_rotor.oei.design_freestream_velocity    = np.sqrt(lift_rotor.oei.design_thrust/(2*1.2*np.pi*(lift_rotor.tip_radius**2)))  
        airfoil                                      = RCAIDE.Components.Airfoils.Airfoil()   
        airfoil.coordinate_file                      = rel_path + 'Airfoils' + separator + 'NACA_4412.txt'
        airfoil.polar_files                          =[rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt' ,
                                                       rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt' ,
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt' ,
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt' ,
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_3500000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_5000000.txt',
                                                        rel_path + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_7500000.txt' ]
        lift_rotor.append_airfoil(airfoil)               
        lift_rotor.airfoil_polar_stations                = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        lift_rotor                                       = lift_rotor_design(lift_rotor)   
          
        lift_rotor_origins                           = [[ 4.440 , 1.950 ,  1.2] ,[ 4.440 , -1.950,  1.2],
                                                       [ 4.196 ,  4.891 , 1.2] ,[4.196, - 4.891 , 1.2]]
        
        # Appending rotors with different origins
        rotations = [-1,1,-1,1] 
        angle_offsets        = np.random.rand(8)*(np.pi)    
        for i in range(4):
            lr                        = deepcopy(lift_rotor)
            lr.tag                    = 'lift_rotor_' + str(i)
            lr.rotation               = rotations[i]
            lr.origin                 = [lift_rotor_origins[i]]
            lr.phase_offset_angle     = angle_offsets[i]
            net.rotors.append(lr)    
            
             
        # 3. Lift Rotor Motors 
        lift_rotor_motor                         = RCAIDE.Components.Energy.Converters.Motor()
        lift_rotor_motor.efficiency              = 0.9
        lift_rotor_motor.nominal_voltage         = bat.pack.max_voltage*3/4  
        lift_rotor_motor.origin                  = lift_rotor.origin 
        lift_rotor_motor.propeller_radius        = lift_rotor.tip_radius   
        lift_rotor_motor.no_load_current         = 0.01  
        lift_rotor_motor.rotor_radius            = lift_rotor.tip_radius
        lift_rotor_motor.design_torque           = lift_rotor.hover.design_torque
        lift_rotor_motor.angular_velocity        = lift_rotor.hover.design_angular_velocity/lift_rotor_motor.gear_ratio  
        lift_rotor_motor                         = size_optimal_motor(lift_rotor_motor)
        lift_rotor_motor.mass_properties.mass    = nasa_motor(lift_rotor_motor.design_torque)    
    
        # Appending motors with different origins
        for i in range(4):
            lr_motor = deepcopy(lift_rotor_motor)
            lr_motor.tag    = 'lift_rotor_motor_' + str(i+1)
            lr_motor.origin = [lift_rotor_origins[i]] 
            net.motors.append(lr_motor)   
    
        #------------------------------------------------------------------
        # Design Payload
        #------------------------------------------------------------------
        payload                        = RCAIDE.Components.Energy.Peripherals.Avionics()
        payload.power_draw             = 10. # Watts 
        payload.mass_properties.mass   = 1.0 * Units.kg
        net.payload                    = payload
    
        #------------------------------------------------------------------
        # Design Avionics
        #------------------------------------------------------------------
        avionics                       = RCAIDE.Components.Energy.Peripherals.Avionics()
        avionics.power_draw            = 20. # Watts  
        net.avionics                   = avionics  
    
        #------------------------------------------------------------------
        # Miscellaneous Systems 
        #------------------------------------------------------------------ 
        sys                            = RCAIDE.Components.Systems.System()
        sys.mass_properties.mass       = 5 # kg      
        
           
        # append motor locations to wing 
        rotor_motor_origins                                  = np.array(lift_rotor_origins + prop_rotor_origins ) 
        vehicle.wings['main_wing'].motor_spanwise_locations  = rotor_motor_origins[:,1]/vehicle.wings['main_wing'].spans.projected 
        
        # append motor origin spanwise locations onto wing data structure
        vehicle.append_component(net)
    
        settings = Data()
        converge_evtol_weight(vehicle,settings,contingency_factor = 1.0) 
        breakdown = empty(vehicle,settings,contingency_factor = 1.0 )
        print(breakdown)
        
        vehicle.weight_breakdown  = breakdown
        compute_component_centers_of_gravity(vehicle)
        vehicle.center_of_gravity() 
        
        save_aircraft_geometry(vehicle,vehicle.tag)        
        
    else: 
        vehicle = load_aircraft_geometry(vehicle_name) 
    
    return vehicle 
 
# ---------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle): 

    configs = RCAIDE.Components.Configs.Config.Container() 

    config              = RCAIDE.Components.Configs.Config(vehicle)
    config.tag          = 'vertical_flight' 
    vector_angle        = 90.0 * Units.degrees 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.orientation_euler_angles  = [0,vector_angle,0] 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.pitch_command    = 0
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.pitch_command    = 0
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.pitch_command    = 0
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.pitch_command    = 0
    config.networks.battery_electric_rotor.active_propulsor_groups                      = [True,True]# [prop_rotor,lift_rotor]
    configs.append(config)   


    config              = RCAIDE.Components.Configs.Config(vehicle)
    config.tag          = 'low_speed_transition_flight'
    vector_angle        = 60.0 * Units.degrees 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.orientation_euler_angles  = [0,vector_angle,0] 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.cruise.design_collective_pitch*1/3 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.cruise.design_collective_pitch*1/3 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.cruise.design_collective_pitch*1/3 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.cruise.design_collective_pitch*1/3 
    config.networks.battery_electric_rotor.active_propulsor_groups                      = [True,True]# [prop_rotor,lift_rotor]
    configs.append(config)   


    config              = RCAIDE.Components.Configs.Config(vehicle)
    config.tag          = 'high_speed_transition_flight'
    vector_angle        = 30.0 * Units.degrees 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.orientation_euler_angles  = [0,vector_angle,0] 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.cruise.design_collective_pitch*2/3
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.cruise.design_collective_pitch*2/3
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.cruise.design_collective_pitch*2/3
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.cruise.design_collective_pitch*2/3
    config.networks.battery_electric_rotor.active_propulsor_groups                      = [True,False]# [prop_rotor,lift_rotor]
    configs.append(config)  

    

    config              = RCAIDE.Components.Configs.Config(vehicle)
    config.tag          = 'forward_flight' 
    vector_angle        = 0.0 * Units.degrees 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.orientation_euler_angles  = [0,vector_angle,0]
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.y_axis_rotation  = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.orientation_euler_angles  = [0,vector_angle,0] 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_0.cruise.design_collective_pitch 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_1.cruise.design_collective_pitch 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_2.cruise.design_collective_pitch 
    vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.inputs.pitch_command    = vehicle.networks.battery_electric_rotor.rotors.prop_rotor_3.cruise.design_collective_pitch 
    config.networks.battery_electric_rotor.active_propulsor_groups                      = [True,False]# [prop_rotor,lift_rotor]
    configs.append(config)  
 
    return configs

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