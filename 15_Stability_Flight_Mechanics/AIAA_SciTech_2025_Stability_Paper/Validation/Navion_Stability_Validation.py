# vlm_pertubation_test.py
# 
# Created: May 2024, M. Clarke
 
# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
import RCAIDE 
from RCAIDE.Framework.Core import Units     
from RCAIDE.Library.Plots       import *  
from RCAIDE.Library.Methods.Geometry.Planform  import segment_properties
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor import design_propeller

# python imports  
import pylab as plt
import numpy as np 
import pandas as pd


# local imports 
import sys 
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
    
    #'''
    #Save Data
    #'''
    #df = pd.DataFrame({
    #'Clift_alpha'    : results.segments.cruise.conditions.static_stability.derivatives.Clift_alpha[0], 
    #'Clift_beta'     : results.segments.cruise.conditions.static_stability.derivatives.Clift_beta[0], 
    #'Clift_delta_a'  : results.segments.cruise.conditions.static_stability.derivatives.Clift_delta_a[0], 
    #'Clift_delta_e'  : results.segments.cruise.conditions.static_stability.derivatives.Clift_delta_e[0], 
    #'Clift_delta_r'  : results.segments.cruise.conditions.static_stability.derivatives.Clift_delta_r[0], 
    #'Clift_delta_f'  : results.segments.cruise.conditions.static_stability.derivatives.Clift_delta_f[0], 
    #'Clift_delta_s'  : results.segments.cruise.conditions.static_stability.derivatives.Clift_delta_s[0], 
    #'Cdrag_alpha'    : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_alpha[0], 
    #'Cdrag_beta'     : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_beta[0],  
    #'Cdrag_delta_a'  : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_delta_a[0], 
    #'Cdrag_delta_e'  : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_delta_e[0], 
    #'Cdrag_delta_r'  : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_delta_r[0], 
    #'Cdrag_delta_f'  : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_delta_f[0], 
    #'Cdrag_delta_s'  : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_delta_s[0], 
    #'CX_alpha'       : results.segments.cruise.conditions.static_stability.derivatives.CX_alpha[0], 
    #'CX_beta'        : results.segments.cruise.conditions.static_stability.derivatives.CX_beta[0], 
    #'CX_delta_a'     : results.segments.cruise.conditions.static_stability.derivatives.CX_delta_a[0], 
    #'CX_delta_e'     : results.segments.cruise.conditions.static_stability.derivatives.CX_delta_e[0], 
    #'CX_delta_r'     : results.segments.cruise.conditions.static_stability.derivatives.CX_delta_r[0], 
    #'CX_delta_f'     : results.segments.cruise.conditions.static_stability.derivatives.CX_delta_f[0], 
    #'CX_delta_s'     : results.segments.cruise.conditions.static_stability.derivatives.CX_delta_s[0], 
    #'CY_alpha'       : results.segments.cruise.conditions.static_stability.derivatives.CY_alpha[0], 
    #'CY_beta'        : results.segments.cruise.conditions.static_stability.derivatives.CY_beta[0], 
    #'CY_delta_a'     : results.segments.cruise.conditions.static_stability.derivatives.CY_delta_a[0], 
    #'CY_delta_e'     : results.segments.cruise.conditions.static_stability.derivatives.CY_delta_e[0], 
    #'CY_delta_r'     : results.segments.cruise.conditions.static_stability.derivatives.CY_delta_r[0], 
    #'CY_delta_f'     : results.segments.cruise.conditions.static_stability.derivatives.CY_delta_f[0], 
    #'CY_delta_s'     : results.segments.cruise.conditions.static_stability.derivatives.CY_delta_s[0], 
    #'CZ_alpha'       : results.segments.cruise.conditions.static_stability.derivatives.CZ_alpha[0], 
    #'CZ_beta'        : results.segments.cruise.conditions.static_stability.derivatives.CZ_beta[0], 
    #'CZ_delta_a'     : results.segments.cruise.conditions.static_stability.derivatives.CZ_delta_a[0], 
    #'CZ_delta_e'     : results.segments.cruise.conditions.static_stability.derivatives.CZ_delta_e[0], 
    #'CZ_delta_r'     : results.segments.cruise.conditions.static_stability.derivatives.CZ_delta_r[0], 
    #'CZ_delta_f'     : results.segments.cruise.conditions.static_stability.derivatives.CZ_delta_f[0], 
    #'CZ_delta_s'     : results.segments.cruise.conditions.static_stability.derivatives.CZ_delta_s[0], 
    #'CL_alpha'       : results.segments.cruise.conditions.static_stability.derivatives.CL_alpha[0], 
    #'CL_beta'        : results.segments.cruise.conditions.static_stability.derivatives.CL_beta[0], 
    #'CL_delta_a'     : results.segments.cruise.conditions.static_stability.derivatives.CL_delta_a[0], 
    #'CL_delta_e'     : results.segments.cruise.conditions.static_stability.derivatives.CL_delta_e[0], 
    #'CL_delta_r'     : results.segments.cruise.conditions.static_stability.derivatives.CL_delta_r[0], 
    #'CL_delta_f'     : results.segments.cruise.conditions.static_stability.derivatives.CL_delta_f[0], 
    #'CL_delta_s'     : results.segments.cruise.conditions.static_stability.derivatives.CL_delta_s[0], 
    #'CM_alpha'       : results.segments.cruise.conditions.static_stability.derivatives.CM_alpha[0], 
    #'CM_beta'        : results.segments.cruise.conditions.static_stability.derivatives.CM_beta[0], 
    #'CM_delta_a'     : results.segments.cruise.conditions.static_stability.derivatives.CM_delta_a[0], 
    #'CM_delta_e'     : results.segments.cruise.conditions.static_stability.derivatives.CM_delta_e[0], 
    #'CM_delta_r'     : results.segments.cruise.conditions.static_stability.derivatives.CM_delta_r[0], 
    #'CM_delta_f'     : results.segments.cruise.conditions.static_stability.derivatives.CM_delta_f[0], 
    #'CM_delta_s'     : results.segments.cruise.conditions.static_stability.derivatives.CM_delta_s[0], 
    #'CN_alpha'       : results.segments.cruise.conditions.static_stability.derivatives.CN_alpha[0], 
    #'CN_beta'        : results.segments.cruise.conditions.static_stability.derivatives.CN_beta[0], 
    #'CN_delta_a'     : results.segments.cruise.conditions.static_stability.derivatives.CN_delta_a[0], 
    #'CN_delta_e'     : results.segments.cruise.conditions.static_stability.derivatives.CN_delta_e[0], 
    #'CN_delta_r'     : results.segments.cruise.conditions.static_stability.derivatives.CN_delta_r[0], 
    #'CN_delta_f'     : results.segments.cruise.conditions.static_stability.derivatives.CN_delta_f[0], 
    #'CN_delta_s'     : results.segments.cruise.conditions.static_stability.derivatives.CN_delta_s[0], 
    #'Clift_u'        : results.segments.cruise.conditions.static_stability.derivatives.Clift_u[0], 
    #'Clift_v'        : results.segments.cruise.conditions.static_stability.derivatives.Clift_v[0], 
    #'Clift_w'        : results.segments.cruise.conditions.static_stability.derivatives.Clift_w[0], 
    #'Cdrag_u'        : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_u[0], 
    #'Cdrag_v'        : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_v[0], 
    #'Cdrag_w'        : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_w[0], 
    #'CX_u'           : results.segments.cruise.conditions.static_stability.derivatives.CX_u[0], 
    #'CX_v'           : results.segments.cruise.conditions.static_stability.derivatives.CX_v[0], 
    #'CX_w'           : results.segments.cruise.conditions.static_stability.derivatives.CX_w[0], 
    #'CY_u'           : results.segments.cruise.conditions.static_stability.derivatives.CY_u[0], 
    #'CY_v'           : results.segments.cruise.conditions.static_stability.derivatives.CY_v[0], 
    #'CY_w'           : results.segments.cruise.conditions.static_stability.derivatives.CY_w[0], 
    #'CZ_u'           : results.segments.cruise.conditions.static_stability.derivatives.CZ_u[0], 
    #'CZ_v'           : results.segments.cruise.conditions.static_stability.derivatives.CZ_v[0], 
    #'CZ_w'           : results.segments.cruise.conditions.static_stability.derivatives.CZ_w[0], 
    #'CL_u'           : results.segments.cruise.conditions.static_stability.derivatives.CL_u[0], 
    #'CL_v'           : results.segments.cruise.conditions.static_stability.derivatives.CL_v[0], 
    #'CL_w'           : results.segments.cruise.conditions.static_stability.derivatives.CL_w[0], 
    #'CM_u'           : results.segments.cruise.conditions.static_stability.derivatives.CM_u[0], 
    #'CM_v'           : results.segments.cruise.conditions.static_stability.derivatives.CM_v[0], 
    #'CM_w'           : results.segments.cruise.conditions.static_stability.derivatives.CM_w[0], 
    #'CN_u'           : results.segments.cruise.conditions.static_stability.derivatives.CN_u[0], 
    #'CN_v'           : results.segments.cruise.conditions.static_stability.derivatives.CN_v[0], 
    #'CN_w'           : results.segments.cruise.conditions.static_stability.derivatives.CN_w[0], 
    #'CZ_alpha_dot'   : results.segments.cruise.conditions.static_stability.derivatives.CZ_alpha_dot[0], 
    #'CM_alpha_dot'   : results.segments.cruise.conditions.static_stability.derivatives.CM_alpha_dot[0], 
    #'Clift_p'        : results.segments.cruise.conditions.static_stability.derivatives.Clift_p[0], 
    #'Clift_q'        : results.segments.cruise.conditions.static_stability.derivatives.Clift_q[0], 
    #'Clift_r'        : results.segments.cruise.conditions.static_stability.derivatives.Clift_r[0], 
    #'Cdrag_p'        : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_p[0], 
    #'Cdrag_q'        : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_q[0], 
    #'Cdrag_r'        : results.segments.cruise.conditions.static_stability.derivatives.Cdrag_r[0], 
    #'CX_p'           : results.segments.cruise.conditions.static_stability.derivatives.CX_p[0], 
    #'CX_q'           : results.segments.cruise.conditions.static_stability.derivatives.CX_q[0], 
    #'CX_r'           : results.segments.cruise.conditions.static_stability.derivatives.CX_r[0], 
    #'CY_p'           : results.segments.cruise.conditions.static_stability.derivatives.CY_p[0], 
    #'CY_q'           : results.segments.cruise.conditions.static_stability.derivatives.CY_q[0], 
    #'CY_r'           : results.segments.cruise.conditions.static_stability.derivatives.CY_r[0], 
    #'CZ_p'           : results.segments.cruise.conditions.static_stability.derivatives.CZ_p[0], 
    #'CZ_q'           : results.segments.cruise.conditions.static_stability.derivatives.CZ_q[0], 
    #'CZ_r'           : results.segments.cruise.conditions.static_stability.derivatives.CZ_r[0], 
    #'CL_p'           : results.segments.cruise.conditions.static_stability.derivatives.CL_p[0], 
    #'CL_q'           : results.segments.cruise.conditions.static_stability.derivatives.CL_q[0], 
    #'CL_r'           : results.segments.cruise.conditions.static_stability.derivatives.CL_r[0], 
    #'CM_p'           : results.segments.cruise.conditions.static_stability.derivatives.CM_p[0], 
    #'CM_q'           : results.segments.cruise.conditions.static_stability.derivatives.CM_q[0], 
    #'CM_r'           : results.segments.cruise.conditions.static_stability.derivatives.CM_r[0], 
    #'CN_p'           : results.segments.cruise.conditions.static_stability.derivatives.CN_p[0], 
    #'CN_q'           : results.segments.cruise.conditions.static_stability.derivatives.CN_q[0], 
    #'CN_r'           : results.segments.cruise.conditions.static_stability.derivatives.CN_r[0]})   
    
    #df.to_csv('Output_Navion_Stability_VLM.csv')    

    return  

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
    vehicle.flight_envelope.ultimate_load             = 5.7
    vehicle.flight_envelope.limit_load                = 3.8
    vehicle.reference_area                            = 17.112 
    vehicle.passengers                                = 2
    vehicle.design_dynamic_pressure                   = 1929.16080736607
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------   

    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
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

    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator + '..' + separator+ '..' + separator+ '..' + separator

    tip_airfoil                           = RCAIDE.Library.Components.Airfoils.NACA_4_Series_Airfoil()
    tip_airfoil.NACA_4_Series_code        = '6410'      
    tip_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'NACA_6410.txt' 
   
    root_airfoil                          = RCAIDE.Library.Components.Airfoils.NACA_4_Series_Airfoil()
    root_airfoil.NACA_4_Series_code       = '4415'   
    root_airfoil.coordinate_file          = rel_path + 'Airfoils' + separator + 'NACA_4415.txt' 
    
    # Wing Segments 
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'root_segment'
    segment.percent_span_location         = 0.0
    segment.twist                         = 2 * Units.degrees  
    segment.root_chord_percent            = 1.0
    segment.dihedral_outboard             = 7.5 * Units.degrees  
    segment.sweeps.quarter_chord          = 0.165 * Units.degrees  
    segment.thickness_to_chord            = .15 
    wing.append_segment(segment)  
         
    segment                               = RCAIDE.Library.Components.Wings.Segment()
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
    flap                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.5
    flap.deflection               = 0.0 * Units.degrees 
    flap.chord_fraction           = 0.20
    wing.append_control_surface(flap)  
    

    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
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
    wing                                  = RCAIDE.Library.Components.Wings.Wing()
    wing.tag                              = 'horizontal_stabilizer'  
    wing.sweeps.leading_edge              = 6 * Units.degrees 
    wing.thickness_to_chord               = 0.12
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
    
    elevator                              = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                          = 'elevator'
    elevator.span_fraction_start          = 0.1
    elevator.span_fraction_end            = 0.9
    elevator.deflection                   = 0.0  * Units.deg
    elevator.chord_fraction               = 0.3
    wing.append_control_surface(elevator)       

    RCAIDE.Library.Methods.Geometry.Planform.wing_planform(wing)     

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------ 
    wing                                  = RCAIDE.Library.Components.Wings.Wing()
    wing.tag                              = 'vertical_stabilizer'   
    wing.sweeps.leading_edge              = 20 * Units.degrees 
    wing.thickness_to_chord               = 0.125
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
    
    rudder                                = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
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
    fuselage = RCAIDE.Library.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage'
    fuselage.seats_abreast                      = 2
    fuselage.lengths.total                      = 8.349950916 
    fuselage.width                              = 1.22028016 
    fuselage.heights.maximum                    = 1.634415138  
    fuselage.areas.wetted                       = 12. # ESTIMATED 
    fuselage.areas.front_projected              = fuselage.width*fuselage.heights.maximum
    fuselage.effective_diameter                 = 1.22028016 

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0
    segment.percent_z_location                  = 0
    segment.height                              = 0.529255748
    segment.width                               = 0.575603849
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  =  0.028527593
    segment.percent_z_location                  =  0
    segment.height                              =  0.737072721
    segment.width                               =  0.921265952 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.187342754 
    segment.percent_z_location                  = 0 
    segment.height                              = 1.174231852 
    segment.width                               = 1.196956212
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_3'
    segment.percent_x_location                  = 0.242034847 
    segment.percent_z_location                  = 0.011503528 
    segment.height                              = 1.450221906 
    segment.width                               = 1.173932059 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_4'
    segment.percent_x_location                  = 0.296715183 
    segment.percent_z_location                  = 0.015984303 
    segment.height                              = 1.634415138 
    segment.width                               = 1.22028016 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_5'
    segment.percent_x_location                  = 0.510275342 
    segment.percent_z_location                  = -0.005
    segment.height                              = 1.082135236 
    segment.width                               = 1.013062774 
    fuselage.Segments.append(segment)

    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_6'
    segment.percent_x_location                  = 0.833284347 
    segment.percent_z_location                  = 0.014138855 
    segment.height                              = 0.621652157 
    segment.width                               = 0.414134978
    fuselage.Segments.append(segment)
 
    # Segment
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment()
    segment.tag                                 = 'segment_7'
    segment.percent_x_location                  = 1
    segment.percent_z_location                  = 0.018978667 
    segment.height                              = 0.092096616 
    segment.width                               = 0.046048308 
    fuselage.Segments.append(segment)
    
    # add to vehicle
    vehicle.append_component(fuselage) 

    # ########################################################  Energy Network  #########################################################  
    net                                         = RCAIDE.Framework.Networks.Fuel()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                   = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Fuel Tank & Fuel
    #------------------------------------------------------------------------------------------------------------------------------------       
    fuel_tank                                             = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank() 
    fuel_tank.origin                                      = vehicle.wings.main_wing.origin  
    fuel_tank.fuel                                        = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline() 
    fuel_tank.fuel.mass_properties.mass                   = 319 *Units.lbs 
    fuel_tank.fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel_tank.volume                                      = fuel_tank.fuel.mass_properties.mass/fuel_tank.fuel.density   
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
    design_propeller(prop)    
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
     
    return configs 

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config, configs)
        analyses[tag] = analysis

    return analyses


def base_analysis(vehicle, configs):

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
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle                             = vehicle
    aerodynamics.settings.number_spanwise_vortices   = 30
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    aerodynamics.settings.model_fuselage             = True                
    aerodynamics.settings.model_nacelle              = True
    aerodynamics.settings.use_surrogate              = False  
    analyses.append(aerodynamics) 
      
    stability                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method() 
    stability.settings.discretize_control_surfaces  = True
    stability.settings.model_fuselage               = True                
    stability.settings.model_nacelle                = True
    stability.settings.use_surrogate                = False  
        
    stability.configuration                         = configs
    stability.vehicle                               = vehicle
    analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Framework.Analyses.Energy.Energy()
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

    # done!
    return analyses  
 
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
    base_segment.state.numerics.number_of_control_points    = 3

    # ------------------------------------------------------------------
    #   Climb Segment : Constant Speed Constant Rate
    # ------------------------------------------------------------------ 
    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment) 
    segment.tag = "cruise"        
    segment.analyses.extend( analyses.base )      
    segment.altitude                                                            = 3000 * Units.feet
    segment.air_speed                                                           = 120 * Units['mph']
    segment.sideslip_angle                                                      = 0 * Units.degrees
                          
    # define flight dynamics to model                       
    segment.flight_dynamics.force_x                                             = True    
    segment.flight_dynamics.force_z                                             = True    
                
    # define flight controls               
    segment.assigned_control_variables.throttle.active                          = True           
    segment.assigned_control_variables.throttle.assigned_propulsors             = [['ice_propeller']]    
    segment.assigned_control_variables.body_angle.active                        = True   
    
    # Longidinal Flight Mechanics
    segment.flight_dynamics.moment_y                                            = True 
    segment.assigned_control_variables.elevator_deflection.active               = True    
    segment.assigned_control_variables.elevator_deflection.assigned_surfaces    = [['elevator']]
    segment.assigned_control_variables.elevator_deflection.initial_guess_values = [[0]]     
   
    # Lateral Flight Mechanics 
    segment.flight_dynamics.force_y                                             = True     
    segment.flight_dynamics.moment_x                                            = True
    segment.flight_dynamics.moment_z                                            = True 
    segment.assigned_control_variables.aileron_deflection.active                = True    
    segment.assigned_control_variables.aileron_deflection.assigned_surfaces     = [['aileron']]
    segment.assigned_control_variables.aileron_deflection.initial_guess_values  = [[0]] 
    segment.assigned_control_variables.rudder_deflection.active                 = True    
    segment.assigned_control_variables.rudder_deflection.assigned_surfaces      = [['rudder']]
    segment.assigned_control_variables.rudder_deflection.initial_guess_values   = [[0]]
    segment.assigned_control_variables.bank_angle.active                        = True    
    segment.assigned_control_variables.bank_angle.initial_guess_values          = [[0]]  
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