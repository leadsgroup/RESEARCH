# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units   
from RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor   import design_turbofan
from RCAIDE.Library.Methods.Stability.Center_of_Gravity            import compute_component_centers_of_gravity
from RCAIDE.Library.Methods.Geometry.Planform                     import segment_properties
from RCAIDE.Library.Plots                 import *     

# python imports 
import numpy as np  
from copy import deepcopy
import matplotlib.pyplot as plt  
import os   

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    # Step 1 design a vehicle
    vehicle  = vehicle_setup()
    
    # plot vehicle 
    plot_3d_vehicle(vehicle,
                    min_x_axis_limit            = -1,
                    max_x_axis_limit            = 100,
                    min_y_axis_limit            = -50,
                    max_y_axis_limit            = 50,
                    min_z_axis_limit            = -50,
                    max_z_axis_limit            = 50)      
    
    # Step 2 create aircraft configuration based on vehicle 
    configs  = configs_setup(vehicle)
    
    # Step 3 set up analysis
    analyses = analyses_setup(configs)
    
    # Step 4 set up a flight mission
    mission = mission_setup(analyses)
    missions = missions_setup(mission) 
    
    # Step 5 execute flight profile
    results = missions.base_mission.evaluate()  
    
    # Step 6 plot results 
    plot_mission(results)
    
     
    
    return


def vehicle_setup():
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Boeing_777-300'

    
    # ################################################# Vehicle-level Properties #################################################   
    vehicle.mass_properties.max_takeoff               = 299370 * Units.kilogram  
    vehicle.mass_properties.takeoff                   = 299370 * Units.kilogram    
    vehicle.mass_properties.operating_empty           = 160530 * Units.kilogram  
    vehicle.mass_properties.max_zero_fuel             = 224530 * Units.kilogram 
    vehicle.mass_properties.cargo                     = 64000 * Units.kilogram  
    vehicle.envelope.ultimate_load                    = 3.75 
    vehicle.envelope.limit_load                       = 2.5  
    vehicle.reference_area                            = 436.80 * Units['meters**2']   
    vehicle.passengers                                = 368 
    vehicle.systems.control                           = "fully powered" 
    vehicle.systems.accessories                       = "medium range"


    # ################################################# Landing Gear #############################################################   
    # ------------------------------------------------------------------        
    #  Landing Gear
    # ------------------------------------------------------------------  
    landing_gear                    = RCAIDE.Library.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag                = "landing_gear" 
    landing_gear.main_tire_diameter = 1.12000 * Units.m  
    landing_gear.nose_tire_diameter = 0.6858 * Units.m  
    landing_gear.main_strut_length  = 1.8 * Units.m 
    landing_gear.nose_strut_length  = 1.3 * Units.m 
    landing_gear.main_units         = 2 
    landing_gear.nose_units         = 1 
    landing_gear.main_wheels        = 6 
    landing_gear.nose_wheels        = 2     
    vehicle.landing_gear            = landing_gear    
    
    # ################################################# Wings ##################################################################### 
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------

    wing                                  = RCAIDE.Library.Components.Wings.Main_Wing()
    wing.tag                              = 'main_wing'
    wing.sweeps.leading_edge              = 35.214 * Units.deg
    wing.thickness_to_chord               = 0.10  
    wing.spans.projected                  = 64.8 
    wing.chords.root                      = 15.1335 * Units.meter
    wing.chords.tip                       = 1.8055 * Units.meter
    wing.taper                            = wing.chords.tip / wing.chords.root
    wing.chords.mean_aerodynamic          = 15.1335 / 2 # INCORRECT 
    wing.areas.reference                  = 436.80 
    wing.areas.exposed                    = 2*wing.areas.reference *0.8
    wing.areas.wetted                     = 2*wing.areas.reference *0.8 
    wing.total_length                     = wing.chords.root 
    wing.aspect_ratio                     = (wing.spans.projected ** 2) /  wing.areas.reference 
    wing.twists.root                      = 4.0 * Units.degrees 
    wing.twists.tip                       = 0.0 * Units.degrees 
    wing.exposed_root_chord_offset        = 6.2 / 2
    wing.origin                           = [[25 ,0, -0.75]]
    wing.aerodynamic_center               = [25+ 0.25*wing.chords.root ,0, -0.75]  
    wing.vertical                         = False
    wing.dihedral                         = 7.5 * Units.degrees 
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
    segment.twist                         = 3 * Units.degrees 
    segment.root_chord_percent            = 1.0
    segment.thickness_to_chord            = 0.10
    segment.dihedral_outboard             = 7.5 * Units.degrees 
    segment.sweeps.leading_edge           = 35.214 * Units.deg 
    segment.append_airfoil(root_airfoil)
    wing.append_segment(segment)

    inboard_airfoil                       = RCAIDE.Library.Components.Airfoils.Airfoil()
    inboard_airfoil.coordinate_file       = rel_path+ 'Airfoils' + separator + 'B737b.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'inboard'
    segment.percent_span_location         = 0.26382
    segment.twist                         = 2 * Units.degrees 
    segment.root_chord_percent            = 0.635424
    segment.thickness_to_chord            = 0.10 
    segment.dihedral_outboard             = 7.5 * Units.degrees 
    segment.sweeps.leading_edge           = 35.214 * Units.deg 
    segment.append_airfoil(inboard_airfoil)
    wing.append_segment(segment)
    

    outboard_airfoil                      = RCAIDE.Library.Components.Airfoils.Airfoil()
    outboard_airfoil.coordinate_file      = rel_path+ 'Airfoils' + separator + 'B737c.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'outboard'
    segment.percent_span_location         = 0.3772
    segment.twist                         = 2 * Units.degrees 
    segment.root_chord_percent            = 0.5154
    segment.thickness_to_chord            = 0.10 
    segment.dihedral_outboard             = 7.5 * Units.degrees 
    segment.sweeps.leading_edge           = 35.214 * Units.deg 
    segment.append_airfoil(outboard_airfoil)
    wing.append_segment(segment) 
 

    tip_airfoil                           =  RCAIDE.Library.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = rel_path + 'Airfoils' + separator + 'B737d.txt'
    segment                               = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                           = 'Tip'
    segment.percent_span_location         = 1.0 
    segment.twist                         = 0 * Units.degrees  
    segment.root_chord_percent            = 0.11930 
    segment.thickness_to_chord            = 0.10   
    segment.dihedral_outboard             = 0 * Units.degrees      
    segment.sweeps.quarter_chord          = 0 * Units.degrees     
    segment.append_airfoil(tip_airfoil) 
    wing.append_segment(segment)
    

    # control surfaces -------------------------------------------
    slat                                   = RCAIDE.Library.Components.Wings.Control_Surfaces.Slat()
    slat.tag                               = 'slat'
    slat.span_fraction_start               = 0.2
    slat.span_fraction_end                 = 0.963
    slat.deflection                        = 0.0 * Units.degrees
    slat.chord_fraction                    = 0.075
    wing.append_control_surface(slat)

    #inboard_flap                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    #inboard_flap.tag                      = 'inboard_flap'
    #inboard_flap.span_fraction_start      = 0.05
    #inboard_flap.span_fraction_end        = 0.2
    #inboard_flap.deflection               = 0.0 * Units.degrees
    #inboard_flap.configuration_type       = 'tripple_slotted'
    #inboard_flap.chord_fraction           = 0.30
    #wing.append_control_surface(inboard_flap)
    

    flap                                  = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    flap.tag                              = 'flap'
    flap.span_fraction_start              = 0.4
    flap.span_fraction_end                = 0.8
    flap.deflection                       = 0.0 * Units.degrees
    flap.configuration_type               = 'double_slotted'
    flap.chord_fraction                   = 0.30
    wing.append_control_surface(flap) 

    aileron                               = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                           = 'aileron'
    aileron.span_fraction_start           = 0.8
    aileron.span_fraction_end             = 0.963
    aileron.deflection                    = 0.0 * Units.degrees
    aileron.chord_fraction                = 0.16
    wing.append_control_surface(aileron)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing,update_wet_areas=True, update_ref_areas=True)   
    
    # add to vehicle
    vehicle.append_component(wing)
    

    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

    wing                         = RCAIDE.Library.Components.Wings.Horizontal_Tail()
    wing.tag                     = 'horizontal_stabilizer' 
    wing.sweeps.leading_edge     = 38.6598082 * Units.degrees 
    wing.thickness_to_chord      = 0.12 
    wing.spans.projected         = 10.79076233 *2 
    wing.chords.root             = 6.9423519 
    wing.chords.tip              = 2.36018797  
    wing.taper                   = wing.chords.tip  / wing.chords.root
    wing.chords.mean_aerodynamic = 6.942 /2 # Incorrect  
    wing.areas.reference         = ( wing.chords.root +  wing.chords.tip) *wing.spans.projected / 2   
    wing.areas.exposed           = 2*wing.areas.reference *0.9
    wing.areas.wetted            = 2*wing.areas.reference *0.9 
    wing.total_length            = wing.chords.root 
    wing.aspect_ratio            = (wing.spans.projected ** 2) /  wing.areas.reference 
    wing.twists.root             = 1.0 * Units.degrees
    wing.twists.tip              = 1.0 * Units.degrees  
    wing.exposed_root_chord_offset   = 6.2 / 2    
    wing.dihedral                    = 7.6 * Units.degrees 
    wing.origin                      = [[63.03027968  , 0  ,  2.1745719]] 
    wing.aerodynamic_center          = [63.03027968 + 0.25*wing.chords.root, 0 , 2.1745719]
    wing.vertical                    = False
    wing.symmetric                   = True 
    wing.dynamic_pressure_ratio      = 1.0 
 
    # Wing Segments
    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'root_segment'
    segment.percent_span_location  = 0.0
    segment.twist                  = 1.0 * Units.degrees
    segment.root_chord_percent     = 1.0
    segment.dihedral_outboard      = 7.6 * Units.degrees 
    segment.sweeps.leading_edge    = 38.6598082  * Units.degrees 
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)

    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'tip_segment'
    segment.percent_span_location  = 1.
    segment.twist                  = 1.0 * Units.degrees
    segment.root_chord_percent     = wing.taper            
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 0 * Units.degrees  
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)
    
    # control surfaces -------------------------------------------
    elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.1 
    elevator.span_fraction_end     = 0.9  
    elevator.deflection            = 0 
    elevator.chord_fraction        = 0.33  
    wing.append_control_surface(elevator)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing,update_wet_areas=True, update_ref_areas=True)   

    # add to vehicle
    vehicle.append_component(wing)
    

    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------ 
    wing                         =  RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag                     =  'vertical_stabilizer' 
    wing.sweeps.leading_edge     = 48.37  * Units.degrees 
    wing.thickness_to_chord      = 0.12  
    wing.spans.projected         = 10.92735426 
    wing.chords.root             = 9.43826157 
    wing.chords.tip              = 2.36081054 
    wing.total_length            = 9.43826157  
    wing.areas.reference         = ( wing.chords.root +  wing.chords.tip) *wing.spans.projected / 2   
    wing.aspect_ratio            = (wing.spans.projected ** 2) /  wing.areas.reference 
    wing.taper                   = wing.chords.tip  / wing.chords.root
    wing.chords.mean_aerodynamic = 9.43826157 /2  
    wing.areas.exposed           = 2*wing.areas.reference *0.9
    wing.areas.wetted            = 2*wing.areas.reference *0.9  
    wing.exposed_root_chord_offset   = 6.2 / 2  
    wing.aspect_ratio            = (wing.spans.projected ** 2) /  wing.areas.reference 
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees 
    wing.origin                  = [[59.90743535, 0 , 2.45865471 ]]  
    wing.aerodynamic_center      = [59.90743535 , 0 , 2.45865471 +  0.25*wing.chords.root]  
    wing.vertical                = True
    wing.symmetric               = False
    wing.t_tail                  = False 
    wing.dynamic_pressure_ratio  = 1.0
     

    # Wing Segments
    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'root_segment'
    segment.percent_span_location  = 0.0
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 1.0
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.leading_edge    = 48.37  * Units.degrees 
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)

    segment                        = RCAIDE.Library.Components.Wings.Segment()
    segment.tag                    = 'tip_segment'
    segment.percent_span_location  = 1.
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = wing.taper            
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 0 * Units.degrees  
    segment.thickness_to_chord     = .1
    wing.append_segment(segment) 
        

    # control surfaces -------------------------------------------
    rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
    rudder.tag                   = 'rudder'
    rudder.span_fraction_start   = 0.1 
    rudder.span_fraction_end     = 0.9  
    rudder.deflection            = 0 
    rudder.chord_fraction        = 0.33  
    wing.append_control_surface(rudder)

    # Fill out more segment properties automatically
    wing = segment_properties(wing,update_wet_areas=True, update_ref_areas=True)   

    # add to vehicle
    vehicle.append_component(wing)    
    
    
    # ################################################# Fuselage ################################################################ 
        
    fuselage                                    = RCAIDE.Library.Components.Fuselages.Tube_Fuselage() 
    fuselage.number_coach_seats                 = vehicle.passengers

    fuselage.seats_abreast                      = 6
    fuselage.seat_pitch                         = 1     * Units.meter 
    fuselage.fineness.nose                      = 1.6
    fuselage.fineness.tail                      = 2. 
    fuselage.lengths.nose                       = 8.9  * Units.meter
    fuselage.lengths.tail                       = 19   * Units.meter
    fuselage.lengths.total                      = 73.08 * Units.meter  
    fuselage.lengths.fore_space                 = 6.    * Units.meter
    fuselage.lengths.aft_space                  = 5.    * Units.meter 
    fuselage.width                              = 6.2   * Units.meter
    fuselage.heights.maximum                    = 6.2  * Units.meter
    fuselage.effective_diameter                 = 6.2   * Units.meter
    fuselage.areas.side_projected               = fuselage.heights.maximum * fuselage.lengths.total * Units['meters**2'] 
    fuselage.areas.wetted                       = np.pi * fuselage.width/2 * fuselage.lengths.total * Units['meters**2'] 
    fuselage.areas.front_projected              = np.pi * fuselage.width/2      * Units['meters**2']  
    fuselage.differential_pressure              = 5.0e4 * Units.pascal
    fuselage.heights.at_quarter_length          = fuselage.heights.maximum * Units.meter
    fuselage.heights.at_three_quarters_length   = fuselage.heights.maximum * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = fuselage.heights.maximum* Units.meter
    
    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_1'    
    segment.percent_x_location                  = 0.0 
    segment.percent_z_location                  = -4.5718274092838756e-05 
    segment.height                              = 0.0 
    segment.width                               = 0.0 
    fuselage.Segments.append(segment)   
     

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_2'    
    segment.percent_x_location                  = 0.0006849315068493151 
    segment.percent_z_location                  = 0.00023201116822421539 
    segment.height                              = 0.506612607796538 
    segment.width                               = 0.29486232497524273 
    fuselage.Segments.append(segment)
    

    ## Segment  
    #segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    #segment.tag                                 = 'segment_3'    
    #segment.percent_x_location                  = 0.0013698630136986301 
    #segment.percent_z_location                  = 0.0005097406105412695 
    #segment.height                              = 1.0105126311743573 
    #segment.width                               = 0.5860983155721387 
    #fuselage.Segments.append(segment)
    

    ## Segment  
    #segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    #segment.tag                                 = 'segment_4'    
    #segment.percent_x_location                  = 0.0027397260273972603 
    #segment.percent_z_location                  = 0.00106519949517538 
    #segment.height                              = 1.1902153719637647 
    #segment.width                               = 0.9369298914718769 
    #fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_5'    
    segment.percent_x_location                  = 0.00410958904109589 
    segment.percent_z_location                  = 0.0015206583798094884 
    segment.height                              = 1.3699181127531718 
    segment.width                               = 1.1328207660995513 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_6'    
    segment.percent_x_location                  = 0.00684931506849315 
    segment.percent_z_location                  = 0.0016650708841355363 
    segment.height                              = 1.7736138256504295 
    segment.width                               = 1.5246025153549005 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_7'    
    segment.percent_x_location                  = 0.00958904109589041 
    segment.percent_z_location                  = 0.0014382635613239057 
    segment.height                              = 1.9377114437855858 
    segment.width                               = 1.8057694327532068 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_8'    
    segment.percent_x_location                  = 0.0136986301369863 
    segment.percent_z_location                  = 0.002076051294030489 
    segment.height                              = 2.2895051006880056 
    segment.width                               = 2.1590037827079334 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_9'    
    segment.percent_x_location                  = 0.0410958904109589 
    segment.percent_z_location                  = 0.0054206755627905015 
    segment.height                              = 4.137475539696899 
    segment.width                               = 3.7836714348574154 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_10'    
    segment.percent_x_location                  = 0.0684931506849315 
    segment.percent_z_location                  = 0.007361830512116725 
    segment.height                              = 5.180529836017472 
    segment.width                               = 5.038003416722632 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_11'    
    segment.percent_x_location                  = 0.0958904109589041 
    segment.percent_z_location                  = 0.008342605799103458 
    segment.height                              = 5.779624874430869 
    segment.width                               = 5.926942646043528 
    fuselage.Segments.append(segment)

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_12'    
    segment.percent_x_location                  = 0.1095890410958904 
    segment.percent_z_location                  = 0.008852308594930632 
    segment.height                              = 5.936562615398982 
    segment.width                               = 6.198721376573906 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_13'    
    segment.percent_x_location                  = 0.1232876712328767 
    segment.percent_z_location                  = 0.009418759029858625 
    segment.height                              = 6.019264378898468 
    segment.width                               = 6.283228699551572 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_14'    
    segment.percent_x_location                  = 0.136986301369863 
    segment.percent_z_location                  = 0.009773342446019904 
    segment.height                              = 6.071033557658016 
    segment.width                               = 6.283228699551572 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_15'    
    segment.percent_x_location                  = 0.15753424657534246 
    segment.percent_z_location                  = 0.010291172676447542 
    segment.height                              = 6.14663677130045 
    segment.width                               = 6.283228699551572 
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_16'    
    segment.percent_x_location                  = 0.7534246575342466
    segment.percent_z_location                  = 0.010385204131418418
    segment.height                              = 6.132908178874702
    segment.width                               = 6.283228699551572
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_17'    
    segment.percent_x_location                  =  0.7808219178082192
    segment.percent_z_location                  =  0.011034594257023024
    segment.height                              =  6.009141025548721
    segment.width                               =  6.179361462712953
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_18'    
    segment.percent_x_location                  =  0.7945205479452054
    segment.percent_z_location                  =  0.011784853626199271
    segment.height                              =  5.796015915164266
    segment.width                               =  6.063584270320871
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_19'    
    segment.percent_x_location                  = 0.8493150684931506
    segment.percent_z_location                  = 0.017830396113581042
    segment.height                              = 4.328719183639479
    segment.width                               = 5.218170689274808
    fuselage.Segments.append(segment)
    

    # Segment  
    segment                                     = RCAIDE.Library.Components.Fuselages.Segment() 
    segment.tag                                 = 'segment_20'    
    segment.percent_x_location                  = 1.0
    segment.percent_z_location                  = 0.03287176935416379
    segment.height                              = 0.6643988043360767
    segment.width                               = 0.3578777817066474
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
     
    #------------------------------------------------------------------------------------------------------------------------- 
    #  Energy Source: Fuel Tank
    #------------------------------------------------------------------------------------------------------------------------- 
    # fuel tank
    fuel_tank                                   = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.tag                               = 'b777_fuel_tank'
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


    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------         
    turbofan                                    = RCAIDE.Library.Components.Propulsors.Turbofan() 
    turbofan.tag                                = 'starboard_ge90_propulsor'
    turbofan.active_fuel_tanks                  = ['b777_fuel_tank']   
    turbofan.origin                             = [[ 25.72797886 , 9.69802 , -2.04  ]]
    turbofan.mass_properties.mass               = 7893
    turbofan.engine_length                      = 7.29
    turbofan.bypass_ratio                       = 9
    turbofan.design_altitude                    = 35000.0*Units.ft
    turbofan.design_mach_number                 = 0.78   
    turbofan.design_thrust                      = 80000  * Units.N  


    # fan                
    fan                                         = RCAIDE.Library.Components.Propulsors.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.7   
    turbofan.fan                                = fan        

    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Library.Attributes.Gases.Air() 

    
    # Ram inlet 
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

    ## high pressure compressor  
    #medium_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    #medium_pressure_compressor.tag                   = 'hpc'
    #medium_pressure_compressor.polytropic_efficiency = 0.91
    #medium_pressure_compressor.pressure_ratio        = 12.38 
    #turbofan.high_pressure_compressor              = high_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 12.38 
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
    combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A()  
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
    # append propulsor to distribution line
    

    # Nacelle 
    nacelle                                     = RCAIDE.Library.Components.Nacelles.Body_of_Revolution_Nacelle()
    nacelle.diameter                            = 3.87
    nacelle.length                              = 5.27571429 
    nacelle.tag                                 = 'nacelle_1'
    nacelle.inlet_diameter                      = 3.1
    nacelle.origin                              = [[26.72797886 , 9.69802 , -2.04 ]] 
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
    turbofan_2.active_fuel_tanks                = ['b777_fuel_tank'] 
    turbofan_2.tag                              = 'port_ge90_propulsor' 
    turbofan_2.origin                           = [[ 25.72797886 , -9.69802 , -2.04  ]]   # change origin 
    turbofan_2.nacelle.origin                   = [[26.72797886 , -9.69802 , -2.04 ]]  
         
    # append propulsor to distribution line 
    fuel_line.propulsors.append(turbofan_2)
  
    # Append fuel line to Network      
    net.fuel_lines.append(fuel_line)   

    # Append energy network to aircraft 
    vehicle.append_energy_network(net)

    
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
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base' 
    base_config.landing_gear.gear_condition                      = 'up'
    configs.append(base_config)

    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'cruise'
    configs.append(config)


    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'takeoff'
    #config.wings['main_wing'].control_surfaces.flap.deflection  = 20. * Units.deg
    #config.wings['main_wing'].control_surfaces.slat.deflection  = 25. * Units.deg 
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['starboard_ge90_propulsor'].fan.angular_velocity =  3470. * Units.rpm
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['port_ge90_propulsor'].fan.angular_velocity      =  3470. * Units.rpm
    config.landing_gear.gear_condition                          = 'down'       
    config.V2_VS_ratio = 1.21
    configs.append(config)

    
    # ------------------------------------------------------------------
    #   Cutback Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'cutback'
    #config.wings['main_wing'].control_surfaces.flap.deflection  = 20. * Units.deg
    #config.wings['main_wing'].control_surfaces.slat.deflection  = 20. * Units.deg
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['starboard_ge90_propulsor'].fan.angular_velocity =  2780. * Units.rpm
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['port_ge90_propulsor'].fan.angular_velocity      =  2780. * Units.rpm
    config.landing_gear.gear_condition                          = 'up'       
    configs.append(config)   
    
        
    
    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'landing'
    #config.wings['main_wing'].control_surfaces.flap.deflection  = 30. * Units.deg
    #config.wings['main_wing'].control_surfaces.slat.deflection  = 25. * Units.deg
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['starboard_ge90_propulsor'].fan.angular_velocity =  2030. * Units.rpm
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['port_ge90_propulsor'].fan.angular_velocity      =  2030. * Units.rpm
    config.landing_gear.gear_condition                          = 'down'   
    config.Vref_VS_ratio = 1.23
    configs.append(config)   
     
    # ------------------------------------------------------------------
    #   Short Field Takeoff Configuration
    # ------------------------------------------------------------------ 

    config = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag = 'short_field_takeoff'    
    #config.wings['main_wing'].control_surfaces.flap.deflection  = 20. * Units.deg
    #config.wings['main_wing'].control_surfaces.slat.deflection  = 25. * Units.deg
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['starboard_ge90_propulsor'].fan.angular_velocity =  3470. * Units.rpm
    #config.networks.turbofan_engine.fuel_lines['fuel_line'].propulsors['port_ge90_propulsor'].fan.angular_velocity      =  3470. * Units.rpm
    config.landing_gear.gear_condition                          = 'down'   
    config.V2_VS_ratio = 1.21 
    configs.append(config)        
    
    
    return configs

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

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
    aerodynamics = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method()
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.number_of_spanwise_vortices   = 30
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

    segment = Segments.Ground.Takeoff(base_segment)
    segment.tag = "Takeoff" 
    segment.analyses.extend( analyses.takeoff )
    segment.velocity_start           = 10.* Units.knots
    segment.velocity_end             = 125.0 * Units['m/s']
    segment.friction_coefficient     = 0.04
    segment.altitude                 = 0.0   
    mission.append_segment(segment)

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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)


    ## ------------------------------------------------------------------    
    ##   Cruise Segment: Constant Speed Constant Altitude
    ## ------------------------------------------------------------------    

    #segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    #segment.tag = "cruise" 
    #segment.analyses.extend( analyses.cruise ) 
    #segment.altitude                                      = 10.668 * Units.km  
    #segment.air_speed                                     = 230.412 * Units['m/s']
    #segment.distance                                      = 2000 * Units.nmi   
    
    ## define flight dynamics to model 
    #segment.flight_dynamics.force_x                       = True  
    #segment.flight_dynamics.force_z                       = True     
    
    ## define flight controls 
    #segment.assigned_control_variables.throttle.active               = True           
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
    #segment.assigned_control_variables.body_angle.active             = True                
    
    #mission.append_segment(segment)


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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
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
    #segment.assigned_control_variables.throttle.assigned_propulsors  = [['starboard_ge90_propulsor','port_ge90_propulsor']] 
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