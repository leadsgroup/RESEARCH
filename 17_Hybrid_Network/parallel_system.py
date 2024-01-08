'''

The script below documents how to set up and plot the results of a flight analysis of a transonic 
passenger carrying aircraft. Here, the Boeing 737-800 model is used. 

''' 
 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------e

# RCAIDE imports
import RCAIDE
from RCAIDE.Core    import Units
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform   import segment_properties
#from RCAIDE.Methods.Propulsion.Design import design_turbofan , size_optimal_motor

# python imports
import numpy as np
import os

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    motor_work = 5E3
    shaft_power_offtake = 4E4
    
    # vehicle data
    vehicle = vehicle_setup()


def vehicle_setup():
    """This is the full physical definition of the vehicle, and is designed to be independent of the
    analyses that are selected"""
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------ 
    
    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Boeing_737-800'
    
    # ################################################# Vehicle-level Properties ################################################# 
    vehicle.mass_properties.max.max_takeoff     = 79015.8 * Units.kilogram
    vehicle.mass_properties.takeoff             = 79015.8 * Units.kilogram
    vehicle.mass_properties.operating_empty     = 62746.4 * Units.kilogram
    vehicle.mass_properties.max_zero_fuel       = 62732.0 * Units.kilogram
    vehicle.mass_properties.cargo               = 10000.  * Units.kilogram
    vehicle.envelope.ultimate_load              = 3.75 #Load Factor Limits: maximum G-force an aircraft can withstand
    vehicle.envelope.limit_load                 = 2.5 #Load factor with no damage on aircraft
    vehicle.reference_area                      = 124.862 * Units['meters**2']
    vehicle.passengers                          = 170
    vehicle.systems.control                     = "fully powered"
    vehicle.systems.accessories                 = "medium range"
    
    # ################################################# Landing Gear #############################################################   
    # ------------------------------------------------------------------        
    #  Landing Gear
    # ------------------------------------------------------------------
    landing_gear                    = RCAIDE.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag                = "main_landing_gear"
    landing_gear.main_tire_diameter = 1.12000 * Units.m
    landing_gear.nose_tire_diamter  = 0.6858 * Units.m
    landing_gear.main_strut_length  = 1.8 * Units.m
    landing_gear.nose_strut_length  = 1.3 * Units.m
    landing_gear.main_units         = 2 #Number of main landing gear
    landing_gear.nose_units         = 1
    landing_gear.main_wheels        = 2 #Number of wheels on the main landing gear
    landing_gear.nose_wheels        = 2
    vehicle.landing_gear            = landing_gear
    
    # ################################################# Wings ##################################################################### 
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------
    
    wing        = RCAIDE.Components.Wings.Main_Wing()
    wing.tag    = 'main_wing'
    wing.aspect_ratio   = 10.18
    wing.sweeps.quarter_chord   = 25 * Units.deg
    wing.thickness_to_chord     = 0.1
    wing.taper                  = 0.1
    wing.spans.projected        = 34.32
    wing.chords.root            = 7.760 * Units.meter
    wing.chord.tip              = 0.782 * Units.meter
    wing.chords.mean_aerodynamic= 4.235 * Units.meter
    wing.areas.reference        = 124.862
    wing.areas.wetted           = 225.08
    wing.twists.root            = 4.0 * Units.degrees
    wing.twists.tip             = 0.0 * Units.degrees
    wing.origin                 = [[13.61,0,-0.5]]
    wing.aerodynamic_center     = [0,0,0]
    wing.vertical               = False
    wing.symmetric              = True
    wing.high_lift              = True
    wing.dynamic_pressure_ratio = 1.0
    
    # Wing Segments
    root_airfoil            = RCAIDE.Components.Airfoils.Airfoil()
    ospath                  = os.path.abspath(__file__)  #Find absolute path of the __file__
    separator               = os.path.sep
    rel_path                = os.path.dirname(ospath) + separator
    root_airfoil.coordinate_file    = rel_path + '..' +separator + 'Airfoils' + separator + 'B737a.txt'
    
    segment                 = RCAIDE.Components.Wings.Segment()
    segment.tag             = 'Root'
    segment.percent_span_location   = 0.0
    segment.twist                   = 4. * Units.deg
    segment.root_chord_percent      = 1.
    segment.thickness_to_chord      = 0.1
    segment.dihedral_outboard       = 2.5 * Units.degrees
    segment.sweeps.quarter_chord    = 28.225 * Units.degrees
    segment.thickness_to_chord      = .1
    segment.append_airfoil(root_airfoil)
    wing.append_segment(segment)
    
    yehudi_airfoil                  = RCAIDE.Components.Airfoils.Airfoil()
    yehudi_airfoil.coordinate_file  = rel_path + '..' +separator + 'Airfoils' + separator + 'B737b.txt'
    segment                         = RCAIDE.Components.Wings.Segment()
    segment.tag                     = 'Yehudi'
    segment.percent_span_location   = 0.324
    segment.twist                   = 0.047193 * Units.deg
    segment.root_chord_percent      = 0.5
    segment.thickness_to_chord      = 0.1
    segment.dihedral_outboard       = 5.5
    segment.sweeps.quarter_chord    = 25. * Units.degrees
    segment.thickness_to_chord      = .1
    segment.appened_airfoil(yehudi_airfoil)
    wing.append_segment(segment)
    
    mid_airfoil                     = RCAIDE.Components.Airfoils.Airfoil()
    mid_airfoil.coordinate_file     = rel_path + '..' +separator + 'Airfoils' + separator + 'B737c.txt'
    segment                         = RCAIDE.Components.Wings.Segment()
    segment.tag                     = 'Section_2'
    segment.percent_span_location   = 0.963
    segment.twist                   = 0.00258 * Units.deg
    segment.root_chord_percent      = 0.220
    segment.thickness_to_chord      = 0.1
    segment.dihedral_outboard       = 5.5 * Units.degrees
    segment.sweeps.quarter_chord    = 56.75 * Units.degrees
    segment.thickness_to_chord      = .1
    segment.append_airfoil(mid_airfoil)
    wing.append_segment(segment)
    
    tip_airfoil                     = RCAIDE.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file     = rel_path + '..' +separator + 'Airfoils' + separator + 'B737d.txt'
    segment                         = RCAIDE.Components.Wings.Segment()
    segment.tag                     = 'Tip'
    segment.percent_span_location   = 1.
    segment.twist                   = 0. * Units.degrees
    segment.root_chord_percent      = 0.10077
    segment.thickness_to_chord      = 0.1
    segment.dihedral_outboard       = 0.
    segment.sweeps.quarter_chord    = 0.
    segment.thickness_to_chord      = .1
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)
    
    # control surfaces -------------------------------------------
    slat                    = RCAIDE.Components.Wings.Control_Surfaces.Slat()
    slat.tag                = 'slat'
    slat.span_fraction_start= 0.2
    slat.span_fraction_end  = 0.963
    slat.deflection         = 0.0 * Units.degrees
    slat.chrod_fraction     = 0.075
    wing.append_control_surface(slat)
    
    flap                    = RCAIDE.Components.Wings.Control_Surfaces.Flap()
    flap.tag                = 'flap'
    flap.span_fraction_start= 0.2
    flap.span_fraction_end  = 0.7
    flap.deflection         = 0.0 * Units.degrees
    flap.chord_fraction     = 0.30
    wing.append_control_surface(flap)
    
    aileron                 = RCAIDE.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag             ='aileron'
    aileron.span_fraction_start = 0.7
    aileron.span_fraction_end   = 0.963
    aileron.deflection          = 0.0 * Units.degrees
    aileron.chord_fraction      = 0.16
    wing.append_control_surface(aileron)
    
    
    # add to vehicle
    vehicle.append_component(wing)
    
    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------
    
    wing    = RCAIDE.Components.Wings.Horizontal_Tail()
    wing.tag= 'horizontal_stabilizer'
    
    wing.aspect_ratio           = 4.99
    wing.sweeps.quarter_chord   = 28.2250 * Units.deg
    wing.thickness_to_chord     = 0.08
    wing.taper                  = 0.3333
    wing.spans.projected        = 14.4
    wing.chords.root            = 4.2731
    wing.chords.tip             = 1.4243
    wing.chords.mean_aerodynamic= 8.0
    wing.area.reference         = 41.49
    wing.areas.exposed          = 59.354    
    wing.areas.wetted           = 71.81
    wing.twist.root             = 3.0 * Units.degrees
    wing.twists.tip             = 3.0 * Units.degrees
    wing.origin                 = [[33.02,0,1.466]]
    wing.aerodynamic_center     = [0,0,0]
    wing.vertical               = False
    wing.symmetric              = True
    wing.dynamic_pressure_ratio = 0.9
    
    
    # Wing Segments
    segment                         = RCAIDE.Components.Wings.Segment()
    segment.tag                     = 'root_segment'
    segment.percent_span_location   = 0.0
    segment.twist                   = 0. * Units.deg
    segment.root_chord_percent      = 1.0
    segment.diherdral_outboard      = 8.63 * Units.degrees
    segment.sweeps.quarter_chord    = 28.2250 * Units.degrees
    segment.thickness_to_chord      = .1
    wing.append_segment(segment)
    
    segment                        = RCAIDE.Components.Wings.Segment()
    segment.tag                    = 'tip_segment'
    segment.percent_span_location  = 1.
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 0.3333               
    segment.dihedral_outboard      = 0 * Units.degrees
    segment.sweeps.quarter_chord   = 0 * Units.degrees  
    segment.thickness_to_chord     = .1
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)
    
    # control surfaces ------------------------------------------- 
    elevator                    = RCAIDE.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                = 'elevator'
    elevator.span_fraction_start= 0.09
    elevator.span_fraction_end  = 0.92
    elevator.deflection         = 0.0 * Units.deg
    elevator.chord_fraction     = 0.3
    wing.append_control_surface(elevator)
    
    # add to vehicle
    vehicle.append_component(wing)
    
    
    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------
    
    wing = RCAIDE.Components.Wings.Vertical_Tail()
    wing,tag = 'vertical_stabilizer'
    
    wing.aspect_ratio           = 1.98865
    wing.sweeps.quarter_chord   = 31.2 * Units.deg
    wing.thickness_to_chord     = 0.08
    wing.taper                  = 0.1183
    
    wing.spans.projected        = 8.33
    wing.total_length           = wing.spans.projected
    
    wing.chords.root            = 10.1
    wing.chords.tip             = 1.20
    wing.chords.mean_aerodnyamic= 4.0
    
    wing.areas.reference        = 34.89
    wing.areas.wetted           = 57.25
    
    wing.twists.root            = 0.0 * Units.degrees
    wing.twist.tip              = 0.0 * Units.degrees
    
    wing.origin                 = [[26.944,0,1.54]]
    wing.aerodynamic_center     = [0,0,0]
    
    wing.vertical               = True
    wing.symmetric              = False
    wing.t_tail                 = False
    
    wing.dynamic_pressire_ratop = 1.0
    
    
    
    # Wing Segments
    segment                     = RCAIDE.Components.Wings.Segment()
    segment.tag                 = 'root'
    segment.percent_span_location   = 0.0
    segment.twist               = 0. * Units.deg
    segment.root_chord_percent  = 1.
    segment.dihedral_outboard   = 0 * Units.degrees
    segment.sweeps.quarter_chord= 61.485 * Units.degrees
    segment.thickness_to_chord  = .1
    wing.append_segment(segment)
    
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'segment_1'
    segment.percent_span_location         = 0.2962
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.45
    segment.dihedral_outboard             = 0. * Units.degrees
    segment.sweeps.quarter_chord          = 31.2 * Units.degrees   
    segment.thickness_to_chord            = .1
    wing.append_segment(segment)

    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'segment_2'
    segment.percent_span_location         = 1.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 0.1183 
    segment.dihedral_outboard             = 0.0 * Units.degrees
    segment.sweeps.quarter_chord          = 0.0    
    segment.thickness_to_chord            = .1  
    wing.append_segment(segment)
    
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)
    
    # add to vehicle
    vehicle.append_component(wing)
    
    # ################################################# Fuselage ################################################################
    
    fuselage                                    = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage'
    fuselage.number_coach_seats                 = vehicle.passengers
    fuselage.seats_abreast                      = 6
    fuselage.seat_pitch                         = 1 * Units.meter
    fuselage.finess.nose                        = 1.6
    fuselage.finess.tail                        = 2.
    fuselage.lengths.nose                       = 6.4 * Units.meter
    fuselage.lengths.tail                       = 8.0 * Units.meter
    fuselage.lengths.total                      = 38.02 * Units.meter
    fuselage.lenghts.fore_space                 = 6. * Units.meter
    fuselage.lengths.aft_space                  = 5. * Units.meter
    fuselage.width                              = 3.74 * Units.meter
    fuselage.heights.maximum                    = 3.74 * Units.meter
    fuselage.effective_diamater                 = 
    
    
