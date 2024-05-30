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
from RCAIDE.Methods.Energy.Propulsion.Converters.Turbofan.design_turbofan   import design_turbofan, size_optimal_motor
from RCAIDE.Methods.Power.Battery.Sizing                import initialize_from_circuit_configuration
from RCAIDE.Methods.Weights.Correlation_Buildups.Propulsion     import nasa_motor

# python imports
import numpy as np
from copy import deepcopy
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
    fuselage.effective_diamater                 = 3.74 * Units.meter
    fuselage.area.side_projected                = 142.1948 * Units['meters*2']
    fuselage.area.wetted                        = 446.718 * Units['meters**2']
    fuselage.area.front_projected               = 12.57 * Units['meters**2']
    fuselage.differential_pressure              = 5.0e4 * Units.pascal
    fuselage.heights.at_quarter_length          = 3.74 * Units.meter
    fuselage.heights.at_three_quarters_length   = 3.65 * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = 3.74 * Units.meter
    
    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_0'
    segment.percent_x_location                  = 0.0000
    segment.percent_z_location                  = -0.00144
    segment.height                              = 0.0100
    segment.width                               = 0.0100
    fuselage.Segments.append(segment)
    
    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_1'
    segment.percent_x_location                  = 0.00576
    segment.percent_z_location                  = -0.00144
    segment.height                              = 0.7500
    segment.width                               = 0.6500
    
    # Segment
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_2'
    segment.percent_x_location                  = 0.02017 
    segment.percent_z_location                  = 0.00000 
    segment.height                              = 1.52783 
    segment.width                               = 1.20043 
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_3'   
    segment.percent_x_location                  = 0.03170 
    segment.percent_z_location                  = 0.00000 
    segment.height                              = 1.96435 
    segment.width                               = 1.52783 
    fuselage.Segments.append(segment)   

    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_4'   
    segment.percent_x_location                  = 0.04899 	
    segment.percent_z_location                  = 0.00431 
    segment.height                              = 2.72826 
    segment.width                               = 1.96435 
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_5'   
    segment.percent_x_location                  = 0.07781 
    segment.percent_z_location                  = 0.00861 
    segment.height                              = 3.49217 
    segment.width                               = 2.61913 
    fuselage.Segments.append(segment)     
    
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_6'   
    segment.percent_x_location                  = 0.10375 
    segment.percent_z_location                  = 0.01005 
    segment.height                              = 3.70130 
    segment.width                               = 3.05565 
    fuselage.Segments.append(segment)             
     
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_7'   
    segment.percent_x_location                  = 0.16427 
    segment.percent_z_location                  = 0.01148 
    segment.height                              = 3.92870 
    segment.width                               = 3.71043 
    fuselage.Segments.append(segment)    
    
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_8'   
    segment.percent_x_location                  = 0.22478 
    segment.percent_z_location                  = 0.01148 
    segment.height                              = 3.92870 
    segment.width                               = 3.92870 
    fuselage.Segments.append(segment)   
    
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_9'     
    segment.percent_x_location                  = 0.69164 
    segment.percent_z_location                  = 0.01292
    segment.height                              = 3.81957
    segment.width                               = 3.81957
    fuselage.Segments.append(segment)     
        
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_10'     
    segment.percent_x_location                  = 0.71758 
    segment.percent_z_location                  = 0.01292
    segment.height                              = 3.81957
    segment.width                               = 3.81957
    fuselage.Segments.append(segment)   
        
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_11'     
    segment.percent_x_location                  = 0.78098 
    segment.percent_z_location                  = 0.01722
    segment.height                              = 3.49217
    segment.width                               = 3.71043
    fuselage.Segments.append(segment)    
        
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_12'     
    segment.percent_x_location                  = 0.85303
    segment.percent_z_location                  = 0.02296
    segment.height                              = 3.05565
    segment.width                               = 3.16478
    fuselage.Segments.append(segment)             
        
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_13'     
    segment.percent_x_location                  = 0.91931 
    segment.percent_z_location                  = 0.03157
    segment.height                              = 2.40087
    segment.width                               = 1.96435
    fuselage.Segments.append(segment)               
        
    # Segment                                   
    segment                                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    segment.tag                                 = 'segment_14'     
    segment.percent_x_location                  = 1.00 
    segment.percent_z_location                  = 0.04593
    segment.height                              = 1.09130
    segment.width                               = 0.21826
    fuselage.Segments.append(segment)
    
    # add to vehicle
    vehicle.append_component(fuselage)
    
    
    # ################################################# Nacelles ##################################################################        

    # -----------------------------------------------------------------
    # Design the Nacelle
    # -----------------------------------------------------------------
    nacelle                         = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.diameter                 = 2.05
    nacelle.length                  = 2.71
    nacelle.tag                     = 'nacelle_1'
    nacelle.inlet_diameter          = 2.0
    nacelle.origin                  = [[13.5,4.38,-1.1]]
    Awet                            = 1.1*np.pi*nacelle.diameter*nacelle.length # 1.1 is simple coefficient
    nacelle.area.wetted             = Awet
    nacelle.Airfoil.NACA_4_series_flag  = True
    nacelle.Airfoil.coordinate_file = '2410'
    
    nac_segment                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                 = 'segment_1'
    nac_segment.percent_x_location  = 0.0
    nac_segment.height              = 2.05
    nac_segment.width               = 2.05
    nacelle.append_segment(nac_segment)
    
    nac_segment                     = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                 = 'segment_2'
    nac_segment.percent_x_location  = 0.3
    nac_segment.height              = 2.1
    nac_segment.width               = 2.1
    nacelle.append_segment(nac_segment)
    
    nac_segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                       = 'segment_3'
    nac_segment.percent_x_location        = 0.4  
    nac_segment.height                    = 2.05
    nac_segment.width                     = 2.05 
    nacelle.append_segment(nac_segment)  
     
    nac_segment                           = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                       = 'segment_4'
    nac_segment.percent_x_location        = 0.75  
    nac_segment.height                    = 1.9
    nac_segment.width                     = 1.9
    nacelle.append_segment(nac_segment)  
    
    nac_segment                          = RCAIDE.Components.Lofted_Body_Segment.Segment()
    nac_segment.tag                      = 'segment_5'
    nac_segment.percent_x_location       = 1.0
    nac_segment.height                   = 1.7 
    nac_segment.width                    = 1.7
    nacelle.append_segment(nac_segment)    
        
    nacelle_2                             = deepcopy(nacelle)
    nacelle_2.tag                         = 'nacelle_2'
    nacelle_2.origin                      = [[13.5,-4.38,-1.1]]
    
    vehicle.append_component(nacelle)   
    vehicle.append_component(nacelle_2)    
    
    # ################################################# Energy Network #######################################################         
    #------------------------------------------------------------------------------------------------------------------------- 
    #  Turbofan Network
    #-------------------------------------------------------------------------------------------------------------------------     
    net                         = RCAIDE.Energy.Networks.Series_Parallel_Hybrid_Turboelectric_Engine()
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # HPCU(Hybrid Power Control Unit)
    #------------------------------------------------------------------------------------------------------------------------------------
    HPCU                        = RCAIDE.Energy.Distributors.Hybrid_Power_Control_Unit()
    HPCU.fixed_voltage          = True
    
    #------------------------------------------------------------------------------------------------------------------------- 
    #   Fuel
    #-------------------------------------------------------------------------------------------------------------------------
    # fuel tank
    fuel_tank                   = RACIDE.Energy.Storages.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin            = wing.origin
    
    #fuel
    fuel                        = RCAIDE.Attributes.Propellants.Aviation_Gasoline()
    fuel.mass_properties.mass   = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    fuel.origin                 = vehicle.wings.main_wing.mass_properties.center_of_gravity
    fuel.mass_properties.center_of_gravity  = vehicle.wings.main_wing.aerodynamic_center
    fuel.internal_volume        = fuel.mass_properties.mass/fuel.density
    fuel_tank.fuel              = fuel
    HPCU.fuel_tanks.append(fuel_tank)
    
    turbofan                    = RCAIDE.Energy.Converters.Turbofan()
    turbofan.tag                = 'modified_pratt_whitney_jt9d'
    turbofan.origin             = [[13,72, 4.86, -1.1]]
    turbofan.engine_length      = 2.71
    turbofan.bypass_ratio       = 5.4
    turbofan.design_altitude    = 35000.0*Units.ft
    turbofan.design_mach_number = 0.78
    turbofan.design_thrust      = 35000.0* Units.N
    
    # fan
    fan                         = RCAIDE.Energy.Converters.Fan()
    fan.tag                     = 'fan'
    fan.polytropic_efficiency   = 0.93
    fan.pressure_ratio          = 1.7
    fan.diameter                = 2.05
    fan.design_angular_velocity = 3750*Units.rpm
    turbofan.fan                = fan
    
    # working fluid
    turbofan.working_fluid      = RCAIDE.Attributes.Gases.Air()
    ram                         = RCAIDE.Energy.Converters.Ram()
    ram.tag                     = 'ram'
    turbofan.ram                = ram
    
    # inlet nozzle
    inlet_nozzle                = RCAIDE.Energy.Converters.Nozzle()
    inlet_nozzle.tag            = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency  = 0.98
    inlet_nozzle.pressure_ratio = 0.98
    turbofan.inlet_nozzle       = inlet_nozzle
    
    
    # low pressure compressor
    low_pressure_compressor     = RCAIDE.Energy.Converters.Compressor()
    low_pressure_compressor.tag = 'lpc'
    low_pressure_compressor.polytropic_efficiency   = 0.91
    low_pressure_compressor.pressure_ratio  = 1.9
    turbofan.low_pressure_compressor    = low_pressure_compressor
    
    # high pressure compressor
    high_pressure_compressor    = RCAIDE.Energy.Converters.Compressor()
    high_pressure_compressor.tag    = 'hpc'
    high_pressure_compressor.polytropic_efficiency  = 0.91
    high_pressure_compressor.pressure_ratio = 10.0
    turbofan.high_pressure_compressor   = high_pressure_compressor
    
    # low pressure turbine
    low_pressure_turbine        = RCAIDE.Energy.Converters.Turbine()
    low_pressure_turbine.tag    = 'lpt'
    low_pressure_turbine.mechanical_efficiency  = 0.99
    low_pressure_turbine.polytropic_efficiency  = 0.93
    turbofan.low_pressure_turbine   = low_pressure_turbine
    
    # high pressure turbine
    high_pressure_turbine       = RCAIDE.Energy.Converters.Turbine()
    high_pressure_turbine.tag   = 'hpt'
    high_pressure_turbine.mechanical_efficiency = 0.99
    high_pressure_turbine.polytropic_efficiency = 0.93
    turbofan.high_pressure_turbine  = high_pressure_turbine
    
    # combustor
    combustor                   = RCAIDE.Energy.Converters.Combustor()
    combustor.tag               = 'Comb'
    combustor.efficiency        = 0.99
    combustor.alphc             = 1.0
    combustor.turbine_inlet_temperature = 1500
    combustor.pressure_ratio    = 0.95
    combustor.fuel_data         = RCAIDE.Attributes.Propellants.Jet_A()
    turbofan.combustor          = combustor
    
    # core nozzle
    core_nozzle                 = RCAIDE.Energy.Converters.Expansion_Nozzle()
    core_nozzle.tag             = 'core nozzle'
    core_nozzle.polytropic_efficiency   = 0.95
    core_nozzle.pressure_ratio  = 0.99
    turbofan.core_nozzle        = core_nozzle
    
    # fan nozzle
    fan_nozzle                  = RCAIDE.Energy.Converters.Expansion_Nozzle()
    fan_nozzle.tag              = 'fan nozzle'
    fan_nozzle.polytropic_efficiency    = 0.95
    fan_nozzle.pressure_ratio   = 0.99
    turbofan.fan_nozzle         = fan_nozzle
    
    # design turbofan(putting together)
    design_turbofan(turbofan)
    
    HPCU.turbofans.append(turbofan)
    
    turbofan_2          = deepcopy(turbofan)
    turbofan_2.tag      = 'preopeller_2'
    turbofan_2.origin   = [[13.72, -4.86, -1.1]]
    HPCU.turbofans.append(turbofan_2)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Electronic Speed Controller    
    #------------------------------------------------------------------------------------------------------------------------------------
    esc_1           = RCAIDE.Library.Components.Propulsors.Modulators.Electronic_Speed_Controller()
    esc_1.tag       = 'esc_1'
    esc_1.efficiency    = 0.95
    HPCU.electronic_speed_controllers.append(esc_1)
    
    esc_2           = RCAIDE.Library.Components.Propulsors.Modulators.Electronic_Speed_Controller()
    esc_2.tag       = 'esc_2'
    esc_2.efficiency    = 0.95
    HPCU.electronic_speed_controllers.append(esc_2)
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Electronic Speed Controller    
    #------------------------------------------------------------------------------------------------------------------------------------    
    bat                                             = RCAIDE.Library.Components.Energy.Batteries.Lithium_Ion_NMC()
    bat.pack.electrical_configuration.series        = 140
    bat.pack.electrical_configuration.parallel      = 100
    initialize_from_circuit_configuration(bat)
    bat.pack.number_of_modules                    = 14
    bat.module.geometric_configuration.total        = bat.pack.electrical_configuration.total
    bat.module.voltage                              = bat.pack.maximum_voltage/bat.pack.number_of_modules
    bat.module.geomtric_configuration.normal_count  = 24
    bat.module.geomtric_configuration.prallel_count = 40
    bat.thermal_management_system                   = RCAIDE.Energy.Thermal_Management.Batteries.Atmospheric_Air_Convection_Heat_Exchanger()
    HPCU.voltage                                    = bat.pack.maximum_voltage
    HPCU.batteries.append(bat)
    
    #------------------------------------------------------------------------------------------------------------------------------------           
    # Motors 
    #------------------------------------------------------------------------------------------------------------------------------------
    motor                   = RCAIDE.Library.Components.Propulsors.Converters.DC_Motor()
    motor.efficiency        = 0.98
    motor.nominal_voltage   = bat.pack.maximum_voltage
    motor.no_load_current   = 1
    motor.rotor_radius      = 2.05
    motor.design_torque     = 93135 # HP*5252/(fan.design_angular_velocity)
    motor.angular_velocity  = fan.design_angular_velocity
    motor                   = size_optimal_motor(motor)
    motor.mass_properties.mass  = nasa_motor(motor.design_torque)
    HPCU.motors.append(motor)
    motor_2                 = deepcopy(motor)
    motor_2.orgin           = [[2., -2.5, 0.95]]
    HPCU.motors.append(motor_2)
    
    # append bus
    net.hybrid_power_control_units.append(HPCU)
    
    
    vehicle.append_energy_network(net)
    
    return vehicle

    
    

    
    
    
    