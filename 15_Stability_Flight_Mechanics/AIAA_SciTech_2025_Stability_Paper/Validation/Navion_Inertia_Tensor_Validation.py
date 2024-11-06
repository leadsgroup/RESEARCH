# Vehicle.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     

import RCAIDE
from RCAIDE.Framework.Core                              import Units   
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor import design_propeller
from RCAIDE.Library.Methods.Geometry.Planform           import segment_properties
from RCAIDE.Library.Plots                               import *
import numpy as np
import os

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units
from RCAIDE.Library.Plots                 import *     
from RCAIDE.Library.Methods.Weights.Correlation_Buildups import Common
from RCAIDE.Library.Methods.Weights.Moment_of_Inertia.compute_aircraft_moment_of_inertia import compute_aircraft_moment_of_inertia
from RCAIDE.Library.Methods.Weights.Center_of_Gravity     import compute_vehicle_center_of_gravity


# python imports 
import numpy as np  
import matplotlib.pyplot as plt  


def main():
    vehicle = vehicle_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_General_Aviation()
    weight_analysis.vehicle                       = vehicle
    weight_analysis.method                        = 'Raymer'
    #weight_analysis.settings.use_max_fuel_weight  = True 
    results                                       = weight_analysis.evaluate() 
    print("Operating empty weight estimate for Navion: "+str(results))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = vehicle.mass_properties.center_of_gravity
    print("Navion CG location: " + str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI = compute_aircraft_moment_of_inertia(vehicle, CG_location)

    print(MOI)
    sft2     = 1.355817 # 1 slug*ft^2 to 1.355817 kg*m^2
    Navion_true = np.array([[1048.0 , 0, 0], [0, 3000.0, 0], [0, 0, 3530.0]]) * sft2
    error    = (MOI - Navion_true) / Navion_true * 100
    print(error)

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup():

    '''
    This function defines the base vehicle including 
    1) center of gravity (either hard coded or use RCAIDE's built in function)
    2) mass moment of interita (optional)
    
    Key Notes:
    1) The wing that is intended to be the main must be given the tag "main wing". This wing will be used to append 
       a flap and an aileron 
    
    2) If present, the wing that is intended to be the horizontal stabilizer must be given the tag "horizontal_stabilizer" 
       This wing will be used to append an elevator 
    
    3) If present, The wing that is intended to be the  vertical stabilizer must be given the tag "vertical_stabilizer" 
       This wing will be used to append a rudder (optional)
    
    
    '''
    
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
    vehicle.mass_properties.moments_of_inertia.tensor = np.array([[1741,0.0,0.0],[0.0,3759,0.0],[0.0,0.0,4386]])
    vehicle.mass_properties.center_of_gravity         = [[2.239696797,0,-0.131189711 ]]
    vehicle.flight_envelope.ultimate_load             = 5.7
    vehicle.flight_envelope.limit_load                = 3.8
    vehicle.flight_envelope.design_mach_number        = 0.2
    vehicle.design_dynamic_pressure                   = 10e04
    vehicle.reference_area                            = 17.112 
    vehicle.passengers                                = 2
    
    # ------------------------------------------------------------------        
    #   Landing Gear
    # ------------------------------------------------------------------   
    
    landing_gear                          =  RCAIDE.Library.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag                      = "main_landing_gear"
    landing_gear.main_tire_diameter       =  0.3* Units.m
    landing_gear.nose_tire_diameter       =  0.5* Units.m
    landing_gear.main_strut_length        =  0.80* Units.m
    landing_gear.nose_strut_length        =  0.80* Units.m
    landing_gear.main_units               =  2   #number of nose landing gear
    landing_gear.nose_units               =  1   #number of nose landing gear
    landing_gear.main_wheels              =  1  #number of wheels on the main landing gear
    landing_gear.nose_wheels              =  1  #number of wheels on the nose landing gear
    vehicle.landing_gear                  = landing_gear
    
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
    rel_path                              = os.path.dirname(ospath) + separator  + '..' + separator + '..' + separator
     #rel_path                              = os.path.dirname(ospath) + separator + '..' + separator+ '..' + separator+ '..' + separator+ '..' + separator+ '..' + separator

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
  
    # add to vehicle
    vehicle.append_component(wing) 
    

    # ------------------------------------------------------------------        
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------       
    wing                                  = RCAIDE.Library.Components.Wings.Horizontal_Tail()
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
    
    RCAIDE.Library.Methods.Geometry.Planform.wing_planform(wing)     

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------ 
    wing                                  = RCAIDE.Library.Components.Wings.Vertical_Tail()
    wing.tag                              = 'vertical_stabilizer'   
    wing.sweeps.leading_edge              = 20 * Units.degrees 
    wing.thickness_to_chord               = 0.125
    wing.areas.reference                  = 1.163  
    wing.spans.projected                  = 1.4816  
    wing.chords.root                      = 1.2176 
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
    
    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    fuselage = RCAIDE.Library.Components.Fuselages.Fuselage()
    fuselage.tag                                = 'fuselage'
    fuselage.seats_abreast                      = 2
    fuselage.number_coach_seats                 = 4
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
 
    # ------------------------------------------------------------------
    #   Fuel
    # ------------------------------------------------------------------    
    # define fuel weight needed to size fuel system
    fuel                                        = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline()
    fuel.mass_properties                        = RCAIDE.Library.Components.Mass_Properties() 
    fuel.number_of_tanks                        = 1.
    fuel.origin                                 = wing.origin
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density #all of the fuel volume is internal
    fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel.mass_properties.mass                   = 319 *Units.lbs
    vehicle.fuel                                = fuel




    # ########################################################  Energy Network  #########################################################  
    net                                         = RCAIDE.Framework.Networks.Fuel()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    fuel_line                                   = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()   

    #------------------------------------------------------------------------------------------------------------------------------------  
    #  Fuel Tank & Fuel
    #------------------------------------------------------------------------------------------------------------------------------------       
    fuel_tank                                   = RCAIDE.Library.Components.Energy.Sources.Fuel_Tanks.Fuel_Tank()
    fuel_tank.origin                            = wing.origin  
    fuel                                        = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline() 
    fuel.mass_properties.mass                   = 319 *Units.lbs 
    fuel.mass_properties.center_of_gravity      = wing.mass_properties.center_of_gravity
    fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    fuel_tank.fuel                              = fuel     
    fuel_line.fuel_tanks.append(fuel_tank)  
    net.fuel_lines.append(fuel_line)    

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------   
    ice_prop                                   = RCAIDE.Library.Components.Propulsors.ICE_Propeller()     
    ice_prop.active_fuel_tanks                 = ['fuel_tank']
    ice_prop.engine_mass                       = 159
    ice_prop.origin                            = [[0.25, 0, 0]]
                                                     
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
    airfoil.coordinate_file                 =  rel_path + 'Aircraft' + separator + 'Airfoils' + separator + 'NACA_4412.txt'   # absolute path   
    airfoil.polar_files                     =[ rel_path + 'Aircraft' + separator + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_50000.txt',
                                               rel_path + 'Aircraft' + separator + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_100000.txt',
                                               rel_path + 'Aircraft' + separator + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_200000.txt',
                                               rel_path + 'Aircraft' + separator + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_500000.txt',
                                               rel_path + 'Aircraft' + separator + 'Airfoils' + separator + 'Polars' + separator + 'NACA_4412_polar_Re_1000000.txt']  
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

    return vehicle

    
if __name__ == '__main__': 
    main()
    plt.show()    