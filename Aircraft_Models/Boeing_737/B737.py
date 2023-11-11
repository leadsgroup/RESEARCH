# tut_mission_B737.py
# 
# Created:  Aug 2014, RCAIDE Team
# Modified: Jun 2015, RCAIDE Team

""" setup file for a mission with a 737
"""


# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import RCAIDE
from RCAIDE.Core import Units

import numpy as np
import pylab as plt
import pickle
import copy, time
import os 
from RCAIDE.Core import Data, Container 
##import vsp 
#from RCAIDE.Input_Output.OpenVSP.vsp_write import write 

from RCAIDE.Visualization.Performance.Aerodynamics.Vehicle                 import *  
from RCAIDE.Visualization.Performance.Mission                              import *  
from RCAIDE.Visualization.Performance.Aerodynamics.Rotor import *  
from RCAIDE.Visualization.Performance.Energy.Battery                       import *   
from RCAIDE.Visualization.Performance.Noise                                import *  
from RCAIDE.Visualization.Geometry.Three_Dimensional.plot_3d_vehicle       import plot_3d_vehicle 
from RCAIDE.Visualization.Geometry                                         import *
from RCAIDE.Methods.Propulsion.turbofan_sizing                             import turbofan_sizing
from RCAIDE.Methods.Propulsion.turbofan_sizing                             import turbofan_sizing
from RCAIDE.Methods.Geometry.Two_Dimensional.Planform                      import segment_properties
from copy import deepcopy 

from RCAIDE.Methods.Noise.Certification import sideline_noise
from RCAIDE.Methods.Noise.Certification import flyover_noise 
from RCAIDE.Methods.Noise.Certification import approach_noise

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    ti = time.time()
    configs, analyses = full_setup()

    simple_sizing(configs)

    configs.finalize()
    analyses.finalize()

    # weight analysis
    weights = analyses.configs.base.weights
    breakdown = weights.evaluate()      

    # mission analysis
    mission = analyses.missions.base
    results = mission.evaluate() 
    
    # certification calcula 
    #sideline_SPL  = sideline_noise(analyses,configs) 
    #flyover_SPL   = flyover_noise(analyses,configs)  
    #approach_SPL  = approach_noise(analyses,configs)
      
    tf = time.time()
    print ('Time taken: ' + str(round((tf-ti)/60,3))  + ' min')
    
    #save_results(results)
    plot_mission(results)


    return


# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup():

    # vehicle data
    vehicle  = vehicle_setup()
    #write(vehicle, "B737-800")
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(configs_analyses)
    missions_analyses = missions_setup(mission,configs_analyses)

    analyses = RCAIDE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    analyses.missions = missions_analyses

    return configs, analyses

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = RCAIDE.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)

    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Analyses.Weights.Weights_Transport()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Analyses.Aerodynamics.Subsonic_VLM()
    aerodynamics.settings.plot_vortex_distribution  = True     
    aerodynamics.geometry = vehicle
   
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)

    #aerodynamics = RCAIDE.Analyses.Aerodynamics.AVL()
    #aerodynamics.process.compute.lift.inviscid.settings.filenames.avl_bin_name =  '/Users/matthewclarke/Documents/AVL/avl3.35'
    #aerodynamics.settings.print_output = True 
    #aerodynamics.settings.keep_files = True 
    #aerodynamics.geometry = vehicle
    #analyses.append(aerodynamics)
    

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = RCAIDE.Analyses.Stability.Fidelity_Zero()
    stability.geometry = vehicle
    analyses.append(stability)
    
    #stability = RCAIDE.Analyses.Stability.AVL() 
    #stability.settings.filenames.avl_bin_name = '/Users/matthewclarke/Documents/AVL/avl3.35'
    #stability.geometry = vehicle
    #analyses.append(stability)

    # ------------------------------------------------------------------
    #  Noise Analysis 
    noise                   = RCAIDE.Analyses.Noise.Fidelity_Zero()  
    noise.geometry          = vehicle 
    analyses.append(noise)
    
    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Analyses.Energy.Energy()
    energy.network = vehicle.networks
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses    

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup():

    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------

    vehicle = RCAIDE.Vehicle()
    vehicle.tag = 'Boeing_737800'

    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------

    # mass properties
    vehicle.mass_properties.max_takeoff               = 79015.8   # kg
    vehicle.mass_properties.takeoff                   = 79015.8   # kg
    vehicle.mass_properties.operating_empty           = 62746.4   # kg
    vehicle.mass_properties.takeoff                   = 79015.8   # kg
    vehicle.mass_properties.max_zero_fuel             = 62732.0   # kg
    vehicle.mass_properties.cargo                     = 10000.  * Units.kilogram
    vehicle.mass_properties.center_of_gravity         = [[ 15.30987849,   0.        ,  -0.48023939]]
    vehicle.mass_properties.moments_of_inertia.tensor = [[3173074.17, 0 , 28752.77565],[0 , 3019041.443, 0],[0, 0, 5730017.433]] # estimated, not correct
    vehicle.design_mach_number                        = 0.78
    vehicle.design_range                              = 3582 * Units.miles
    vehicle.design_cruise_alt                         = 35000.0 * Units.ft

    # envelope properties
    vehicle.envelope.ultimate_load = 3.75
    vehicle.envelope.limit_load    = 1.5

    # basic parameters
    vehicle.reference_area         = 124.862
    vehicle.passengers             = 170
    vehicle.systems.control        = "fully powered"
    vehicle.systems.accessories    = "medium range"
  
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------

    wing = RCAIDE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'

    wing.aspect_ratio            = 10.18
    wing.sweeps.quarter_chord    = 25 * Units.deg
    wing.thickness_to_chord      = 0.1
    wing.taper                   = 0.1

    wing.spans.projected         = 34.32

    wing.chords.root             = 7.760 * Units.meter
    wing.chords.tip              = 0.782 * Units.meter
    wing.chords.mean_aerodynamic = 4.235 * Units.meter

    wing.areas.reference         = 124.862
    wing.areas.wetted            = 225.08
    
    wing.twists.root             = 4.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees

    wing.origin                  = [[13.61,0,-0.93]]
    wing.aerodynamic_center      = [0,0,0]   

    wing.vertical                = False
    wing.symmetric               = True
    wing.high_lift               = True

    wing.dynamic_pressure_ratio  = 1.0


    # Wing Segments
    root_airfoil                          = RCAIDE.Components.Airfoils.Airfoil()
    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = os.path.dirname(ospath) + separator   
    root_airfoil.coordinate_file          = '..' + separator + 'Airfoils' + separator + 'B737a.txt'
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'Root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 4. * Units.deg
    segment.root_chord_percent            = 1.
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 2.5 * Units.degrees
    segment.sweeps.quarter_chord          = 28.225 * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(root_airfoil)
    wing.append_segment(segment)

    yehudi_airfoil                        = RCAIDE.Components.Airfoils.Airfoil()
    yehudi_airfoil.coordinate_file        = '..' + separator + 'Airfoils' + separator + 'B737b.txt'
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'Yehudi'
    segment.percent_span_location         = 0.324
    segment.twist                         = 0.047193 * Units.deg
    segment.root_chord_percent            = 0.5
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 5.5 * Units.degrees
    segment.sweeps.quarter_chord          = 25. * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(yehudi_airfoil)
    wing.append_segment(segment)

    mid_airfoil                           = RCAIDE.Components.Airfoils.Airfoil()
    mid_airfoil.coordinate_file           = '..' + separator + 'Airfoils' + separator + 'B737c.txt'
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'Section_2'
    segment.percent_span_location         = 0.963
    segment.twist                         = 0.00258 * Units.deg
    segment.root_chord_percent            = 0.220
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 5.5 * Units.degrees
    segment.sweeps.quarter_chord          = 56.75 * Units.degrees
    segment.thickness_to_chord            = .1
    segment.append_airfoil(mid_airfoil)
    wing.append_segment(segment)

    tip_airfoil                           =  RCAIDE.Components.Airfoils.Airfoil()
    tip_airfoil.coordinate_file           = '..' + separator + 'Airfoils' + separator + 'B737d.txt'
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'Tip'
    segment.percent_span_location         = 1.
    segment.twist                         = 0. * Units.degrees
    segment.root_chord_percent            = 0.10077
    segment.thickness_to_chord            = 0.1
    segment.dihedral_outboard             = 0.
    segment.sweeps.quarter_chord          = 0.
    segment.thickness_to_chord            = .1
    segment.append_airfoil(tip_airfoil)
    wing.append_segment(segment)
    
    # Fill out more segment properties automatically
    wing = segment_properties(wing)    

    # control surfaces -------------------------------------------
    slat                          = RCAIDE.Components.Wings.Control_Surfaces.Slat()
    slat.tag                      = 'slat'
    slat.span_fraction_start      = 0.2
    slat.span_fraction_end        = 0.963
    slat.deflection               = 0.0 * Units.degrees
    slat.chord_fraction           = 0.075
    wing.append_control_surface(slat)

    flap                          = RCAIDE.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.7
    flap.deflection               = 0.0 * Units.degrees
    flap.configuration_type       = 'double_slotted'
    flap.chord_fraction           = 0.30
    wing.append_control_surface(flap)

    aileron                       = RCAIDE.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.963
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.16
    wing.append_control_surface(aileron)
    


    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

    wing = RCAIDE.Components.Wings.Horizontal_Tail()
    wing.tag = 'horizontal_stabilizer'

    wing.aspect_ratio            = 4.99
    wing.sweeps.quarter_chord    = 28.2250 * Units.deg  
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.3333 

    wing.spans.projected         = 14.4

    wing.chords.root             = 4.2731 
    wing.chords.tip              = 1.4243 
    wing.chords.mean_aerodynamic = 8.0

    wing.areas.reference         = 41.49
    wing.areas.exposed           = 59.354    # Exposed area of the horizontal tail
    wing.areas.wetted            = 71.81     # Wetted area of the horizontal tail
    wing.twists.root             = 3.0 * Units.degrees
    wing.twists.tip              = 3.0 * Units.degrees

    wing.origin                  = [[33.02,0,1.466]]
    wing.aerodynamic_center      = [0,0,0]

    wing.vertical                = False
    wing.symmetric               = True

    wing.dynamic_pressure_ratio  = 0.9


    # Wing Segments
    segment                        = RCAIDE.Components.Wings.Segment()
    segment.tag                    = 'root_segment'
    segment.percent_span_location  = 0.0
    segment.twist                  = 0. * Units.deg
    segment.root_chord_percent     = 1.0
    segment.dihedral_outboard      = 8.63 * Units.degrees
    segment.sweeps.quarter_chord   = 28.2250  * Units.degrees 
    segment.thickness_to_chord     = .1
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
    elevator                       = RCAIDE.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.09
    elevator.span_fraction_end     = 0.92
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    wing.append_control_surface(elevator)

    # add to vehicle
    vehicle.append_component(wing)


    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing = RCAIDE.Components.Wings.Vertical_Tail()
    wing.tag = 'vertical_stabilizer'

    wing.aspect_ratio            = 1.98865
    wing.sweeps.quarter_chord    = 31.2  * Units.deg   
    wing.thickness_to_chord      = 0.08
    wing.taper                   = 0.1183

    wing.spans.projected         = 8.33
    wing.total_length            = wing.spans.projected 
    
    wing.chords.root             = 10.1 
    wing.chords.tip              = 1.20 
    wing.chords.mean_aerodynamic = 4.0

    wing.areas.reference         = 34.89
    wing.areas.wetted            = 57.25 
    
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees

    wing.origin                  = [[26.944,0,1.54]]
    wing.aerodynamic_center      = [0,0,0]

    wing.vertical                = True
    wing.symmetric               = False
    wing.t_tail                  = False

    wing.dynamic_pressure_ratio  = 1.0


    # Wing Segments
    segment                               = RCAIDE.Components.Wings.Segment()
    segment.tag                           = 'root'
    segment.percent_span_location         = 0.0
    segment.twist                         = 0. * Units.deg
    segment.root_chord_percent            = 1.
    segment.dihedral_outboard             = 0 * Units.degrees
    segment.sweeps.quarter_chord          = 61.485 * Units.degrees  
    segment.thickness_to_chord            = .1
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


    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------

    fuselage = RCAIDE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage'

    fuselage.number_coach_seats    = vehicle.passengers
    fuselage.seats_abreast         = 6
    fuselage.seat_pitch            = 31. * Units.inches
    fuselage.fineness.nose         = 1.6
    fuselage.fineness.tail         = 2.

    fuselage.lengths.nose          = 6.4
    fuselage.lengths.tail          = 8.0
    fuselage.lengths.cabin         = 28.85  
    fuselage.lengths.total         = 38.02  
    fuselage.lengths.fore_space    = 6.
    fuselage.lengths.aft_space     = 5.

    fuselage.width                 = 3.74  

    fuselage.heights.maximum       = 3.74 
    fuselage.heights.at_quarter_length          = 3.74 
    fuselage.heights.at_three_quarters_length   = 3.65 
    fuselage.heights.at_wing_root_quarter_chord = 3.74 

    fuselage.areas.side_projected  = 142.1948 
    fuselage.areas.wetted          = 385.51
    fuselage.areas.front_projected = 12.57

    fuselage.effective_diameter    = 3.74 

    fuselage.differential_pressure = 5.0e4 * Units.pascal # Maximum differential pressure
    
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
    fuselage.Segments.append(segment)   
    
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
    
    # ------------------------------------------------------------------
    #   Nacelles
    # ------------------------------------------------------------------ 
    nacelle                               = RCAIDE.Components.Nacelles.Nacelle()
    nacelle.tag                           = 'nacelle_1'
    nacelle.length                        = 2.71
    nacelle.inlet_diameter                = 1.90
    nacelle.diameter                      = 2.05
    nacelle.areas.wetted                  = 1.1*np.pi*nacelle.diameter*nacelle.length
    nacelle.origin                        = [[13.72, -4.86,-1.9]]
    nacelle.flow_through                  = True   
    nacelle.Airfoil.tag                   = 'nacelle_naca_airfoil'    
    nacelle.Airfoil.naca_4_series_airfoil = '2410' 

    nacelle_2                     = deepcopy(nacelle)
    nacelle_2.tag                 = 'nacelle_2'
    nacelle_2.origin              = [[13.72, 4.86,-1.9]]
    
    vehicle.append_component(nacelle)  
    vehicle.append_component(nacelle_2)     
    
    
    # ------------------------------------------------------------------
    #   Turbofan Network
    # ------------------------------------------------------------------

    #instantiate the gas turbine network
    turbofan = RCAIDE.Components.Energy.Networks.Turbofan()
    turbofan.tag = 'turbofan'

    # setup
    turbofan.number_of_engines = 2.0
    turbofan.bypass_ratio      = 5.4
    turbofan.engine_length     = 2.71

    # This origin is overwritten by compute_component_centers_of_gravity(base,compute_propulsor_origin=True)
    turbofan.origin            = [[13.72, 4.86,-1.9],[13.72, -4.86,-1.9]]

    # working fluid
    turbofan.working_fluid = RCAIDE.Attributes.Gases.Air()


    # ------------------------------------------------------------------
    #   Component 1 - Ram

    # to convert freestream static to stagnation quantities

    # instantiate
    ram = RCAIDE.Components.Energy.Converters.Ram()
    ram.tag = 'ram'

    # add to the network
    turbofan.append(ram)


    # ------------------------------------------------------------------
    #  Component 2 - Inlet Nozzle

    # instantiate
    inlet_nozzle = RCAIDE.Components.Energy.Converters.Compression_Nozzle()
    inlet_nozzle.tag = 'inlet_nozzle'

    # setup
    inlet_nozzle.polytropic_efficiency = 0.98
    inlet_nozzle.pressure_ratio        = 0.98

    # add to network
    turbofan.append(inlet_nozzle)


    # ------------------------------------------------------------------
    #  Component 3 - Low Pressure Compressor

    # instantiate
    compressor = RCAIDE.Components.Energy.Converters.Compressor()
    compressor.tag = 'low_pressure_compressor'

    # setup
    compressor.polytropic_efficiency = 0.91
    compressor.pressure_ratio        = 1.14

    # add to network
    turbofan.append(compressor)


    # ------------------------------------------------------------------
    #  Component 4 - High Pressure Compressor

    # instantiate
    compressor = RCAIDE.Components.Energy.Converters.Compressor()
    compressor.tag = 'high_pressure_compressor'

    # setup
    compressor.polytropic_efficiency = 0.91
    compressor.pressure_ratio        = 13.415

    # add to network
    turbofan.append(compressor)


    # ------------------------------------------------------------------
    #  Component 5 - Low Pressure Turbine

    # instantiate
    turbine = RCAIDE.Components.Energy.Converters.Turbine()
    turbine.tag='low_pressure_turbine'

    # setup
    turbine.mechanical_efficiency = 0.99
    turbine.polytropic_efficiency = 0.93

    # add to network
    turbofan.append(turbine)


    # ------------------------------------------------------------------
    #  Component 6 - High Pressure Turbine

    # instantiate
    turbine = RCAIDE.Components.Energy.Converters.Turbine()
    turbine.tag='high_pressure_turbine'

    # setup
    turbine.mechanical_efficiency = 0.99
    turbine.polytropic_efficiency = 0.93

    # add to network
    turbofan.append(turbine)


    # ------------------------------------------------------------------
    #  Component 7 - Combustor

    # instantiate
    combustor = RCAIDE.Components.Energy.Converters.Combustor()
    combustor.tag = 'combustor'

    # setup
    combustor.efficiency                = 0.99
    combustor.alphac                    = 1.0
    combustor.turbine_inlet_temperature = 1450
    combustor.pressure_ratio            = 0.95
    combustor.fuel_data                 = RCAIDE.Attributes.Propellants.Jet_A()

    # add to network
    turbofan.append(combustor)


    # ------------------------------------------------------------------
    #  Component 8 - Core Nozzle

    # instantiate
    nozzle = RCAIDE.Components.Energy.Converters.Expansion_Nozzle()
    nozzle.tag = 'core_nozzle'

    # setup
    nozzle.polytropic_efficiency = 0.95
    nozzle.pressure_ratio        = 0.99

    # add to network
    turbofan.append(nozzle)


    # ------------------------------------------------------------------
    #  Component 9 - Fan Nozzle

    # instantiate
    nozzle = RCAIDE.Components.Energy.Converters.Expansion_Nozzle()
    nozzle.tag = 'fan_nozzle'

    # setup
    nozzle.polytropic_efficiency = 0.95
    nozzle.pressure_ratio        = 0.99

    # add to network
    turbofan.append(nozzle)


    # ------------------------------------------------------------------
    #  Component 10 - Fan

    # instantiate
    fan = RCAIDE.Components.Energy.Converters.Fan()
    fan.tag = 'fan'

    # setup
    fan.polytropic_efficiency = 0.93
    fan.pressure_ratio        = 1.7

    # add to network
    turbofan.append(fan)


    # ------------------------------------------------------------------
    #Component 10 : thrust (to compute the thrust)
    thrust = RCAIDE.Components.Energy.Processes.Thrust()
    thrust.tag ='compute_thrust'

    #total design thrust (includes all the engines)
    thrust.total_design             = 2*24000. * Units.N #Newtons

    #design sizing conditions
    altitude      = 35000.0*Units.ft
    mach_number   = 0.78
    isa_deviation = 0.

    #Engine setup for noise module


    # add to network
    turbofan.thrust = thrust

    turbofan.core_nozzle_diameter = 0.92
    turbofan.fan_nozzle_diameter  = 1.659
    turbofan.engine_height        = 0.5  #Engine centerline heigh above the ground plane
    turbofan.exa                  = 1    #distance from fan face to fan exit/ fan diameter)
    turbofan.plug_diameter        = 0.1  #dimater of the engine plug
    turbofan.geometry_xe          = 1. # Geometry information for the installation effects function
    turbofan.geometry_ye          = 1. # Geometry information for the installation effects function
    turbofan.geometry_Ce          = 2. # Geometry information for the installation effects function





    #size the turbofan
    turbofan_sizing(turbofan,mach_number,altitude)

    # add  gas turbine network turbofan to the vehicle
    vehicle.append_component(turbofan)

    # ------------------------------------------------------------------
    #  Fuel
    # ------------------------------------------------------------------
    fuel                                  = RCAIDE.Components.Physical_Component()
    vehicle.fuel                          = fuel
    fuel.mass_properties.mass             = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    fuel.origin                           = vehicle.wings.main_wing.mass_properties.center_of_gravity
    fuel.mass_properties.center_of_gravity= vehicle.wings.main_wing.aerodynamic_center

    # ------------------------------------------------------------------
    #  Landing Gear
    # ------------------------------------------------------------------
    landing_gear                          = RCAIDE.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag                      = "main_landing_gear"
    landing_gear.main_tire_diameter       = 1.12000 * Units.m
    landing_gear.nose_tire_diameter       = 0.6858 * Units.m
    landing_gear.main_strut_length        = 1.8 * Units.m
    landing_gear.nose_strut_length        = 1.3 * Units.m
    landing_gear.main_units               = 1    #number of nose landing gear
    landing_gear.nose_units               = 1    #number of nose landing gear
    landing_gear.main_wheels              = 2    #number of wheels on the main landing gear
    landing_gear.nose_wheels              = 2    #number of wheels on the nose landing gear
    vehicle.landing_gear                  = landing_gear

    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------
 
    # plot vehicle 
    #plot_3d_vehicle(vehicle) 
    #plt.show()  


    return vehicle


# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):

    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------
    configs = RCAIDE.Components.Configs.Config.Container()

    base_config = RCAIDE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    base_config.landing_gear.gear_condition                               = 'up'
    configs.append(base_config)

    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag = 'cruise' 
    configs.append(config)
    config.wings['main_wing'].control_surfaces.flap.deflection       = 0. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection       = 0. * Units.deg

    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag                                                       = 'takeoff'
    config.wings['main_wing'].control_surfaces.flap.deflection       = 20. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection       = 25. * Units.deg
    #Noise input for the landing gear                                
    config.landing_gear.gear_condition                               = 'up'       
    config.output_filename                                           = 'Flyover_'

    config.networks['turbofan'].fan.rotation            = 3470. #N1 speed
    config.networks['turbofan'].fan_nozzle.noise_speed  = 315.
    config.networks['turbofan'].core_nozzle.noise_speed = 415.

    configs.append(config)

    # ------------------------------------------------------------------
    #   Cutback Configuration
    # ------------------------------------------------------------------
    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag                                                       = 'cutback'
    config.wings['main_wing'].control_surfaces.flap.deflection       = 20. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection       = 20. * Units.deg
    #Noise input for the landing gear                                
    config.landing_gear.gear_condition                               = 'up'       
    config.output_filename                                           = 'Cutback_'

    config.networks['turbofan'].fan.rotation            = 2780. #N1 speed
    config.networks['turbofan'].fan_nozzle.noise_speed  = 210.
    config.networks['turbofan'].core_nozzle.noise_speed = 360.

    configs.append(config)

    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag = 'landing'

    config.wings['main_wing'].control_surfaces.flap.deflection       = 30. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection       = 25. * Units.deg  
    #Noise input for the landing gear                              
    config.landing_gear.gear_condition                               = 'down'    
    config.output_filename                                           = 'Approach_'

    config.networks['turbofan'].fan.rotation     = 2030.  #N1 speed
    config.networks['turbofan'].fan_nozzle.noise_speed  = 109.3
    config.networks['turbofan'].core_nozzle.noise_speed = 92.

    configs.append(config)

    # ------------------------------------------------------------------
    #   Short Field Takeoff Configuration
    # ------------------------------------------------------------------

    config = RCAIDE.Components.Configs.Config(base_config)
    config.tag = 'short_field_takeoff'

    config.wings['main_wing'].control_surfaces.flap.deflection       = 20. * Units.deg
    config.wings['main_wing'].control_surfaces.slat.deflection       = 20. * Units.deg
    configs.append(config)

    return configs

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------
def plot_mission(results,line_style='bo-'):

    
    # Plot Flight Conditions 
    plot_flight_conditions(results, line_style)
    
    # Plot Aerodynamic Forces 
    plot_aerodynamic_forces(results, line_style)
    
    # Plot Aerodynamic Coefficients 
    plot_aerodynamic_coefficients(results, line_style)
    
    # Plot Static Stability Coefficients 
    plot_stability_coefficients(results, line_style)    
    
    # Drag Components
    plot_drag_components(results, line_style)
    
    # Plot Altitude, sfc, vehicle weight 
    plot_altitude_sfc_weight(results, line_style)
    
    # Plot Velocities 
    plot_aircraft_velocities(results, line_style)  
    
    # Plot Trajectory
    plot_flight_trajectory(results)

    return 

# ----------------------------------------------------------------------
#   Sizing for the Vehicle Configs
# ----------------------------------------------------------------------
def simple_sizing(configs):

    base = configs.base
    base.pull_base()

    # zero fuel weight
    base.mass_properties.max_zero_fuel = 0.9 * base.mass_properties.max_takeoff 

    # wing areas
    for wing in base.wings: 
        wing = RCAIDE.Methods.Geometry.Two_Dimensional.Planform.wing_planform(wing) 
        wing.areas.wetted   = 2.0 * wing.areas.reference
        wing.areas.exposed  = 0.8 * wing.areas.wetted
        wing.areas.affected = 0.6 * wing.areas.wetted
        
    # fuselage seats
    base.fuselages['fuselage'].number_coach_seats = base.passengers

    # diff the new data
    base.store_diff()

    # done!
    return

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'base_mission'

    #airport
    airport = RCAIDE.Attributes.Airports.Airport()
    airport.altitude   =  0.0  * Units.ft
    airport.delta_isa  =  0.0
    airport.atmosphere = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976()

    mission.airport = airport    

    # unpack Segments module
    Segments = RCAIDE.Analyses.Mission.Segments

    # base segment
    base_segment = Segments.Segment()


    # ------------------------------------------------------------------
    #   First Climb Segment: constant Mach, constant segment angle 
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1"

    segment.analyses.extend( analyses.takeoff )
    
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg      

    segment.altitude_start = 0.001   * Units.km
    segment.altitude_end   = 3.0   * Units.km
    segment.air_speed      = 125.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']

    # add to misison
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Second Climb Segment: constant Speed, constant segment angle 
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2"

    segment.analyses.extend( analyses.cutback )

    segment.altitude_end   = 8.0   * Units.km
    segment.air_speed      = 190.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']
    segment.true_course_angle    = 20   * Units.degrees 

    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Climb Segment: constant Mach, constant segment angle 
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3"

    segment.analyses.extend( analyses.base )

    segment.altitude_end = 10.668 * Units.km
    segment.air_speed    = 226.0  * Units['m/s']
    segment.climb_rate   = 3.0    * Units['m/s']
    segment.true_course_angle      = 10   * Units.degrees 

    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------    
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise"

    segment.analyses.extend( analyses.base )
    segment.altitude   = 35000*Units.feet
    segment.air_speed  = 230.412 * Units['m/s']
    segment.distance   =(3933.65 + 770 - 92.6) * Units.km
    segment.true_course_angle    = 36.87  * Units.degrees 
    segment.true_course_angle    = 15   * Units.degrees 
    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   First Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1"

    segment.analyses.extend( analyses.base )

    segment.altitude_end = 8.0   * Units.km
    segment.air_speed    = 220.0 * Units['m/s']
    segment.descent_rate = 4.5   * Units['m/s']
    segment.true_course_angle  = 135   * Units.degrees 

    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Second Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_2"

    segment.analyses.extend( analyses.base )

    segment.altitude_end = 6.0   * Units.km
    segment.air_speed    = 195.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']
    segment.true_course_angle  = 160   * Units.degrees 

    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"

    segment.analyses.extend( analyses.base )

    segment.altitude_end = 4.0   * Units.km
    segment.air_speed    = 170.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']
    segment.true_course_angle      = 180   * Units.degrees 

    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Fourth Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_4"

    segment.analyses.extend( analyses.base )

    segment.altitude_end = 2.0   * Units.km
    segment.air_speed    = 150.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']
    segment.true_course_angle      = 180   * Units.degrees 


    # add to mission
    mission.append_segment(segment)



    # ------------------------------------------------------------------
    #   Fifth Descent Segment: consant speed, constant segment rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_5"

    segment.analyses.extend( analyses.landing)

    segment.altitude_end = 0.0   * Units.km
    segment.air_speed    = 145.0 * Units['m/s']
    segment.descent_rate = 3.0   * Units['m/s']
    segment.true_course_angle      = 180   * Units.degrees 


    # append to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------

    return mission

def missions_setup(base_mission,analyses):

    # the mission container
    missions = RCAIDE.Analyses.Mission.Mission.Container()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------ 
    missions.base = base_mission


    # ------------------------------------------------------------------
    #   Mission for Constrained Fuel
    # ------------------------------------------------------------------    
    fuel_mission           = RCAIDE.Analyses.Mission.Mission() 
    fuel_mission.tag       = 'fuel'
    fuel_mission.range     = 1277. * Units.nautical_mile
    fuel_mission.payload   = 19000.
    missions.append(fuel_mission)    


    # ------------------------------------------------------------------
    #   Mission for Constrained Short Field
    # ------------------------------------------------------------------    
    short_field            = RCAIDE.Analyses.Mission.Mission(base_mission) 
    short_field.tag        = 'short_field'  
    
    airport                = RCAIDE.Attributes.Airports.Airport()
    airport.altitude       =  0.0  * Units.ft
    airport.delta_isa      =  0.0
    airport.atmosphere     = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976()
    airport.available_tofl = 1500.
    short_field.airport    = airport    
    missions.append(short_field)
    
    
    # ------------------------------------------------------------------
    #   Mission for Fixed Payload
    # ------------------------------------------------------------------    
    payload         = RCAIDE.Analyses.Mission.Mission()  
    payload.tag     = 'payload'
    payload.range   = 2316. * Units.nautical_mile
    payload.payload = 19000.
    missions.append(payload)

    
    # ------------------------------------------------------------------
    #   Mission for Takeoff Noise
    # ------------------------------------------------------------------    
    takeoff                           = RCAIDE.Analyses.Mission.Sequential_Segments()
    takeoff.tag                       = 'takeoff'   
                                      
    # airport                          
    airport                           = RCAIDE.Attributes.Airports.Airport()
    airport.altitude                  =  0.0  * Units.ft
    airport.delta_isa                 =  0.0
    airport.atmosphere                = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    takeoff.airport                   = airport    

    # unpack Segments module
    Segments                          = RCAIDE.Analyses.Mission.Segments 
    base_segment                      = Segments.Segment()  
    atmosphere                        = RCAIDE.Attributes.Atmospheres.Earth.US_Standard_1976()
    planet                            = RCAIDE.Attributes.Planets.Earth() 
    
    # Climb Segment: Constant throttle, constant speed
    segment                           = Segments.Climb.Constant_Throttle_Constant_Speed(base_segment)
    segment.tag                       = "climb"    
    segment.analyses.extend(analyses.base ) 
    segment.atmosphere                = atmosphere    
    segment.planet                    = planet 
    segment.altitude_start            =  35. *  Units.fts
    segment.altitude_end              = 304.8 *  Units.meter
    segment.air_speed                 = 85.4 * Units['m/s']
    segment.throttle                  = 1.  
    ones_row                          = segment.state.ones_row  
    takeoff.append_segment(segment)

    # Cutback Segment: Constant speed, constant segment angle
    segment                           = Segments.Climb.Constant_Speed_Constant_Angle_Noise(base_segment)
    segment.tag                       = "cutback"   
    segment.atmosphere                = atmosphere    
    segment.planet                    = planet     
    segment.analyses.extend(analyses.base )
    segment.air_speed                 = 85.4 * Units['m/s']
    segment.climb_angle               = 2.86  * Units.degrees 
    takeoff.append_segment(segment)  
    
    # append mission 
    missions.append(takeoff)

    # ------------------------------------------------------------------
    #   Mission for Sideline Noise
    # ------------------------------------------------------------------     
    sideline_takeoff                  = RCAIDE.Analyses.Mission.Sequential_Segments()
    sideline_takeoff.tag              = 'sideline_takeoff'   
    sideline_takeoff.airport          = airport  
    segment                           = Segments.Climb.Constant_Throttle_Constant_Speed(base_segment)
    segment.tag                       = "climb"    
    segment.analyses.extend(analyses.base)
    segment.atmosphere                = atmosphere    
    segment.planet                    = planet     
    segment.altitude_start            =  35. *  Units.fts
    segment.altitude_end              = 1600 *  Units.fts
    segment.air_speed                 = 85.4 * Units['m/s']
    segment.throttle                  = 1.  
    ones_row                          = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 12. * Units.deg  
    segment.state.unknowns.wind_angle = ones_row(1) * 5. * Units.deg  
    sideline_takeoff.append_segment(segment)   
    
    missions.append(sideline_takeoff)
    
    # -------------------   -----------------------------------------------
    #   Mission for Landing Noise
    # ------------------------------------------------------------------    
    landing                           = RCAIDE.Analyses.Mission.Sequential_Segments()
    landing.tag                       = 'landing'   
    landing.airport                   = airport      
    segment                           = Segments.Descent.Constant_Speed_Constant_Angle(base_segment)
    segment.tag                       = "descent"
    segment.analyses.extend(analyses.base ) 
    segment.atmosphere                = atmosphere    
    segment.planet                    = planet     
    segment.altitude_start            = 2.0   * Units.km
    segment.altitude_end              = 0.
    segment.air_speed                 = 67. * Units['m/s']
    segment.descent_angle             = 3.0   * Units.degrees  
    landing.append_segment(segment)
        
    missions.append(landing)
    
    return missions  


def save_results(results):
 
    # Store data (serialize)
    with open('B737_results.pkl', 'wb') as file:
        pickle.dump(results, file)
        
    return
 
if __name__ == '__main__': 
    main()    
    plt.show()




    