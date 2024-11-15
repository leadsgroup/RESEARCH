# Optimize.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------  
import RCAIDE
from   RCAIDE.Framework.Core                        import Units, Data 
from   RCAIDE.Framework.Optimization.Packages.scipy import scipy_setup
from   RCAIDE.Framework.Optimization.Common         import Nexus
from   Optimization_Data_2_CSV                      import Optimization_Data_2_CSV

# python imports 
import numpy as np
import time
import  pickle

# local imports 
import Vehicles
import Missions
import Analyses
import Procedure
import matplotlib.pyplot as plt 
# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def size_control_surfaces(configuration_CGbatt_MOIbatt,vehicle, cruise_velocity = 120 * Units['mph'], cruise_altitude= 5000*Units.feet): 
    '''
    STICK FIXED (STATIC STABILITY AND DRAG OTIMIZATION
    '''
    ti_0 = time.time()
    ti   = time.time()  
    solver_name       = 'SLSQP' 
    planform_optimization_problem = stick_fixed_stability_and_drag_optimization_setup(vehicle,cruise_velocity,cruise_altitude)
    output_stick_fixed = scipy_setup.SciPy_Solve(planform_optimization_problem,solver=solver_name, sense_step = 1E-2, tolerance = 1E-2)  
    print (output_stick_fixed)    
    tf           = time.time()
    elapsed_time_stick_fixed = round((tf-ti)/60,2)
    print('Stick Fixed Stability and Drag Otimization Simulation Time: ' + str(elapsed_time_stick_fixed))    
    
    '''
    ELEVATOR SIZING
    '''      
    # define vehicle for elevator sizing   
    optimized_vehicle_v1                             = planform_optimization_problem.vehicle_configurations.stick_fixed_cruise
    optimized_vehicle_v1.maxiumum_load_factor        = 3.0
    optimized_vehicle_v1.minimum_load_factor         = -1
    
    ti = time.time()   
    solver_name       = 'SLSQP'  
    elevator_sizing_optimization_problem = elevator_sizing_optimization_setup(optimized_vehicle_v1,cruise_velocity,cruise_altitude)
    output_elevator_sizing = scipy_setup.SciPy_Solve(elevator_sizing_optimization_problem,solver=solver_name, sense_step = 1E-2, tolerance = 1E-2) 
    print (output_elevator_sizing)     
    tf           = time.time()
    elapsed_time_elevator_sizing = round((tf-ti)/60,2)
    print('Elevator Sizing Simulation Time: ' + str(elapsed_time_elevator_sizing))   
     
    '''
    AILERON AND RUDDER SIZING
    '''      
    # define vehicle for aileron and rudder sizing
    optimized_vehicle    = elevator_sizing_optimization_problem.vehicle_configurations.elevator_sizing 
    optimized_vehicle.rudder_flag                       = True  
    optimized_vehicle.crosswind_velocity                = 20 * Units.knots
    
    ti = time.time()   
    solver_name       = 'SLSQP'  
    aileron_rudder_sizing_optimization_problem = aileron_rudder_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude)
    output_aileron_and_rudder_sizing = scipy_setup.SciPy_Solve(aileron_rudder_sizing_optimization_problem,solver=solver_name, sense_step = 1E-1, tolerance = 1E-1) 
    print (output_aileron_and_rudder_sizing)     
    tf           = time.time()
    elapsed_time_aileron_and_rudder_sizing = round((tf-ti)/60,2)
    print('Aileron and Rudder Sizing Simulation Time: ' + str(elapsed_time_aileron_and_rudder_sizing))   
    
    '''
    FLAP SIZING
    '''      
    # define vehicle for flap sizing     
    optimized_vehicle_v3 = aileron_rudder_sizing_optimization_problem.vehicle_configurations.aileron_rudder_sizing 
    
    ti = time.time()   
    solver_name       = 'SLSQP'  
    flap_sizing_optimization_problem = flap_sizing_optimization_setup(optimized_vehicle_v3,cruise_velocity,cruise_altitude)
    output_flap_sizing = scipy_setup.SciPy_Solve(flap_sizing_optimization_problem,solver=solver_name, sense_step = 1E-4, tolerance = 1E-4) 
    print (output_flap_sizing)     
    tf           = time.time()
    elapsed_time_flap_sizing = round((tf-ti)/60,2)
    print('Flap Sizing Simulation Time: ' + str(elapsed_time_flap_sizing))   
    
    '''
    PRINT VEHICLE CONTROL SURFACES
    '''          
    optimized_vehicle_v4  = flap_sizing_optimization_problem.vehicle_configurations.flap_sizing 
    print_vehicle_control_surface_geoemtry(optimized_vehicle_v4)
     
    tf_0           = time.time()
    total_elapsed_time = round((tf_0-ti_0)/60,2)    
    print('Total Control Surface Sizing Time: ' + str(total_elapsed_time))
    
    Optimization_Data_2_CSV(configuration_CGbatt_MOIbatt,
                            optimized_vehicle_v4.tag, 
                            output_stick_fixed, 
                            output_elevator_sizing, 
                            output_aileron_and_rudder_sizing, 
                            output_flap_sizing, 
                            planform_optimization_problem, 
                            elevator_sizing_optimization_problem, 
                            aileron_rudder_sizing_optimization_problem, 
                            flap_sizing_optimization_problem,
                            elapsed_time_stick_fixed,
                            elapsed_time_elevator_sizing,
                            elapsed_time_aileron_and_rudder_sizing,
                            elapsed_time_flap_sizing)
    
    # Save Vehicle!!!
    
    return
  
def stick_fixed_stability_and_drag_optimization_setup(vehicle,cruise_velocity,cruise_altitude): 
    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    
    
    scaling_factor = 0.05

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #             [ tag,                          initial,                                       (lb , ub) , scaling , units ]  
    problem.inputs = np.array([       
                  #[ 'mw_span'                     , 11.82855  , 10  , 13   , 1.0  ,  1*Units.less],    
                  #[ 'mw_AR'                       , 8.95198 , 7 , 10  , 10. ,  1*Units.meter**2],                                                                                                         
                  [ 'mw_root_twist'               , vehicle.wings.main_wing.twists.root,          vehicle.wings.main_wing.twists.root*(1 - scaling_factor),          vehicle.wings.main_wing.twists.root*(1 + scaling_factor),           1.,   1*Units.degree], 
                  [ 'mw_tip_twist'                , vehicle.wings.main_wing.twists.tip,           vehicle.wings.main_wing.twists.tip*(1 - scaling_factor),           vehicle.wings.main_wing.twists.tip*(1 + scaling_factor),            1.,   1*Units.degree], 
                  [ 'vt_span'                     , vehicle.wings.vertical_tail.spans.projected,  vehicle.wings.vertical_tail.spans.projected*(1 - scaling_factor),  vehicle.wings.vertical_tail.spans.projected*(1 + scaling_factor),   1.,   1*Units.meter],  
                  [ 'vt_AR'                       , vehicle.wings.vertical_tail.aspect_ratio,     vehicle.wings.vertical_tail.aspect_ratio*(1 - scaling_factor),     vehicle.wings.vertical_tail.aspect_ratio*(1 + scaling_factor),      100., 1*Units.meter**2],    
                  [ 'ht_span'                     , vehicle.wings.horizontal_tail.spans.projected,vehicle.wings.horizontal_tail.spans.projected*(1 - scaling_factor),vehicle.wings.horizontal_tail.spans.projected*(1 + scaling_factor), 10.,  1*Units.less], 
                  [ 'ht_AR'                       , vehicle.wings.horizontal_tail.aspect_ratio,   vehicle.wings.horizontal_tail.aspect_ratio*(1 - scaling_factor),   vehicle.wings.horizontal_tail.aspect_ratio*(1 + scaling_factor),    10.,  1*Units.meter**2], 
                  [ 'AoA'                         , 5     , -10  , 10    , 1    ,  1*Units.degree],  
                  
    ],dtype=object)   

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                 [  'CD'  ,  1.0  ,    1*Units.less] 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ] See http://everyspec.com/MIL-SPECS/MIL-SPECS-MIL-F/MIL-F-8785C_5295/
    # Level 1, Category B
    problem.constraints = np.array([
        [ 'CL_residual'               ,   '<' ,   1E-2  ,   1E-2  , 1*Units.less], # close to zero 2 works 
        [ 'CM_residual'               ,   '<' ,   1E-2  ,   1E-2  , 1*Units.less], # close to zero 2 works 
        [ 'static_margin'             ,   '>' ,   0.1   ,   0.1   , 1*Units.less],  # checked 
        [ 'CM_alpha'                  ,   '<' ,   0.0   ,   1.0   , 1*Units.less],  # checked 
        [ 'phugoid_damping_ratio'     ,   '>' ,   0.04  ,   1.0   , 1*Units.less],  # checked 
        [ 'short_period_damping_ratio',   '<' ,   2.0   ,   1.0   , 1*Units.less],  # checked 
        [ 'short_period_damping_ratio',   '>' ,   0.3   ,   1.0   , 1*Units.less], # checked    
        [ 'dutch_roll_frequency'      ,   '>' ,   0.4   ,   1.0   , 1*Units.less],  # checked   frequency in rad/sec
        [ 'dutch_roll_damping_ratio'  ,   '>' ,   0.08  ,   1.0   , 1*Units.less],  # checked   
        [ 'spiral_doubling_time'      ,   '>' ,   20.0  ,   1.0   , 1*Units.less],  # checked   
        [ 'spiral_criteria'           ,   '>' ,   1.0   ,   1.0   , 1*Units.less],  # checked   
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'CL_residual'                       , 'summary.CL_stick_fixed_residual' ],
        [ 'CD'                                , 'summary.CD' ],
        [ 'CM_residual'                       , 'summary.CM_residual' ],  
        [ 'CM_alpha'                          , 'summary.CM_alpha' ],    
        [ 'static_margin'                     , 'summary.static_margin' ], 
        [ 'phugoid_damping_ratio'             , 'summary.phugoid_damping_ratio' ],  
        [ 'short_period_damping_ratio'        , 'summary.short_period_damping_ratio' ],  
        [ 'dutch_roll_frequency'              , 'summary.dutch_roll_frequency' ],  
        [ 'dutch_roll_damping_ratio'          , 'summary.dutch_roll_damping_ratio' ],  
        [ 'spiral_doubling_time'              , 'summary.spiral_doubling_time' ],   
        [ 'spiral_criteria'                   , 'summary.spiral_criteria' ],       
        #[ 'mw_span'                           , 'vehicle_configurations.*.wings.main_wing.spans.projected'],
        #[ 'mw_AR'                             , 'vehicle_configurations.*.wings.main_wing.aspect_ratio'],         
        [ 'mw_root_twist'                     , 'vehicle_configurations.*.wings.main_wing.twists.root' ], 
        [ 'mw_tip_twist'                      , 'vehicle_configurations.*.wings.main_wing.twists.tip'  ],     
        [ 'vt_span'                           , 'vehicle_configurations.*.wings.vertical_tail.spans.projected'],
        [ 'vt_AR'                             , 'vehicle_configurations.*.wings.vertical_tail.aspect_ratio'],         
        [ 'ht_AR'                             , 'vehicle_configurations.*.wings.horizontal_tail.aspect_ratio'],     
        [ 'ht_span'                           , 'vehicle_configurations.*.wings.horizontal_tail.spans.projected'], 
        [ 'AoA'                               , 'missions.stick_fixed_cruise.segments.cruise.state.conditions.aerodynamics.angles.alpha[0][0]' ], 
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.stick_fixed_stability_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.stick_fixed_stability_setup(nexus.analyses,nexus.vehicle_configurations.stick_fixed_cruise,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.stick_fixed_stability_and_drag_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus 

def elevator_sizing_optimization_setup(vehicle,cruise_velocity,cruise_altitude):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    problem.inputs = np.array([            
                  [ 'ht_elevator_chord_fraction' , 0.35    , 0.15 , 0.45   ,  1.0  ,  1*Units.less],
                  [ 'ht_elevator_span_frac_start', 0.1   , 0.05 , 0.5   ,  1.0 ,  1*Units.less], 
                  [ 'ht_elevator_span_frac_end'  , 0.9   , 0.6  , 0.95  ,  1.0  ,  1*Units.less],     
                  
    ],dtype=object)   

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'elevator_surface_area', 1. , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([ 
        [ 'elevator_push_deflection'           ,   '<' ,  30 , 1.0   , 1*Units.less], 
        [ 'elevator_pull_deflection'           ,   '<' ,  30 , 1.0    , 1*Units.less], 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'elevator_surface_area'             , 'summary.elevator_surface_area' ], 
        [ 'elevator_push_deflection'          , 'summary.elevator_push_deflection' ],   
        [ 'elevator_pull_deflection'          , 'summary.elevator_pull_deflection' ],     
        [ 'ht_elevator_chord_fraction'        , 'vehicle_configurations.*.wings.horizontal_tail.control_surfaces.elevator.chord_fraction'],    
        [ 'ht_elevator_span_frac_start'       , 'vehicle_configurations.*.wings.horizontal_tail.control_surfaces.elevator.span_fraction_start'],    
        [ 'ht_elevator_span_frac_end'         , 'vehicle_configurations.*.wings.horizontal_tail.control_surfaces.elevator.span_fraction_end'],  
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.elevator_sizing_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.elevator_sizing_setup(nexus.analyses,nexus.vehicle_configurations.elevator_sizing,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.elevator_sizing_setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  


 
 
def aileron_rudder_sizing_optimization_setup(vehicle,cruise_velocity,cruise_altitude):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    if vehicle.rudder_flag:
        problem.inputs = np.array([             
                      [ 'mw_aileron_chord_fraction'  , 0.2    , 0.15 , 0.3  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_start' , 0.7   , 0.55 , 0.8  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_end'   , 0.9    , 0.85 , 0.95 ,  1.0 ,  1*Units.less],  
                      [ 'vs_rudder_chord_fraction'   , 0.4    , 0.15 , 0.5  ,  1.0 ,  1*Units.less],
                      [ 'vs_rudder_span_frac_start'  , 0.1   , 0.05 , 0.35 ,  1.0 ,  1*Units.less],
                      [ 'vs_rudder_span_frac_end'    , 0.9   , 0.5  , 0.95 ,  1.0 ,  1*Units.less] ],dtype=object)   
    else:
        problem.inputs = np.array([             
                      [ 'mw_aileron_chord_fraction'  , 0.2    , 0.15 , 0.3  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_start' , 0.7   , 0.55 , 0.8  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_end'   , 0.9    , 0.85 , 0.95 ,  1.0 ,  1*Units.less],   
                      
        ],dtype=object)      

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'aileron_rudder_surface_area', 1. , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    if vehicle.rudder_flag:
        problem.constraints = np.array([
            [ 'aileron_roll_deflection'                ,   '<' ,   30 , 1.0     , 1*Units.less], 
            [ 'rudder_roll_deflection'                 ,   '<' ,   30 , 1.0     , 1*Units.less],  
            [ 'aileron_crosswind_deflection'           ,   '<' ,   30 , 1.0     , 1*Units.less], 
            [ 'rudder_crosswind_deflection'            ,   '<' ,   30 , 1.0     , 1*Units.less],  
        ],dtype=object)
    else:
        problem.constraints = np.array([
            [ 'aileron_roll_deflection'       ,   '<' ,   30 , 1.0     , 1*Units.less],  
            [ 'aileron_crosswind_deflection'  ,   '<' ,   30 , 1.0     , 1*Units.less],  
        ],dtype=object)
        
        
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    if vehicle.rudder_flag:
        problem.aliases = [ 
            [ 'aileron_rudder_surface_area'            , 'summary.aileron_rudder_surface_area' ],  
            [ 'aileron_roll_deflection'                , 'summary.aileron_roll_deflection' ],  
            [ 'rudder_roll_deflection'                 , 'summary.rudder_roll_deflection' ], 
            [ 'aileron_crosswind_deflection'           , 'summary.aileron_crosswind_deflection' ],      
            [ 'rudder_crosswind_deflection'            , 'summary.rudder_crosswind_deflection' ],         
            [ 'mw_aileron_chord_fraction'              , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
            [ 'mw_aileron_span_frac_start'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
            [ 'mw_aileron_span_frac_end'               , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end'],     
            [ 'vs_rudder_chord_fraction'               , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.chord_fraction'],    
            [ 'vs_rudder_span_frac_start'              , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.span_fraction_start'],    
            [ 'vs_rudder_span_frac_end'                , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.span_fraction_end']]      
    else:
        problem.aliases = [ 
            [ 'aileron_rudder_surface_area'            , 'summary.aileron_rudder_surface_area' ],  
            [ 'aileron_roll_deflection_residual'       , 'summary.aileron_roll_deflection_residual' ],          
            [ 'aileron_crosswind_deflection_residual'  , 'summary.aileron_crosswind_deflection_residual' ],   
            [ 'mw_aileron_chord_fraction'              , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
            [ 'mw_aileron_span_frac_start'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
            [ 'mw_aileron_span_frac_end'               , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end']] 
        
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.aileron_rudder_sizing_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.aileron_rudder_sizing_setup(nexus.analyses,nexus.vehicle_configurations.aileron_rudder_sizing,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.aileron_rudder_sizing_setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  


def flap_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    problem.inputs = np.array([           
                  [ 'mw_flap_chord_fraction'     , 0.2    , 0.15 , 0.4  ,  1.0 ,  1*Units.less],
                  [ 'mw_flap_span_frac_start'    , 0.2    , 0.05 , 0.25 ,  1.0 ,  1*Units.less],
                  [ 'mw_flap_span_frac_end'      , 0.5    , 0.3  , 0.6  ,  1.0 ,  1*Units.less],   
                  
    ],dtype=object)   

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'flap_surface_area', 1. , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([
        [ 'flap_criteria'    ,   '>' ,  0.   ,  1.0   , 1*Units.less],  
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'flap_surface_area'                 , 'summary.flap_surface_area' ], 
        [ 'flap_criteria'                     , 'summary.flap_criteria' ],   
        [ 'mw_flap_chord_fraction'            , 'vehicle_configurations.*.wings.main_wing.control_surfaces.flap.chord_fraction'],    
        [ 'mw_flap_span_frac_start'           , 'vehicle_configurations.*.wings.main_wing.control_surfaces.flap.span_fraction_start'],    
        [ 'mw_flap_span_frac_end'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.flap.span_fraction_end'],    
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.flap_sizing_setup(optimized_vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.flap_sizing_setup(nexus.analyses,nexus.vehicle_configurations.flap_sizing,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.flap_sizing_setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  
    
def print_vehicle_control_surface_geoemtry(vehicle): 
   
    for wing in vehicle.wings:
        if 'control_surfaces' in wing:  
            for CS in wing.control_surfaces:  
                print('Wing                : ' + wing.tag)
                print('Control Surface     : ' + CS.tag)
                print('Span Fraction Start : ' + str(CS.span_fraction_start))
                print('Span Fraction End   : ' + str(CS.span_fraction_end)) 
                print('Chord Fraction      : ' + str(CS.chord_fraction)) 
                print("\n\n")     

    return

if __name__ == '__main__':
    main()
    plt.show()
